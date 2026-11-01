"""Datasheet-driven q_0 / tau extraction for the DiodeRecov body-diode model.

Protocol (see out/RECOV_EXTRACTION_PROTOCOL.md):
  1. Build a FAST single-event IEC-style recovery bench (seeded .ic states)
     matching the datasheet test conditions (I_F, V_R, di/dt). The netlist is
     a scratch file: generated per run, deleted when the run is parsed.
     dTmax = t_rr_target / 10 (snap-region step control; see protocol §2).
  2. One LTspice -b process with nested .step over a coarse (q_0 x tau) grid;
     .meas directives return t_rr (10% and 25% of I_RM) and Qrr.
  3. Fit bilinear surrogates Qrr(q_0,tau), t_rr(q_0,tau); solve the 2x2 inverse
     problem (exact if the target curves cross, else constrained least squares).
  4. Three variants are reported and plotted (recovery current waveforms):
       A. LS optimum        - balanced (Qrr, t_rr) least squares
       B. Qrr-strict        - Qrr pinned exactly, tau = 1 (timing-strict)
       C. shape-strict      - Qrr pinned exactly, tau = 5 (soft-tail shape)
  5. One check run at the low-di/dt datasheet point; residuals reported as the
     model-family limitation.

Total cost: ~10 fast runs, ~0.1 s each -> seconds, not minutes.
"""
from __future__ import annotations

import glob
import os
import re
import subprocess
import sys

import numpy as np
from scipy.optimize import minimize

# ---------------------------------------------------------------------------
# bench generation
# ---------------------------------------------------------------------------

NETLIST_TMPL = """* extract_recov bench — {label}
.param i_test={i_f} didt={didt} tr=i_test*2/didt*1e-6 tf=tr
.param k_chg=5 r_max=0
Vv N001 0 {v_r}
D5 SW N001 IdealD
C3 SW N001 1p
I1 0 SW PULSE({{-i_test}} {{i_test}} {td} {{tr}} {{tf}} 10m 20m)
X1 0 N002 DiodeRecov params: Is=1n Rs=5m N=1.4
R4 N002 SW 1n
.ic V(X1:q_j)={qj_seed} V(X1:Ki)=-0.9

.subckt DiodeRecov A K
D1 Ai Ki MyDiode
B2 Ki Ai I=max(V(q_j)*(V(Ki)-V(Ai))/tau,0)
C1 q_j 0 {{q_0/tau}} Rser=0
B1 0 q_j I=-I(R2)/tau
D2 0 q_j ClampD
C2 Ki A {cds}
R1 Ki N001 {{Rs}}
R2 Ai A 1n
B3 N001 K V=-I(R2)*(r_max/(1+k_chg*V(q_j))+1u)
.model ClampD D(Ron=1u,Roff={roff},Vfwd=0)
.model MyDiode D(Is=Is,Rs=0,N=N)
.ends DiodeRecov

.model D D
.lib {diolib}
.tran 0 {tstop} 0 {dtmax}
.model IdealSW SW(Ron=1m,Roff=10Meg,vt=5,vh=-2)
.model IdealD D(Ron=1m,Roff=1Meg)
.meas i_pk MIN I(R4) FROM {td_meas}
.meas t_rr10 TRIG I(R4)=0 TD={td_meas} FALL=1 TARG I(R4)={{i_pk*0.1}} TD={td_meas} RISE=1
.meas t_rr25 TRIG I(R4)=0 TD={td_meas} FALL=1 TARG I(R4)={{i_pk*0.25}} TD={td_meas} RISE=1
.meas q_rr INTEG -I(R4)*1e9 TRIG I(R4)=0 TD={td_meas} FALL=1 TARG I(R4)={{i_pk*0.1}} TD={td_meas} RISE=1
.meas v_rail MAX V(SW) FROM {td_meas}
.options reltol=1e-4 method=gear plotwinsize=0 solver=alt
{steps}.end
"""


def find_ltspice() -> str:
    cands = []
    for base in (r"C:\Program Files\ADI\LTspice", r"C:\Program Files (x86)\ADI\LTspice",
                 os.path.expandvars(r"%LOCALAPPDATA%\LTspice")):
        cands += glob.glob(os.path.join(base, "LTspice.exe"))
    if not cands:
        raise RuntimeError("LTspice.exe not found — set the path in find_ltspice()")
    return cands[0]


LTSPICE = find_ltspice()
DIOLIB = os.path.expandvars(r"%LOCALAPPDATA%\LTspice\lib\cmp\standard.dio")

# measured C_ds charge-master ramp-sum (user bench edit; from
# out/ISC040N10NM8_caps.subckt via cds_expr.txt): "Cds D S Q=<expr>"
CDS_EXPR_FILE = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))),
    "data", "ISC040N10NM8_cds_expr.txt")


def _cds_clause():
    try:
        with open(CDS_EXPR_FILE) as f:
            return "Q=" + f.read().strip().split("Q=", 1)[1]
    except (OSError, IndexError):
        return "500p"  # placeholder fallback


def make_netlist(i_f, v_r, didt, q_list, tau_list, dtmax, label="bench",
                 tm_list=None):
    """Fast single-event bench netlist. V_R is a stiff source (Vv N001 0 v_r)
    with D5 SW->N001 as the clamp: the branch reverse voltage is exactly V_R
    (was a self-bias rail ~I_F*R2s/2 with C2/R2s; the source removes the
    ~1.3% rail error and the N001 .ic seed entirely).

    tm_list=None: two-knob model (ClampD Roff=tau, V_q,ss=I_F).
    tm_list given: three-knob model, Roff=tau*t_minor/q_0 -> storage time
    constant = t_minor, V_q,ss = I_F*t_minor/q_0 (identity with two-knob at
    t_minor=q_0; verified by probe)."""
    three = tm_list is not None
    roff = "tau*t_minor/q_0" if three else "tau"
    qj_seed = "{i_test*t_minor/q_0}" if three else "{i_test}"
    td = 0.5e-6
    tstop = td + 2 * (2 * i_f / didt * 1e-6) + 1.1e-6
    steps = ""
    def emit(name, vals, fmt="{:g}n"):
        nonlocal steps
        if len(vals) == 1:  # single-value .step is rejected by LTspice
            steps += f".param {name}={fmt.format(vals[0])}\n"
        else:
            steps += (f".step param {name} list " +
                      " ".join(fmt.format(v) for v in vals) + "\n")
    emit("q_0", q_list)
    emit("tau", tau_list, "{:g}")
    if three:
        emit("t_minor", tm_list)
    return NETLIST_TMPL.format(
        label=label, i_f=i_f, v_r=v_r, didt=didt, td=td, tstop=tstop,
        dtmax=f"{dtmax:g}", td_meas=td * 0.8, cds=_cds_clause(),
        diolib=DIOLIB, steps=steps, roff=roff, qj_seed=qj_seed)


# ---------------------------------------------------------------------------
# LTspice run + scratch cleanup
# ---------------------------------------------------------------------------

def run_ltspice(netlist, workdir, want_waveform=False):
    """Run LTspice -b on a scratch netlist; parse the log (and the ASCII raw
    when want_waveform), then DELETE the netlist and all run artifacts."""
    os.makedirs(workdir, exist_ok=True)
    path = os.path.join(workdir, "bench.net")
    with open(path, "w") as f:
        f.write(netlist)
    cmd = [LTSPICE, "-b"] + (["-ascii"] if want_waveform else []) + [path]
    subprocess.run(cmd, cwd=workdir, timeout=300, check=False)
    log = path[:-4] + ".log"
    with open(log, "r", errors="ignore") as f:
        text = f.read()
    if "Total elapsed time" not in text:
        raise RuntimeError("LTspice run failed:\n" + text[-2000:])
    wave = parse_ascii_raw(path[:-4] + ".raw", ["time", "I(R4)"]) if want_waveform else None
    # scratch cleanup: netlist + raw/log/db artifacts
    for leftover in glob.glob(path[:-4] + ".*"):
        os.remove(leftover)
    return text, wave


def parse_ascii_raw(raw_path, traces):
    """Parse an LTspice ASCII .raw (from `ltspice -b -ascii`). Returns dict
    name -> np.array, keys lowercased as given."""
    with open(raw_path, "r", errors="ignore") as f:
        lines = f.read().splitlines()
    nv = int(next(l for l in lines if l.startswith("No. Variables")).split(":")[1])
    vi = next(i for i, l in enumerate(lines) if l.startswith("Variables:"))
    names = [lines[i].split()[1].strip('"') for i in range(vi + 1, vi + 1 + nv)]
    vals_i = next(i for i, l in enumerate(lines) if l.strip() == "Values:")
    # LTspice ASCII raw: one point per group — index line, then one tab-
    # prefixed value per line. Tokenize the whole block and walk it.
    tokens = [t for l in lines[vals_i + 1:] for t in l.split()]
    rows = []
    i = 0
    while i + nv < len(tokens) + 1:
        if not re.fullmatch(r"\d+", tokens[i]):
            break
        rows.append([float(x) for x in tokens[i + 1:i + 1 + nv]])
        i += 1 + nv
    arr = np.array(rows)
    # value column k corresponds to variable k (the point index column is dropped)
    want = {t.lower(): arr[:, names.index(t)] for t in traces}
    return want


def parse_stepped(text):
    """Return {meas_name: [(step, value), ...]} from stepped .measure tables."""
    out = {}
    for m in re.finditer(r"Measurement: (\w+)\r?\n[^\r\n]*\r?\n((?:\s+\d+\s+[^\r\n]+\r?\n)+)", text):
        name, rows = m.group(1), m.group(2)
        vals = []
        for line in rows.strip().splitlines():
            parts = line.split()
            try:
                vals.append((int(parts[0]), float(parts[1])))
            except ValueError:  # failed .meas rows print 'failed', not a number
                continue
        out[name] = vals
    return out


# ---------------------------------------------------------------------------
# surrogate + inverse solve
# ---------------------------------------------------------------------------

def bilinear(q, tau, z):
    A = np.column_stack([np.ones(len(q)), q, tau, q * tau])
    c, *_ = np.linalg.lstsq(A, z, rcond=None)
    return c, np.abs(A @ c - z).max()


def surrogate(c, q, tau):
    return c[0] + c[1] * q + c[2] * tau + c[3] * q * tau


def _meas_value(text, name):
    m = re.search(rf'^{name}[:=]\s*(?:[A-Za-z_()0-9*,-]*=)?([\d.eE+-]+)', text, re.M)
    if not m:
        m = re.search(rf'Measurement: {name}[\s\S]*?\s+1\s+([\d.eE+-]+)', text)
    return float(m.group(1)) if m else None


# ---------------------------------------------------------------------------
# three-knob (q_0, tau, t_minor) surrogates
# ---------------------------------------------------------------------------

def grid_points(q_list, tau_list, tm_list):
    """(q_0, tau, t_minor) per LTspice step index. Nested .step directives
    nest with the LAST directive as the outer loop (verified): q_0 fastest,
    tau middle, t_minor outermost."""
    nq, nt, nm = len(q_list), len(tau_list), len(tm_list)
    return [(q_list[(s - 1) % nq], tau_list[(s - 1) // nq % nt],
             tm_list[(s - 1) // (nq * nt)])
            for s in range(1, nq * nt * nm + 1)]


def trilinear_fit(Q, T, M, z):
    keep = ~np.isnan(z)
    A = np.column_stack([np.ones(keep.sum()), Q[keep], T[keep], M[keep],
                         Q[keep] * T[keep], Q[keep] * M[keep],
                         T[keep] * M[keep], Q[keep] * T[keep] * M[keep]])
    c, *_ = np.linalg.lstsq(A, z[keep], rcond=None)
    return c, np.abs(A @ c - z[keep]).max()


def trilinear_eval(c, q, t, m):
    return (c[0] + c[1] * q + c[2] * t + c[3] * m + c[4] * q * t +
            c[5] * q * m + c[6] * t * m + c[7] * q * t * m)


# ---------------------------------------------------------------------------
# plotting
# ---------------------------------------------------------------------------

def plot_waveforms(waves, ttr_ns, plot_path):
    """waves: list of (label, q_0_nC, tau, {'time':..,'i(r4)':..}).
    Zoomed on the recovery event: zero crossing -> tail end."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (ax1, ax2) = plt.subplots(
        2, 1, sharex=True, figsize=(9, 6.5),
        gridspec_kw={"height_ratios": [3, 1]})
    xlo, xhi = -15.0, max(4.0 * ttr_ns, 80.0)   # recovery completes in ~2-3 t_rr
    tail_min = 0.0
    for label, q_nC, tau, w in waves:
        t = w["time"] * 1e9
        i = w["i(r4)"]
        # downward zero crossing after the pre-event plateau
        idx = np.flatnonzero((i[:-1] > 0) & (i[1:] <= 0))
        if len(idx) == 0:
            raise RuntimeError(f"no recovery zero crossing for {label}")
        t0 = t[idx[0]]
        m = (t >= t0 + xlo) & (t <= t0 + xhi)
        ax1.plot(t - t0, i, label=f"{label}  (q_0={q_nC:.1f} n, τ={tau:g})")
        ax2.plot(t - t0, i, label=None)
        tail = i[m & (i < -1e-3) & (i > -1.0)]   # tail region only (post-snap)
        if len(tail):
            tail_min = min(tail_min, tail.min())
    ax1.axhline(0, color="k", lw=0.5)
    ax1.set_ylabel("branch current I(R4) [A]")
    ax1.set_title("Body-diode reverse recovery — extraction variants "
                  "(time relative to I(R4) zero crossing)")
    ax1.legend(loc="lower left", fontsize=8)
    ax1.grid(alpha=0.3)
    ax2.set_xlim(xlo, xhi)
    ax2.set_ylim(1.3 * tail_min if tail_min < 0 else -0.1, 0.02)
    ax2.axhline(0, color="k", lw=0.5)
    ax2.set_xlabel("t − t0 [ns]")
    ax2.set_ylabel("tail [A]")
    ax2.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(plot_path, dpi=150)
    plt.close(fig)
    return plot_path


# ---------------------------------------------------------------------------
# top-level extraction
# ---------------------------------------------------------------------------

def extract(i_f, v_r, didt_hi, qrr_hi_nC, ttr_hi_ns, tdef="25%",
            didt_lo=None, qrr_lo_nC=None, ttr_lo_ns=None,
            tau_grid=(0.5, 1, 2, 5), workdir=None, plot_path=None):
    """Returns (variants, report_text). variants: {label: (q_0_nC, tau)}."""
    workdir = workdir or os.path.join(os.getcwd(), "recov_extract_run")
    dtmax = ttr_hi_ns * 1e-9 / 10.0   # step control: 1/10 of the target t_rr
    q_guess = (qrr_hi_nC + 35.0) / 8.2   # empirical slope ~8.2 nC/nC
    q_grid = (0.85 * q_guess, 1.15 * q_guess)

    # 1) grid run
    text, _ = run_ltspice(make_netlist(i_f, v_r, didt_hi, q_grid, tau_grid,
                                       dtmax, "grid"), workdir)
    tab = parse_stepped(text)
    pts = []
    nq = len(q_grid)
    for s in range(1, nq * len(tau_grid) + 1):
        # LTspice nests .step directives with the LAST one as the outer loop
        pts.append((q_grid[(s - 1) % nq], tau_grid[(s - 1) // nq]))
    Q = np.array([p[0] for p in pts])
    T = np.array([p[1] for p in pts])
    key_t = "t_rr25" if tdef == "25%" else "t_rr10"
    Qrr = np.array([d for _, d in sorted(tab["q_rr"])])        # nC (meas x1e9)
    Ttr = np.array([d for _, d in sorted(tab[key_t])]) * 1e9   # s -> ns
    cQ, rQ = bilinear(Q, T, Qrr)
    cT, rT = bilinear(Q, T, Ttr)
    print(f"grid q_0=({q_grid[0]:.1f},{q_grid[1]:.1f}) n  tau={tau_grid}  "
          f"dTmax={dtmax*1e9:.1f} ns  surrogate max resid: "
          f"Qrr {rQ:.2f} nC, t_{tdef} {rT:.3f} ns")

    def fQ(q, t): return surrogate(cQ, q, t)
    def fT(q, t): return surrogate(cT, q, t)

    # 2) solve 2x2 inverse problem
    taus = np.linspace(0.2, 10, 200)
    qQ = (qrr_hi_nC - cQ[0] - cQ[2] * taus) / (cQ[1] + cQ[3] * taus)
    qT = (ttr_hi_ns - cT[0] - cT[2] * taus) / (cT[1] + cT[3] * taus)
    gap = qQ - qT
    ok = (taus >= 0.2) & (qQ > 0) & (qT > 0)
    cross = np.flatnonzero(np.diff(np.sign(gap[ok])))
    if len(cross):
        i0 = np.flatnonzero(ok)[cross[0]]
        q_star, t_star = qQ[i0], taus[i0]
        mode = "exact 2x2 solve"
    else:
        def cost(x):
            return ((fQ(x[0], x[1]) - qrr_hi_nC) / (0.1 * qrr_hi_nC)) ** 2 + \
                   ((fT(x[0], x[1]) - ttr_hi_ns) / (0.1 * ttr_hi_ns)) ** 2
        r = minimize(cost, [q_guess, 1.0], bounds=[(5, 100), (0.5, 10)],
                     method="L-BFGS-B")
        q_star, t_star = r.x
        mode = "constrained least squares (no exact solution — model-family limit)"

    # 3) three variants (B/C pin Qrr exactly at fixed tau, re-pinned via surrogate)
    variants = {
        "A: LS optimum": (float(q_star), float(t_star)),
        "B: Qrr-strict (tau=1)": (float((qrr_hi_nC - cQ[0] - cQ[2] * 1.0)
                                        / (cQ[1] + cQ[3] * 1.0)), 1.0),
        "C: shape-strict (tau=5)": (float((qrr_hi_nC - cQ[0] - cQ[2] * 5.0)
                                          / (cQ[1] + cQ[3] * 5.0)), 5.0),
    }

    # 4) waveform runs (ASCII raw parsed, scratch deleted) + measurements
    lines = [
        f"== q_0/tau extraction ({mode}) ==",
        f"targets @ {didt_hi} A/us: Qrr={qrr_hi_nC} nC, t_{tdef}={ttr_hi_ns} ns "
        f"(I_F={i_f} A, V_R={v_r} V), dTmax={dtmax*1e9:.1f} ns",
    ]
    waves = []
    for label, (q_nC, tau) in variants.items():
        text, wave = run_ltspice(make_netlist(i_f, v_r, didt_hi, [q_nC], [tau],
                                              dtmax, label), workdir,
                                 want_waveform=True)
        Qrr_v = _meas_value(text, "q_rr")
        Ttr_v = _meas_value(text, key_t) * 1e9
        lines.append(
            f"{label}: q_0={q_nC:.2f} n, tau={tau:.3f} -> "
            f"Qrr={Qrr_v:.1f} nC ({(Qrr_v-qrr_hi_nC)/qrr_hi_nC*100:+.1f}%), "
            f"t_{tdef}={Ttr_v:.2f} ns ({(Ttr_v-ttr_hi_ns)/ttr_hi_ns*100:+.1f}%)")
        waves.append((label, q_nC, tau, wave))

    # 5) low-di/dt residual check (model-family limitation, reported not fitted)
    if didt_lo:
        text, _ = run_ltspice(make_netlist(i_f, v_r, didt_lo, [q_star], [t_star],
                                           dtmax, "didt_lo"), workdir)
        Qrr_l = _meas_value(text, "q_rr")
        Ttr_l = _meas_value(text, key_t) * 1e9
        lines.append(
            f"low-di/dt check @ {didt_lo} A/us (variant A): Qrr={Qrr_l:.1f} nC "
            f"(datasheet {qrr_lo_nC}), t_{tdef}={Ttr_l:.2f} ns "
            f"(datasheet {ttr_lo_ns}) — residual is the two-knob stored-charge "
            f"limitation, do not chase with extra knobs unless the waveform "
            f"demands it")

    report = "\n".join(lines)
    print(report)
    os.makedirs(workdir, exist_ok=True)
    with open(os.path.join(workdir, "extract_report.txt"), "w") as f:
        f.write(report + f"\n\nworkdir: {workdir}\n")

    if plot_path:
        plot_waveforms(waves, ttr_hi_ns, plot_path)
        print(f"waveform plot: {plot_path}")
    return variants, report


# ---------------------------------------------------------------------------
# three-knob extraction: fit (q_0, tau, t_minor) to BOTH di/dt points
# ---------------------------------------------------------------------------

def extract3k(i_f, v_r, didt_hi, qrr_hi_nC, ttr_hi_ns,
              didt_lo, qrr_lo_nC, ttr_lo_ns, tdef="25%",
              workdir=None, plot_path=None, weights=(1.0, 1.0, 1.0, 1.0)):
    """Fit q_0, tau, t_minor (ClampD Roff = tau*t_minor/q_0) to four targets:
    (Qrr, t_rr) at the high-di/dt AND the low-di/dt datasheet points.
    weights: per-term multipliers (Qrr_hi, t_hi, Qrr_lo, t_lo) on the squared
    normalized residuals; e.g. (1000, 10, 1, 1) pins the hard-switching point.
    Trust-region: local 2x2x2 grid (exact trilinear surrogate) around the
    current best -> solve inside the box -> accept on REAL cost only.
    Returns ((q_nC, tau, tm_ns), report)."""
    workdir = workdir or os.path.join(os.getcwd(), "recov_extract3k_run")
    dtmax = ttr_hi_ns * 1e-9 / 10.0
    key_t = "t_rr25" if tdef == "25%" else "t_rr10"

    def grid(qg, taug, tmg):
        pts = grid_points(qg, taug, tmg)
        n = len(pts)
        Q = np.array([p[0] for p in pts])
        T = np.array([p[1] for p in pts])
        M = np.array([p[2] for p in pts])
        resp = {}
        for tag, didt in (("hi", didt_hi), ("lo", didt_lo)):
            text, _ = run_ltspice(make_netlist(i_f, v_r, didt, qg, taug,
                                               dtmax, f"grid_{tag}",
                                               tm_list=tmg), workdir)
            tab = parse_stepped(text)
            for meas, out, scale in (("q_rr", "Qrr", 1.0),
                                     (key_t, "Ttr", 1e9),
                                     ("i_pk", "i_pk", 1.0)):
                d = dict(tab.get(meas, []))
                resp[(tag, out)] = np.array([d.get(s, np.nan) * scale
                                             for s in range(1, n + 1)])
        return Q, T, M, resp

    def actual(x):
        """Run the real bench at x=(q_nC, tau, tm_ns) for both di/dt; return
        (Qrr_hi, T_hi, Qrr_lo, T_lo, I_RM_lo) with None for failed meas."""
        out = {}
        for tag, didt in (("hi", didt_hi), ("lo", didt_lo)):
            text, _ = run_ltspice(make_netlist(i_f, v_r, didt, [x[0]], [x[1]],
                                               dtmax, f"eval_{tag}",
                                               tm_list=[x[2]]), workdir)
            qv = _meas_value(text, "q_rr")
            tv = _meas_value(text, key_t)
            iv = _meas_value(text, "i_pk")
            out[tag] = (qv, tv * 1e9 if tv is not None else None, iv)
        return out

    def cost_from(vals):
        c = 0.0
        w = dict(zip(("Qh", "Th", "Ql", "Tl"), weights))
        for tag, tq, tt, kq, kt in (("hi", qrr_hi_nC, ttr_hi_ns, "Qh", "Th"),
                                    ("lo", qrr_lo_nC, ttr_lo_ns, "Ql", "Tl")):
            qv, tv, iv = vals[tag]
            c += w[kq] * (100.0 if qv is None
                          else ((qv - tq) / (0.1 * tq)) ** 2)
            if tv is None:  # t_rr target never reached (I_RM below target):
                c += w[kt] * (25.0 + ((max(0.0, 0.25 * i_f + 0.1 - (-iv))) / 1.0) ** 2)
            else:
                c += w[kt] * ((tv - tt) / (0.1 * tt)) ** 2
        return c

    # trust-region iteration: local 2x2x2 grid (exact trilinear surrogate,
    # small residual) around the current best; solve inside the box; accept
    # only if the REAL cost improves. Start at the two-knob identity point
    # (t_minor = q_0 reproduces the old model).
    x = np.array([35.0, 5.0, 35.0])
    rad = np.array([0.15, 0.40, 0.50])   # fractional radii (multiplicative)
    best_cost = cost_from(actual(x))
    print(f"start x={x} actual cost {best_cost:.3f}")
    for it in range(4):
        qg = (x[0] * (1 - rad[0]), x[0] * (1 + rad[0]))
        taug = (max(0.3, x[1] * (1 - rad[1])), x[1] * (1 + rad[1]))
        tmg = (max(3.0, x[2] * (1 - rad[2])), x[2] * (1 + rad[2]))
        Q, T, M, resp = grid(qg, taug, tmg)
        cs = {}
        for k, v in resp.items():
            nvalid = int((~np.isnan(v)).sum())
            if nvalid >= 6:
                cs[k] = trilinear_fit(Q, T, M, v)
        box = [(qg[0], qg[1]), (taug[0], taug[1]), (tmg[0], tmg[1])]

        def scost(xx):
            c = 0.0
            w = dict(zip(("Qh", "Th", "Ql", "Tl"), weights))
            for (tag, out), tq, kq in ((("hi", "Qrr"), qrr_hi_nC, "Qh"),
                                       (("lo", "Qrr"), qrr_lo_nC, "Ql")):
                if (tag, out) in cs:
                    c += w[kq] * ((trilinear_eval(cs[(tag, out)][0], *xx) - tq)
                                  / (0.1 * tq)) ** 2
            for (tag, out), tt, kt in ((("hi", "Ttr"), ttr_hi_ns, "Th"),
                                       (("lo", "Ttr"), ttr_lo_ns, "Tl")):
                if (tag, out) in cs:
                    c += w[kt] * ((trilinear_eval(cs[(tag, out)][0], *xx) - tt)
                                  / (0.1 * tt)) ** 2
            return c

        corners = ((0, 0, 0), (1, 1, 1), (0, 1, 0), (1, 0, 1))
        starts = [tuple(x)] + [tuple(b[c] for b, c in zip(box, cn))
                               for cn in corners]
        best = None
        for s in starts:
            r = minimize(scost, s, bounds=box, method="L-BFGS-B")
            if best is None or r.fun < best.fun:
                best = r
        xnew = np.clip(best.x, [b[0] for b in box], [b[1] for b in box])
        c_new = cost_from(actual(xnew))
        print(f"iter {it}: box q[{qg[0]:.1f},{qg[1]:.1f}] "
              f"tau[{taug[0]:.2f},{taug[1]:.2f}] tm[{tmg[0]:.1f},{tmg[1]:.1f}] "
              f"-> x={np.round(xnew, 3)} actual cost {c_new:.3f} "
              f"(best {best_cost:.3f})")
        if c_new < best_cost:
            x, best_cost = xnew, c_new
            rad *= 0.7
        else:
            rad *= 0.5
    q_nC, tau, tm_ns = x

    # 3) verification: waveforms + measurements at both di/dt for the 3-knob
    #    optimum and the two-knob references (A: LS, C: shape-strict)
    lines = [
        "== q_0/tau/t_minor extraction, both di/dt points fitted ==",
        f"targets: {didt_hi} A/us -> Qrr={qrr_hi_nC} nC, t_{tdef}={ttr_hi_ns} ns; "
        f"{didt_lo} A/us -> Qrr={qrr_lo_nC} nC, t_{tdef}={ttr_lo_ns} ns "
        f"(I_F={i_f} A, V_R={v_r} V), dTmax={dtmax*1e9:.1f} ns",
        f"3-knob fit: q_0={q_nC:.2f} n, tau={tau:.3f}, t_minor={tm_ns:.2f} ns "
        f"(cost {best_cost:.3f})",
    ]
    panels = {"hi": [], "lo": []}
    refs = [("2-knob C (34.8n/5)", 34.79, 5.0, None),
            ("3-knob fit", q_nC, tau, tm_ns)]
    for tag, didt, ttr_ns in (("hi", didt_hi, ttr_hi_ns),
                              ("lo", didt_lo, ttr_lo_ns)):
        for label, q, t, m in refs:
            tm_list = [m] if m is not None else None
            text, wave = run_ltspice(make_netlist(i_f, v_r, didt, [q], [t],
                                                  dtmax, label, tm_list=tm_list),
                                     workdir, want_waveform=True)
            Qrr_v = _meas_value(text, "q_rr")
            Ttr_v = _meas_value(text, key_t)
            ipk = _meas_value(text, "i_pk")
            Ttr_v = Ttr_v * 1e9 if Ttr_v is not None else None
            tgt_q = qrr_hi_nC if tag == "hi" else qrr_lo_nC
            tgt_t = ttr_hi_ns if tag == "hi" else ttr_lo_ns
            qs = f"{Qrr_v:.1f} nC ({(Qrr_v-tgt_q)/tgt_q*100:+.1f}%)" if Qrr_v else "FAIL"
            ts = (f"{Ttr_v:.2f} ns ({(Ttr_v-tgt_t)/tgt_t*100:+.1f}%)"
                  if Ttr_v else "FAIL (I_RM below target)")
            lines.append(f"  {didt} A/us {label}: q_0={q:.2f} n tau={t:.3f} "
                         f"t_minor={m if m else '-'} -> Qrr {qs}, "
                         f"t_{tdef} {ts}, I_RM={-ipk:.2f} A")
            panels[tag].append((label, q, t, m, wave))

    report = "\n".join(lines)
    print(report)
    os.makedirs(workdir, exist_ok=True)
    with open(os.path.join(workdir, "extract3k_report.txt"), "w") as f:
        f.write(report + f"\n\nworkdir: {workdir}\n")

    if plot_path:
        plot_waveforms_didt(panels, ttr_hi_ns, ttr_lo_ns, didt_hi, didt_lo,
                            plot_path)
        print(f"waveform plot: {plot_path}")
    return (q_nC, tau, tm_ns), report


def plot_waveforms_didt(panels, ttr_hi_ns, ttr_lo_ns, didt_hi, didt_lo,
                        plot_path):
    """Two stacked panels (high- and low-di/dt), each overlaying the
    recovery-current waveforms of the compared variants."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 1, figsize=(9, 8))
    for ax, tag, ttr_ns, didt in ((axes[0], "hi", ttr_hi_ns, didt_hi),
                                  (axes[1], "lo", ttr_lo_ns, didt_lo)):
        xhi = max(4.0 * ttr_ns, 80.0)
        for label, q, t, m, w in panels[tag]:
            tt = w["time"] * 1e9
            i = w["i(r4)"]
            idx = np.flatnonzero((i[:-1] > 0) & (i[1:] <= 0))
            if len(idx) == 0:
                raise RuntimeError(f"no recovery zero crossing for {label}")
            t0 = tt[idx[0]]
            par = f"q_0={q:.1f} n, τ={t:g}" + (f", t_minor={m:g} n" if m else "")
            ax.plot(tt - t0, i, label=f"{label} ({par})")
        ax.axhline(0, color="k", lw=0.5)
        ax.set_xlim(-15, xhi)
        ax.set_ylim(-14, 14)
        ax.set_title(f"{didt} A/µs (t_{ttr_ns:g} ns target)", fontsize=10)
        ax.set_ylabel("I(R4) [A]")
        ax.set_xlabel("t − t0 [ns]")
        ax.legend(loc="lower left", fontsize=8)
        ax.grid(alpha=0.3)
    fig.suptitle("Body-diode reverse recovery — two-knob vs three-knob "
                 "(t_minor) fits, both datasheet di/dt points")
    fig.tight_layout()
    fig.savefig(plot_path, dpi=150)
    plt.close(fig)
    return plot_path


if __name__ == "__main__":
    root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    out = os.path.join(root, "out")
    # ISC040N10NM8: fit q_0/tau/t_minor to BOTH datasheet recovery points
    extract3k(i_f=12.5, v_r=50.0,
              didt_hi=1000, qrr_hi_nC=247.0, ttr_hi_ns=25.0,
              didt_lo=100, qrr_lo_nC=36.0, ttr_lo_ns=29.0, tdef="25%",
              workdir=os.path.join(out, "recov_extract3k_run"),
              plot_path=os.path.join(out, "RECOV_3k_waveforms.png"))
