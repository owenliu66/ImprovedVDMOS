"""Datasheet-driven extraction for the DiodeRecov body-diode recovery model.

Model (the "transit-time pole" form, session 25b — see
out/RECOV_POLE_FIT.md):

    C1 q_j0 0 {t_minor}      storage node; Roff=1 leak -> storage time
                             constant = t_minor, V(q_j0)_ss = I_F volts,
                             stored charge = I_F * t_minor
    B1 0 q_j0 I=-I(R2)       injection mirror (branch current)
    D2 0 q_j0 ClampD         clamps V(q_j0) >= 0  (charge exhaustion)
    E1 + R3 + C3             ONE POLE: V(q_j) = filtered V(q_j0),
                             time constant t_pole = R3*C3 (R3 = 1 ohm)
                             -> finite carrier transit/sweep-out memory
    B2 Ki Ai I=max(0, V(q_j)*(V(Ki)-V(Ai))/tau)
                             conduction drive from the FILTERED node;
                             tau is the snap-threshold parameter

Three knobs, all free: t_minor (carrier lifetime), t_pole (transit memory),
tau (snap threshold). q_0 is GONE -- it was a pure gauge (observables
depended only on t_minor and on V_R/(q_0*tau)).

Measurement doctrine (see out/RECOV_INTRINSIC_FIT.md): the datasheet Qrr
does NOT include the C_oss excursion, so fit the INTRINSIC branch current
I(x1:R2), which carries no C2 displacement current. Terminal-path q_rr/t_rr
are reported for reference only.

Protocol
  1. One LTspice -b process per di/dt point: a single `.step param k` index
     sweep with mixed-radix decomposition into (tau, t_minor, t_pole).
     Deterministic index->parameter mapping, no `.step` nesting ambiguity.
  2. Log-log multilinear surrogates for (Qrr, t_rr) at each di/dt point.
  3. Solve the inverse problem on the surrogate (multi-start, instant).
  4. Polish on real runs: Levenberg-Marquardt/TRF in log space with
     diff_step = 0.03 (see "Finite-difference step" below) and a generous
     step budget; the tau direction is a flat valley coupled to both other
     knobs, so the budget matters.
  5. Verify at every datasheet di/dt point; report residuals, the
     fix-and-refit identifiability penalties, and the terminal-path numbers.
  6. Plot intrinsic vs terminal recovery current at each point.

Finite-difference step (measured, session 25c): diff_step = 1e-4 makes the
fit WORSE (tau stuck at 0.119, RMS 5.57 %) because the cost differences
across a 1e-4 relative parameter step (~1e-6) fall below the simulator's own
measurement noise (~1e-4 relative on q_rri/t_rri), so the Jacobian is
noise-dominated. diff_step = 0.03 (3 % in log space) puts the cost
difference ~1e-2, far above the noise, and reaches the valley bottom.
Tighter ftol/xtol/gtol do not help for the same reason.

dTmax is NOT used (measured, session 25c): with the pole in the circuit the
adaptive controller resolves the snap on its own. dTmax = none / 2.5n / 1n /
0.5n / 0.2n agree to <= 0.01 % on q_rri and t_rri, and dropping dTmax is
also faster (4088 vs 8368 iterations at 100 A/us). The session-24 "dTmax
mandatory" finding belonged to the pole-free model, where the snap was a
hard discontinuity.

Total cost: two batched grid runs (~30 s each) + ~40 polish runs + 2
verification runs -> a few minutes.
"""
from __future__ import annotations

import glob
import json
import math
import os
import re
import subprocess
from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import least_squares

# ---------------------------------------------------------------------------
# device specification
# ---------------------------------------------------------------------------

@dataclass
class RecovSpec:
    """Datasheet recovery inputs for one part + the grid to search."""
    part: str
    i_f: float = 12.5                    # forward current [A]
    v_r: float = 50.0                    # reverse voltage [V]
    # (di/dt [A/us], Qrr [nC], t_rr [ns]) -- t_rr = 10 % of I_RM window
    points: tuple = ((1000.0, 247.0, 25.0), (100.0, 36.0, 29.0))
    # subckt DC branch (already fitted by the diode protocol)
    is_diode: float = 1e-6
    rs_diode: float = 10e-3
    n_diode: float = 1.4
    # forward-recovery branch (not identifiable from recovery data; fixed)
    k_chg: float = 5.0
    r_max: float = 0.0
    # search grid, log-spaced: (n, lo, hi) per knob
    grid_tau: tuple = (5, 1e-5, 1e-2)
    grid_tminor: tuple = (5, 12e-9, 40e-9)
    grid_tpole: tuple = (5, 0.4e-9, 2.0e-9)
    # parameter bounds for the polish fit
    bounds: dict = field(default_factory=lambda: {
        "tau": (1e-6, 5.0), "t_minor": (5e-9, 60e-9), "t_pole": (0.2e-9, 6e-9)})
    cds_expr_file: str = "ISC040N10NM8_cds_expr.txt"


ISC040 = RecovSpec(part="ISC040N10NM8")

# ---------------------------------------------------------------------------
# bench generation
# ---------------------------------------------------------------------------

NETLIST_TMPL = """* extract_recov bench -- {label}
* antiparallel DiodeRecov pair + stiff V_R: the DC operating point IS the
* pre-event state, so no .ic seeding is needed (no DC current in R3 =>
* V(q_j) = V(q_j0)).
.param i_test={i_f} didt={didt} tr=i_test*2/didt*1e-6 tf=tr
.param k_chg={k_chg} r_max={r_max}
V1 N001 0 {v_r}
I1 0 SW PULSE({{-i_test}} {{i_test*9}} {{tr}} {{tr*5}} {{tf}} 10m 20m)
X1 0 SW DiodeRecov params: Is={is_d} Rs={rs_d} N={n_d} tau={{tau_p}} t_minor={{tm_p}} \
k_chg={{k_chg}} r_max={{r_max}}
X2 SW N001 DiodeRecov params: Is={is_d} Rs={rs_d} N={n_d} tau={{tau_p}} t_minor={{tm_p}} \
k_chg={{k_chg}} r_max={{r_max}}

.subckt DiodeRecov A K
B2 Ki Ai I=max(0,V(q_j)*(V(Ki)-V(Ai))/tau)
C1 q_j0 0 {{tm_p}} Rser=0
B1 0 q_j0 I=-I(R2)
D2 0 q_j0 ClampD
R1 N001 Ki {{Rs}}
R2 Ai A 1u
B3 N001 K V=if(r_max>0,-I(R2)*(r_max/(1+k_chg*V(q_j))+1u),0)
C2 Ki A {cds}
D1 Ai Ki MyDiode
E1 N002 0 q_j0 0 1
R3 q_j N002 1
C3 q_j 0 {{tp_p}}
.model ClampD D(Ron=1u,Roff=1,Vfwd=0)
.model MyDiode D(Is=Is,Rs=0,N=N)
.ends DiodeRecov

.model D D
.lib {diolib}
.tran 0 {{tr*5+50n}} 0 convreport
.options reltol=1e-4 method=gear solver=alt
.model IdealSW SW(Ron=1m,Roff=10Meg,vt=5,vh=-2)
.model IdealD D(Ron=1m,Roff=1Meg)
.meas didt_max MAX d(Ix(x1:a))/d(time)*1e-9
.meas i_f_min MIN Ix(x1:a)
.meas i_fi_min MIN -I(x1:R2)
.meas t_rr TRIG Ix(x1:a)=0 TD=0 FALL=1 TARG Ix(x1:a)={{i_f_min*0.1}} TD=0 RISE=1
.meas t_rri TRIG Ix(x1:a)=0 TD=0 FALL=1 TARG -I(x1:R2)={{i_fi_min*0.1}} TD=0 RISE=1
.meas t_rri25 TRIG Ix(x1:a)=0 TD=0 FALL=1 TARG -I(x1:R2)={{i_fi_min*0.25}} \
TD=0 RISE=1
.meas q_rr INTEG -Ix(x1:a)*1e9 TRIG Ix(x1:a)=0 TD=0 FALL=1 TARG Ix(x1:a)={{i_f_min*0.1}} TD=0 RISE=1
.meas q_rri INTEG I(x1:R2)*1e9 TRIG Ix(x1:a)=0 TD=0 FALL=1 TARG \
-I(x1:R2)={{i_fi_min*0.1}} TD=0 RISE=1
{steps}.end
"""

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA = os.path.join(ROOT, "data")
OUT = os.path.join(ROOT, "out")


def find_ltspice() -> str:
    cands = []
    for base in (r"C:\Program Files\ADI\LTspice", r"C:\Program Files (x86)\ADI\LTspice",
                 os.path.expandvars(r"%LOCALAPPDATA%\LTspice")):
        cands += glob.glob(os.path.join(base, "LTspice.exe"))
    if not cands:
        raise RuntimeError("LTspice.exe not found -- set the path in find_ltspice()")
    return cands[0]


LTSPICE = find_ltspice()
DIOLIB = os.path.expandvars(r"%LOCALAPPDATA%\LTspice\lib\cmp\standard.dio")


def _cds_clause(fname: str) -> str:
    """Measured C_ds charge-master ramp-sum (from the caps record)."""
    try:
        with open(os.path.join(DATA, fname)) as f:
            return "Q=" + f.read().strip().split("Q=", 1)[1]
    except (OSError, IndexError):
        return "500p"          # placeholder; recovery numbers move < 1 %


def _spice_num(v: float) -> str:
    """SPICE-safe scientific literal (no suffixes inside {} expressions)."""
    return f"{v:.10e}"


def _grid_vals(grid: tuple) -> np.ndarray:
    """(n, lo, hi) -> n log-spaced values."""
    n, lo, hi = grid
    return np.exp(np.linspace(math.log(lo), math.log(hi), n))


def index_step_block(spec: RecovSpec) -> str:
    """Single `.step param k` index sweep with mixed-radix decomposition
    into (tau, t_minor, t_pole). Verified in outputs/vdmos_recov/decomp.net:
    global step params DO drive values inside the subckt, including the
    C3 q_j 0 {tp_p} value. GOTCHA: every referenced param must be defined
    or LTspice dies at parse with 'No such parameter defined.'

    Index mapping (k is 0-based; log tables restart at every step):
      nA = k % na            (tau, fastest)
      nB = (k // na) % nb    (t_minor)
      nC = k // (na*nb)      (t_pole, slowest)
    """
    tau_g = _grid_vals(spec.grid_tau)
    tm_g = _grid_vals(spec.grid_tminor)
    tp_g = _grid_vals(spec.grid_tpole)
    na, nb, nc = len(tau_g), len(tm_g), len(tp_g)

    def dln(g):
        return math.log(g[-1] / g[0]) / (len(g) - 1)

    lines = [
        ".step param k list " + " ".join(str(i) for i in range(na * nb * nc)),
        f".param nA=mod(k,{na})",
        f".param nB=mod(int(k/{na}),{nb})",
        f".param nC=int(k/{na * nb})",
        f".param tau_p=exp({_spice_num(math.log(tau_g[0]))}+nA*{_spice_num(dln(tau_g))})",
        f".param tm_p=exp({_spice_num(math.log(tm_g[0]))}+nB*{_spice_num(dln(tm_g))})",
        f".param tp_p=exp({_spice_num(math.log(tp_g[0]))}+nC*{_spice_num(dln(tp_g))})",
    ]
    return "\n".join(lines) + "\n"


def grid_points(spec: RecovSpec):
    """(tau, t_minor, t_pole) per LTspice step, in step order (step 1 = k 0)."""
    tau_g = _grid_vals(spec.grid_tau)
    tm_g = _grid_vals(spec.grid_tminor)
    tp_g = _grid_vals(spec.grid_tpole)
    na, nb, nc = len(tau_g), len(tm_g), len(tp_g)
    pts = []
    for k in range(na * nb * nc):
        pts.append((tau_g[k % na], tm_g[(k // na) % nb], tp_g[k // (na * nb)]))
    return pts


def make_netlist(spec: RecovSpec, didt: float, *, tau: float | None = None,
                 t_minor: float | None = None, t_pole: float | None = None,
                 label: str = "bench", steps: str | None = None) -> str:
    """One bench netlist at one di/dt point.

    steps=None: single run, knob values passed literally via
    .param tau_p/tm_p/tp_p (tau/t_minor/t_pole required).
    steps=<text>: the step/param block is appended verbatim (the single
    index sweep from index_step_block); knob literals are ignored.
    """
    if steps is None:
        if tau is None or t_minor is None or t_pole is None:
            raise TypeError("use make_netlist(spec, didt, tau=..., "
                            "t_minor=..., t_pole=...)")
        steps = (f".param tau_p={_spice_num(tau)}\n"
                 f".param tm_p={_spice_num(t_minor)}\n"
                 f".param tp_p={_spice_num(t_pole)}\n")
    return NETLIST_TMPL.format(
        label=label, i_f=spec.i_f, didt=didt, v_r=spec.v_r,
        k_chg=spec.k_chg, r_max=spec.r_max,
        is_d=_spice_num(spec.is_diode), rs_d=_spice_num(spec.rs_diode),
        n_d=spec.n_diode, cds=_cds_clause(spec.cds_expr_file),
        diolib=DIOLIB, steps=steps)


# ---------------------------------------------------------------------------
# LTspice run + log parsing
# ---------------------------------------------------------------------------

def run_ltspice(netlist: str, workdir: str, want_waveform=False):
    """Run LTspice -b on a scratch netlist; parse the log (and the ASCII raw
    when want_waveform), then DELETE the netlist and all run artifacts."""
    os.makedirs(workdir, exist_ok=True)
    path = os.path.join(workdir, "bench.net")
    with open(path, "w") as f:
        f.write(netlist)
    cmd = [LTSPICE, "-b"] + (["-ascii"] if want_waveform else []) + [path]
    subprocess.run(cmd, cwd=workdir, timeout=600, check=False)
    log = path[:-4] + ".log"
    with open(log, errors="ignore") as f:
        text = f.read()
    if "Total elapsed time" not in text:
        raise RuntimeError("LTspice run failed:\n" + text[-2000:])
    wave = parse_ascii_raw(path[:-4] + ".raw",
                           ["time", "Ix(x1:a)", "I(x1:R2)"]) if want_waveform else None
    for leftover in glob.glob(path[:-4] + ".*"):
        os.remove(leftover)
    return text, wave


def parse_ascii_raw(raw_path, traces):
    """Parse an LTspice ASCII .raw (from `ltspice -b -ascii`). Returns dict
    name -> np.array, keys lowercased as given; name matching is
    case-insensitive (raw writes I(X1:R2)-style names)."""
    with open(raw_path, errors="ignore") as f:
        lines = f.read().splitlines()
    nv = int(next(ln for ln in lines
                  if ln.startswith("No. Variables")).split(":")[1])
    vi = next(i for i, ln in enumerate(lines) if ln.startswith("Variables:"))
    names = [lines[i].split()[1].strip('"') for i in range(vi + 1, vi + 1 + nv)]
    lower = {n.lower(): n for n in names}
    vals_i = next(i for i, ln in enumerate(lines) if ln.strip() == "Values:")
    # LTspice ASCII raw: one point per group -- index line, then one tab-
    # prefixed value per line. Tokenize the whole block and walk it.
    tokens = [t for ln in lines[vals_i + 1:] for t in ln.split()]
    rows = []
    i = 0
    while i + nv < len(tokens) + 1:
        if not re.fullmatch(r"\d+", tokens[i]):
            break
        rows.append([float(x) for x in tokens[i + 1:i + 1 + nv]])
        i += 1 + nv
    arr = np.array(rows)
    # value column k corresponds to variable k (the point index column is dropped)
    want = {}
    for t in traces:
        key = lower.get(t.lower())
        if key is not None:
            want[t.lower()] = arr[:, names.index(key)]
    return want


def parse_stepped(text):
    """Return {meas_name: [(step, value), ...]} from stepped .measure tables.
    Chained {i_fi_min} targets DO resolve in stepped runs (session-25c
    correction: the session-25a 'dependent block dropped' claim was a
    bookkeeping error -- the netlist never contained the line)."""
    out = {}
    pat = r"Measurement: (\w+)\r?\n[^\r\n]*\r?\n((?:\s+\d+\s+[^\r\n]+\r?\n)+)"
    for m in re.finditer(pat, text):
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


def _meas_value(text, name):
    """Single-run .meas value. The value token is the one immediately before
    FROM, then split('=')[-1] (the 'INTEG(I(x1:R2)*1e9)=value' form glues
    them; the ':' inside I(x1:R2) breaks naive 'name[:=]value' regexes)."""
    m = re.search(rf"^{name}\b[^\n]*", text, re.M)
    if not m:
        return None
    toks = m.group(0).split()
    v = toks[toks.index("FROM") - 1] if "FROM" in toks else toks[-1]
    try:
        return float(v.split("=")[-1])
    except ValueError:      # 'failed'
        return None


# ---------------------------------------------------------------------------
# surrogates (log-log multilinear, full 3-way interactions)
# ---------------------------------------------------------------------------

def _terms(L):
    """Design matrix for multilinear-in-log(params) interpolation."""
    a, b, c = L[:, 0], L[:, 1], L[:, 2]
    return np.column_stack([np.ones_like(a), a, b, c, a * b, a * c, b * c,
                            a * b * c])


def log_surrogate(X, z):
    """Fit log(z) multilinearly in log(X); X is (n, 3) params, z responses.
    Returns (coeffs, max |relative| residual). NaN z rows are dropped."""
    keep = np.isfinite(z) & np.all(np.isfinite(X), axis=1) & (z > 0)
    L = np.log(X[keep])
    T = _terms(L)
    c, *_ = np.linalg.lstsq(T, np.log(z[keep]), rcond=None)
    rel = np.abs(np.exp(T @ c - np.log(z[keep])) - 1.0).max()
    return c, float(rel)


def surrogate_eval(c, x_log):
    return float(np.exp(_terms(np.asarray(x_log, dtype=float).reshape(1, 3)) @ c)[0])


# ---------------------------------------------------------------------------
# real-run evaluation (cached)
# ---------------------------------------------------------------------------

class Bench:
    """One real bench evaluation at one di/dt point; results cached."""

    def __init__(self, spec: RecovSpec, didt: float, workdir: str):
        self.spec, self.didt, self.workdir = spec, didt, workdir
        self.cache = {}

    def measure(self, x_log):
        key = tuple(round(v, 10) for v in x_log)
        if key in self.cache:
            return self.cache[key]
        tau, tm, tp = (math.exp(v) for v in x_log)
        text, _ = run_ltspice(make_netlist(self.spec, self.didt, tau=tau,
                                           t_minor=tm, t_pole=tp,
                                           label=f"b{self.didt:g}"),
                              self.workdir)
        out = {
            "q_rri": _meas_value(text, "q_rri"),
            "t_rri": _meas_value(text, "t_rri"),
            "t_rri25": _meas_value(text, "t_rri25"),
            "q_rr": _meas_value(text, "q_rr"),
            "t_rr": _meas_value(text, "t_rr"),
            "irm": _meas_value(text, "i_fi_min"),
            "didt_max": _meas_value(text, "didt_max"),
        }
        # a recovery event that never reaches the 10 % target fails t_rri;
        # a degenerate I_RM (< 0.2 A) is not a usable snap either
        if out["q_rri"] is None or out["t_rri"] is None or \
           out["irm"] is None or out["irm"] > -0.2:
            out = None
        self.cache[key] = out
        return out


# ---------------------------------------------------------------------------
# inverse solve on the surrogate, then polish on real runs
# ---------------------------------------------------------------------------

def _targets(spec: RecovSpec):
    """[(bench index, didt, q_nC, t_ns), ...] in residual order."""
    return [(i, p[0], p[1], p[2]) for i, p in enumerate(spec.points)]


def surrogate_solve(spec, surrogates, x_lo, x_hi, n_starts=12):
    """Multi-start least squares on the surrogate surface (instant).
    surrogates: {(bench, 'q'): c, (bench, 't'): c} from log_surrogate."""
    tg = _targets(spec)

    def res(x):
        r = []
        for i, _didt, q, t in tg:
            r.append(surrogate_eval(surrogates[(i, "q")], x) / q - 1.0)
            r.append(surrogate_eval(surrogates[(i, "t")], x) / t - 1.0)
        return r

    rng = np.random.default_rng(0)
    starts = []
    for _j in range(n_starts):
        starts.append(np.array([x_lo[k] + (x_hi[k] - x_lo[k]) * rng.random()
                                for k in range(3)]))
    best = None
    for x0 in starts:
        r = least_squares(res, x0, bounds=(x_lo, x_hi), method="trf")
        if best is None or r.cost < best.cost:
            best = r
    return best.x, best.cost


def polish(spec, benches, x0, weights=(1.0, 1.0, 1.0, 1.0), max_nfev=200):
    """TRF in log space on REAL runs. diff_step=0.03: measured optimum
    (module docstring); generous max_nfev because tau is a flat valley
    coupled to both other knobs. status==3 with tau far from the surrogate
    prediction means 'stopped on tolerance', not 'converged' -- the caller
    reports it, the multi-start + surrogate seed makes it unlikely."""
    tg = _targets(spec)

    def res(x):
        r = []
        for (i, _didt, q, t), w in zip(tg, weights, strict=True):
            m = benches[i].measure(x)
            if m is None:
                r.extend([5.0 * w, 5.0 * w])
                continue
            r.append((m["q_rri"] / q - 1.0) * w)
            r.append((m["t_rri"] * 1e9 / t - 1.0) * w)
        return r

    x_lo = [math.log(spec.bounds["tau"][0]), math.log(spec.bounds["t_minor"][0]),
            math.log(spec.bounds["t_pole"][0])]
    x_hi = [math.log(spec.bounds["tau"][1]), math.log(spec.bounds["t_minor"][1]),
            math.log(spec.bounds["t_pole"][1])]
    return least_squares(res, np.clip(x0, x_lo, x_hi), bounds=(x_lo, x_hi),
                         method="trf", x_scale="jac", diff_step=0.03,
                         max_nfev=max_nfev)


def polish_fix_one(spec, benches, x_opt, fixed_knob, fixed_val,
                   weights=(1.0, 1.0, 1.0, 1.0)):
    """Fix-and-refit identifiability (session-25b doctrine): hold knob
    `fixed_knob` at `fixed_val` (log-space; pass optimum x perturbation --
    fixing AT the optimum is tautological) and refit the two free knobs.
    Returns (rms_rel_pct, x)."""
    tg = _targets(spec)
    fixed = np.array(x_opt, dtype=float)
    fixed[fixed_knob] = fixed_val
    free_idx = [j for j in range(3) if j != fixed_knob]

    def res(u):
        x = fixed.copy()
        for j, k in enumerate(free_idx):
            x[k] = u[j]
        r = []
        for (i, _didt, q, t), w in zip(tg, weights, strict=True):
            m = benches[i].measure(x)
            if m is None:
                r.extend([5.0 * w, 5.0 * w])
                continue
            r.append((m["q_rri"] / q - 1.0) * w)
            r.append((m["t_rri"] * 1e9 / t - 1.0) * w)
        return r

    u0 = [x_opt[k] for k in free_idx]
    knob = tuple(("tau", "t_minor", "t_pole")[i] for i in free_idx)
    lo = [math.log(spec.bounds[k][0]) for k in knob]
    hi = [math.log(spec.bounds[k][1]) for k in knob]
    r = least_squares(res, u0, bounds=(lo, hi), method="trf", x_scale="jac",
                      diff_step=0.03, max_nfev=120)
    x = fixed.copy()
    for j, k in enumerate(free_idx):
        x[k] = r.x[j]
    rms = 100.0 * math.sqrt(2.0 * r.cost / len(tg) / 2)
    return rms, x


def rms_rel_pct(res_vec):
    return 100.0 * math.sqrt(np.mean(np.square(res_vec)))


# ---------------------------------------------------------------------------
# plotting
# ---------------------------------------------------------------------------

def plot_waveforms(waves, plot_path):
    """waves: [(didt, {'time':..,'ix(x1:a)':..,'i(x1:r2)':..}), ...].
    Two panels (one per di/dt): intrinsic I(x1:R2) vs terminal Ix(x1:a),
    time relative to the first falling zero crossing of the terminal current."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(len(waves), 1, figsize=(9, 4 * len(waves)),
                             squeeze=False)
    for ax, (didt, w) in zip(axes[:, 0], waves, strict=True):
        t = w["time"] * 1e9
        i_term = w["ix(x1:a)"]
        i_int = w["i(x1:r2)"]
        idx = np.flatnonzero((i_term[:-1] > 0) & (i_term[1:] <= 0))
        if len(idx) == 0:
            raise RuntimeError(f"no recovery zero crossing at {didt} A/us")
        t0 = t[idx[0]]
        m = (t >= t0 - 15) & (t <= t0 + max(80.0, 3.0 * 25.0))
        ax.plot(t[m] - t0, i_int[m], label="intrinsic  I(x1:R2)")
        ax.plot(t[m] - t0, i_term[m], "--", label="terminal  Ix(x1:a)")
        ax.axhline(0, color="k", lw=0.5)
        ax.set_title(f"{didt:g} A/\u00b5s")
        ax.set_ylabel("recovery current [A]")
        ax.set_xlabel("t - t0 [ns]")
        ax.legend(loc="lower left", fontsize=8)
        ax.grid(alpha=0.3)
    fig.suptitle("Body-diode reverse recovery -- intrinsic vs terminal path "
                 "(pole-gauge extraction)")
    fig.tight_layout()
    fig.savefig(plot_path, dpi=150)
    plt.close(fig)
    return plot_path


# ---------------------------------------------------------------------------
# top-level extraction
# ---------------------------------------------------------------------------

def extract(spec: RecovSpec = ISC040, workdir: str | None = None,
            plot_path: str | None = None, identifiability: bool = True,
            weights=(1.0, 1.0, 1.0, 1.0)):
    """Full datasheet-driven extraction. Returns (params dict, report text).

    params: {"tau":.., "t_minor":.., "t_pole":.., "rms_rel":..,
             "per_point": [{didt, q_rri, t_rri, t_rri25, q_rr, t_rr, irm}, ...]}
    """
    workdir = workdir or os.path.join(OUT, f"recov_extract_{spec.part.lower()}_run")
    os.makedirs(workdir, exist_ok=True)
    lines = [f"== DiodeRecov pole-gauge extraction: {spec.part} ==",
             f"I_F={spec.i_f} A, V_R={spec.v_r} V; targets "
             + ", ".join(f"{d:g} A/us -> {q} nC / {t} ns"
                         for d, q, t in spec.points)]
    names = ("tau", "t_minor", "t_pole")

    # 1) batched grid runs (one LTspice process per di/dt point)
    tg = _targets(spec)
    pts = grid_points(spec)
    X = np.array(pts)
    Y = {}
    for i, didt, _, _ in tg:
        text, _ = run_ltspice(make_netlist(spec, didt,
                                           steps=index_step_block(spec),
                                           label=f"grid{didt:g}"), workdir)
        tab = parse_stepped(text)
        n = len(pts)
        for meas, key, scale in (("q_rri", "q", 1.0), ("t_rri", "t", 1e9)):
            d = dict(tab.get(meas, []))
            Y[(i, key)] = np.array([d.get(s, np.nan) * scale
                                    for s in range(1, n + 1)])
        miss = int(np.isnan(Y[(i, "q")]).sum() + np.isnan(Y[(i, "t")]).sum())
        if miss:
            lines.append(f"grid @ {didt:g} A/us: {miss} failed .meas rows")

    # 2) log-log surrogates
    surrogates = {}
    for (i, key), z in Y.items():
        c, rel = log_surrogate(X, z)
        surrogates[(i, key)] = c
        lines.append(f"surrogate {('q','t')[key == 't']} @ "
                     f"{spec.points[i][0]:g} A/us: max rel resid {rel*100:.2f} %")

    # 3) inverse solve on the surrogate
    x_lo = [math.log(spec.bounds[n][0]) for n in names]
    x_hi = [math.log(spec.bounds[n][1]) for n in names]
    x_surr, cost = surrogate_solve(spec, surrogates, x_lo, x_hi)

    # best grid point (real cost) as an extra seed
    def real_cost(x):
        r = 0.0
        for i, _didt, q, t in tg:
            m = benches[i].measure(x)
            if m is None:
                return 1e9
            r += (m["q_rri"] / q - 1) ** 2 + (m["t_rri"] * 1e9 / t - 1) ** 2
        return r

    benches = [Bench(spec, p[0], workdir) for p in spec.points]
    grid_costs = [real_cost(p) for p in pts]
    x_grid = np.log(pts[int(np.argmin(grid_costs))])
    x_seed = np.array([0.1, 17e-9, 2e-9])   # session-25c verified seed

    # 4) polish on real runs, multi-start
    best = None
    for tag, x0 in (("surr", x_surr), ("grid", x_grid), ("seed", x_seed)):
        r = polish(spec, benches, x0, weights=weights)
        rms = 100.0 * math.sqrt(2.0 * r.cost / (2 * len(tg)))
        lines.append(f"polish[{tag}]: tau={math.exp(r.x[0]):.4g} "
                     f"t_minor={math.exp(r.x[1])*1e9:.3f} n "
                     f"t_pole={math.exp(r.x[2])*1e9:.3f} n "
                     f"RMS {rms:.3f} % nfev={r.nfev} status={r.status}")
        if best is None or r.cost < best.cost:
            best = r
    x_opt = best.x
    rms_opt = 100.0 * math.sqrt(2.0 * best.cost / (2 * len(tg)))
    tau, tm, tp = (math.exp(v) for v in x_opt)
    lines.append(f"OPTIMUM: tau={tau:.5g} t_minor={tm*1e9:.3f} n "
                 f"t_pole={tp*1e9:.4f} n  RMS {rms_opt:.3f} % "
                 f"(nfev={best.nfev}, status={best.status})")

    # 5) verification at every datasheet point (intrinsic + terminal path)
    per_point = []
    for i, didt, q_t, t_t in tg:
        m = benches[i].measure(x_opt)
        if m is None:
            lines.append(f"VERIFY FAILED at {didt:g} A/us")
            continue
        per_point.append(dict(didt=didt, q_rri=m["q_rri"],
                              t_rri=m["t_rri"] * 1e9,
                              t_rri25=(m["t_rri25"] or float("nan")) * 1e9,
                              q_rr=m["q_rr"], t_rr=(m["t_rr"] or float("nan")) * 1e9,
                              irm=-m["irm"]))
        lines.append(
            f"@ {didt:g} A/us: q_rri {m['q_rri']:.1f} nC "
            f"({(m['q_rri']/q_t-1)*100:+.1f} %), "
            f"t_rri {m['t_rri']*1e9:.2f} ns ({(m['t_rri']*1e9/t_t-1)*100:+.1f} %), "
            f"t_rri25 {m['t_rri25']*1e9:.2f} ns, I_RM {-m['irm']:.2f} A | "
            f"terminal q_rr {m['q_rr']:.1f} nC, t_rr {m['t_rr']*1e9:.2f} ns")

    # 6) fix-and-refit identifiability at PERTURBED fixed values (quote
    #    penalties, not marginal sigma; fixing at the optimum is tautological)
    if identifiability:
        def disp(k):
            v = math.exp(x_opt[k])
            return f"{v:.4g}" if k == 0 else f"{v*1e9:.3f} n"
        for k, nm in enumerate(names):
            for f in (0.5, 2.0):
                rms_fix, _ = polish_fix_one(spec, benches, x_opt, k,
                                            x_opt[k] + math.log(f),
                                            weights=weights)
                lines.append(
                    f"identifiability: fix {nm} at {disp(k)} x{f:g}, refit "
                    f"other two -> RMS {rms_fix:.2f} % "
                    f"(penalty {rms_fix - rms_opt:+.2f} pp)")

    report = "\n".join(lines)
    print(report, flush=True)
    with open(os.path.join(workdir, "extract_report.txt"), "w") as f:
        f.write(report + f"\n\nworkdir: {workdir}\n")
    with open(os.path.join(workdir, "extract_results.json"), "w") as f:
        json.dump({"part": spec.part, "tau": tau, "t_minor": tm,
                   "t_pole": tp, "rms_rel_pct": rms_opt,
                   "per_point": per_point}, f, indent=1)

    # 7) waveform plot (intrinsic vs terminal at each di/dt)
    if plot_path:
        waves = []
        for _i, didt, _, _ in tg:
            _, w = run_ltspice(make_netlist(spec, didt, tau=tau, t_minor=tm,
                                            t_pole=tp, label=f"wave{didt:g}"),
                               workdir, want_waveform=True)
            waves.append((didt, w))
        plot_waveforms(waves, plot_path)
        print(f"waveform plot: {plot_path}", flush=True)

    return dict(tau=tau, t_minor=tm, t_pole=tp, rms_rel=rms_opt,
                per_point=per_point), report


if __name__ == "__main__":
    extract(ISC040,
            workdir=os.path.join(OUT, "recov_extract_isc040_run"),
            plot_path=os.path.join(OUT, "RECOV_pole_extract_waveforms.png"))
