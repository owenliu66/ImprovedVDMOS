# Curve Digitizer — Status

*Snapshot date: 2026-10-11; updated with two-point calibration and per-axis minor
ticks. Companion to [README.md](README.md); the project-level status lives in
[../../STATUS.md](../../STATUS.md).*

## What exists

```
src/curve_digitizer/        2,312 lines of Python, 13 modules
├── __init__.py       (13)  package exports, `python -m curve_digitizer`
├── __main__.py        (4)  entry point
├── app.py           (205)  main Tk window: PDF load, page nav, crop DPI, rubber-band region select
├── pdf_source.py     (61)  PyMuPDF page/region rendering (region re-rendered, not cropped)
├── transform.py     (182)  AxisTransform: pixel <-> data, linear/log, tick generation, from_points
├── grid.py          (117)  adjustable Cartesian grid overlay derived from the transform
├── snapping.py      (325)  dark-pixel mask -> skeleton -> sub-pixel centerline cloud -> KDTree
├── corners.py       (323)  90-degree tick/axis corner detection for calibration clicks (own mask)
├── calibration.py   (160)  click-two-points state machine: arming, crosshair artists, validation
├── markers.py       (235)  add/drag/select/delete, X/Y key axis-lock, marker artists, input gate
├── digitizer_window.py (475) pop-up: plot, axis panel, grid panel, snap controls, Export
├── table.py         (173)  Index|X|Y|Condition Treeview, in-place cell editing, 2-way sync
└── export.py         (39)  CSV via native Save As, default `Projects/ImprovedVDMOS/data`
tests/test_digitizer.py (350)  corner snapping, anchored transforms, per-axis minors, headless Tk
```

Runtime dependencies are declared in `pyproject.toml` as the `digitizer` extra
(`pip install -e .[digitizer]`: `pymupdf>=1.24`, `scikit-image>=0.22`), and the
console script `curve-digitizer` is registered alongside `vdmos-fit`.

## What it currently does

- **Region selection** — page navigation plus a rubber-band rectangle; each region
  opens its own pop-up, re-rendered from the PDF at crop DPI (default 300) so the
  trace stays sharp instead of being cropped out of the low-res page bitmap.
- **Coordinate model** — the image is always drawn in crop-pixel space
  (`imshow(extent=(0, W, H, 0))`); datasheet units exist only in tick labels, the
  table and the CSV. A log axis is mapped as pixels → `log10(data)`, so a printed
  log axis is never geometrically distorted.
- **Calibration** — pixel anchors of the two typed values are stored per axis
  (`x_anchors`/`y_anchors`); the default is the crop-edge pair, preserving the
  original edge-based workflow. **Calibrate X/Y** arms a two-point mode: click
  the two printed ticks (only the along-axis coordinate is used; click order is
  directional so reversed axes work; Esc or re-pressing the button cancels;
  clicks closer than 5 px are rejected). `AxisTransform.from_points` builds the
  map; for log axes the two points fix decade scale and offset. `_apply_axes`
  always rebuilds from anchors + typed values, so "Apply grid" cannot silently
  discard a calibration. Markers are pixel-anchored: recalibrating re-values
  table/CSV without moving markers.
- **Corner snapping for calibration clicks** — `corners.py` snaps each click to
  the vertex where a tick stroke crosses the axis line: longest-horizontal-run
  axis band, connected components touching it, covariance principal direction
  within 12° of perpendicular (rejects shallow curve crossings, arrowheads,
  text), grayscale threshold-crossing interpolation for sub-pixel centres in the
  `imshow` display convention (pixel `j` spans `[j, j+1]`). Fallbacks: lone
  perpendicular stroke → snap the along-axis coordinate only; no ink → raw
  click. Toggleable via *snap clicks to tick corners*. Its mask is deliberately
  separate from the snapper's (frames/ticks must survive) but shares the
  dark-threshold slider.
- **Axes** — values at the two anchor points (or crop edges), tick interval
  (data step, or decade step on log), linear/log per axis.
- **Grid overlay** — major/minor, style, colour, alpha, **per-axis** minor
  subdivisions and **per-axis** log minor mantissas; re-derived from the axis
  transform on every calibration change.
- **Snapping** — centerline projection (skeleton + chord refinement along the
  structure-tensor normal), toggleable; radius and dark-threshold adjustable;
  *show curve mask* and *ignore printed grid/axes* diagnostics.
- **Axis locks** — hold **X** for horizontal-only drag (Y taken from the curve),
  **Y** for vertical-only; **Delete** removes the selected marker from plot or table.
- **Table** — five internal columns (`index|x|y|gap|condition`); the `gap` column is
  a fixed-width visual separator and is not editable. CSV output stays 4 columns.
- **Export** — native Save As, default folder `Projects/ImprovedVDMOS/data`,
  written as `index,x,y,condition` in datasheet units.

## Verified state (2026-10-11, calibration update)

- **`ruff check src tests` → clean** for all digitizer code (the pre-existing
  `F401` unused `math` import in `tests/test_fit.py` predates this work).
- **`pytest` → 64 passed, 1 skipped** (the skip is pre-existing). New
  `tests/test_digitizer.py` (26 tests) covers:
  - corner vertex accuracy **exact to <0.05 px** on synthetic printed sheets for
    1/2/3 px line and tick thicknesses, in the display convention;
  - the 90° test rejects a 70° decoy and a shallow curve crossing the frame;
  - lone-tick fallback (frame cropped out) and no-ink raw fallback;
  - `from_points` round-trips linear/log, directional click order, identical-pixel
    rejection, and edge calibration as the exact special case of two points;
  - per-axis minor divisions (2 vs 5) and per-axis log mantissas (2,5 vs 2..9);
  - headless Tk: anchors replaced by clicked pixels, stray clicks while armed
    create no marker, Esc/button-toggle cancel, too-close pair rejected,
    "Apply axes" clears the green marks but keeps the anchors, and a marker
    placed after calibration reads exactly 20 A at the pixel midpoint.
- **Snapping accuracy** (synthetic datasheet PDFs, measured):
  - horizontal traces **1–9 px wide: 0.000 px** centerline error — thickness-independent,
    which was the explicit thick-line requirement;
  - steep anti-aliased traces: **≤0.54 px** perpendicular distance to the true curve,
    mean 0.27 px;
  - X/Y locks preserve the constrained coordinate exactly (0.000 px error).
- **Grid suppression keeps real data**: a 97.5 %-span, 4 px saturation plateau survives,
  while a full-width 2 px gridline and the axis frame are removed.
- **Table renders correctly** after the two layout fixes: frame mapped at 194 px,
  `rowheight = 44 px` on the `vista` theme at this display's scaling, cell text black
  (`foreground #101010` on `fieldbackground #ffffff`), selection white-on-`#2a7fd4`.
- End-to-end headless run: PDF → region → 9 markers → condition edits → X-cell edit
  moves the marker (x = 5.5 exactly) → delete reindexes → CSV written correctly.

## Bugs found and fixed

| # | Symptom | Root cause | Fix |
|---|---|---|---|
| 1 | No curve detected for 1 px traces | 1 px diagonal traces are disconnected under 4-connectivity, so the size filter deleted them | components labelled on a 1 px **dilation**, sizes counted on the raw mask |
| 2 | Long saturation plateaus vanished | grid suppression keyed on dark fraction > 0.8 | strict rule: ≥97 % span **and** ≤3 px thick, user-toggleable |
| 3 | `TypeError: fz_pixmap_xres_set` | PyMuPDF `set_dpi()` rejects the float fit-to-window dpi | `_dpi_int()`; also `import pymupdf` (legacy `fitz` fallback) |
| 4 | Matplotlib grid `UserWarning` | line properties passed for a hidden grid layer | properties supplied only when that layer is visible |
| 5 | Table completely invisible | the container `Frame` was never packed into its parent | `frame.pack(fill="both", expand=True)` |
| 6 | Window taller than the screen | canvas requests `figsize × dpi`, so a big crop pushed the table off-screen | figure capped to 900×470 px from the crop aspect; `_fit_to_screen()` clamps to 94 % of the screen |
| 7 | Row text clipped vertically | `height=` sets row *count*; the style's `rowheight` default is below the font's line height at scaled DPI | `rowheight = linespace × 1.45` (min 24 px), DPI-aware |
| 8 | Cell text white on white, readable only when selected | a **derived** ttk style name (`Digitizer.Treeview`) keeps the layout but loses the theme's per-state colours under `vista` | configure the theme's own `Treeview` / `Treeview.Heading` and `style.map` the selected state |
| 9 | Corner vertex off by 0.75 px on even-width bands | the first corner refinement marched the **binary** mask mixing pixel-centre and pixel-edge sampling conventions | measure band centres from **grayscale threshold-crossing interpolation** in the `imshow` display convention (pixel `j` spans `[j, j+1]`) |
| 10 | A 70° decoy stroke accepted as a tick | strokes were grouped per **column**, so each 1-px column of a diagonal looked vertical | strokes are connected **components** of the mask minus the axis band; direction is the covariance principal axis of the whole stroke |

## Known weak spots

1. **Dark-on-light only.** The mask is a luminance threshold; coloured or inverted
   (white-on-black) datasheets need snapping off, or a threshold change. No hue/colour
   segmentation, so overlapping multi-colour traces are one blob. This applies to the
   corner detector's mask as well.
2. **Log axes cannot show zero.** `AxisTransform` requires strictly positive min/max,
   which is correct mathematically but means a datasheet `I_ds` axis labelled from 0
   must be entered as linear or with a small positive floor.
3. **Suppression and corner constants are hard-coded** (`span_frac = 0.97`,
   `max_thickness = 3`; `WINDOW_RADIUS = 14`, `ANGLE_TOL_DEG = 12`, …); only the
   on/off switches are exposed.
4. **No undo/redo, no autosave, no marker-set persistence.** Closing a pop-up discards
   the markers; there is no re-import of a previously exported CSV.
5. **No rotation/skew correction.** Two points per axis assume the crop is
   axis-parallel; a skewed crop needs a 4-point affine calibration (explicitly
   out of scope by decision).

## Suggested next steps

1. Series grouping: `Condition` → per-`Vgs` blocks or long-format CSV for `vdmos_fit`.
2. Pan/zoom for large crops, and undo/redo for marker edits.
3. Extend `tests/test_digitizer.py` to the snapper's centerline accuracy on
   synthetic traces (the ad-hoc measurements from 2026-10-11 are not yet encoded).
4. Optional: snap calibration clicks to gridline/axis-frame junctions when the
   tick itself is degraded.
