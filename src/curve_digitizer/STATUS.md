# Curve Digitizer — Status

*Snapshot date: 2026-10-11. Companion to [README.md](README.md); the project-level
status lives in [../../STATUS.md](../../STATUS.md).*

## What exists

```
src/curve_digitizer/        1,698 lines of Python, 11 modules
├── __init__.py       (13)  package exports, `python -m curve_digitizer`
├── __main__.py        (4)  entry point
├── app.py           (205)  main Tk window: PDF load, page nav, crop DPI, rubber-band region select
├── pdf_source.py     (61)  PyMuPDF page/region rendering (region re-rendered, not cropped)
├── transform.py     (162)  AxisTransform: pixel <-> data, linear/log, tick generation
├── grid.py          (112)  adjustable Cartesian grid overlay derived from the transform
├── snapping.py      (325)  dark-pixel mask -> skeleton -> sub-pixel centerline cloud -> KDTree
├── markers.py       (224)  add/drag/select/delete, X/Y key axis-lock, marker artists
├── digitizer_window.py (380) pop-up: plot, axis panel, grid panel, snap controls, Export
├── table.py         (173)  Index|X|Y|Condition Treeview, in-place cell editing, 2-way sync
└── export.py         (39)  CSV via native Save As, default `Projects/ImprovedVDMOS/data`
```

Not yet committed — `git status` shows `?? src/curve_digitizer/`. Runtime dependencies
are declared in `pyproject.toml` as the `digitizer` extra
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
- **Axes** — X/Y min and max entered at the crop edges, tick interval (data step,
  or decade step on log), linear/log per axis.
- **Grid overlay** — major/minor, style, colour, alpha, linear minor subdivisions,
  log minor mantissas; re-derived from the axis transform on every calibration change.
- **Snapping** — centerline projection (skeleton + chord refinement along the
  structure-tensor normal), toggleable; radius and dark-threshold adjustable;
  *show curve mask* and *ignore printed grid/axes* diagnostics.
- **Axis locks** — hold **X** for horizontal-only drag (Y taken from the curve),
  **Y** for vertical-only; **Delete** removes the selected marker from plot or table.
- **Table** — five internal columns (`index|x|y|gap|condition`); the `gap` column is
  a fixed-width visual separator and is not editable. CSV output stays 4 columns.
- **Export** — native Save As, default folder `Projects/ImprovedVDMOS/data`,
  written as `index,x,y,condition` in datasheet units.

## Verified state (2026-10-11)

- **`ruff check src/curve_digitizer` → clean** (the `I001` import-order error
  introduced with the `tkfont` import is fixed).
- All modules compile; `import curve_digitizer` resolves with `PYTHONPATH=src`.
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

## Known weak spots

1. **Calibration is edge-based.** X/Y min-max are assumed to sit exactly at the crop
   edges, so the box must be drawn on the axis extremes. A box that includes the tick
   labels or axis frame is mis-calibrated. There is no click-two-points-to-calibrate yet.
2. **Dark-on-light only.** The mask is a luminance threshold; coloured or inverted
   (white-on-black) datasheets need snapping off, or a threshold change. No hue/colour
   segmentation, so overlapping multi-colour traces are one blob.
3. **Log axes cannot show zero.** `AxisTransform` requires strictly positive min/max,
   which is correct mathematically but means a datasheet `I_ds` axis labelled from 0
   must be entered as linear or with a small positive floor.
4. **Suppression constants are hard-coded** (`span_frac = 0.97`, `max_thickness = 3`);
   only the on/off switch is exposed.
5. **No undo/redo, no autosave, no marker-set persistence.** Closing a pop-up discards
   the markers; there is no re-import of a previously exported CSV.
6. **No committed tests.** All verification above was done with ad-hoc headless scripts
   (Tk driven without `mainloop`, `ImageGrab` screenshots); nothing lives in `tests/`.

## Suggested next steps

1. Click-two-ticks calibration (pick the two points that define each axis) instead of
   assuming the crop edges are the axis extremes.
2. `tests/test_digitizer.py` covering `AxisTransform` round-trips, tick/minor generation
   and the snapper's centerline accuracy on synthetic traces.
3. Series grouping: `Condition` → per-`Vgs` blocks or long-format CSV for `vdmos_fit`.
4. Pan/zoom for large crops, and undo/redo for marker edits.
