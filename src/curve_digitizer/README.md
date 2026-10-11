# curve_digitizer

Interactive digitizer for MOSFET datasheet plots (typically `I_ds` vs `V_ds` at several
`V_gs`). Loads a PDF, lets you box a plot region, then extracts marker points against a
calibrated Cartesian overlay and exports CSV.

## Run

```bash
python -m curve_digitizer          # from src/, or with src/ on PYTHONPATH
```

Runtime dependencies beyond the project's NumPy/SciPy/Matplotlib: **PyMuPDF** (`pymupdf`)
for PDF rendering and **scikit-image** for curve skeletonization (optional; a
distance-transform ridge is used if it is missing).

## Workflow

1. **Load PDF** → page navigation, crop DPI (default 300). Drag a rectangle over the plot;
   a digitizer pop-up opens for that region (re-rendered at crop DPI, so it stays sharp).
2. **Calibrate the axes** — two options, and they mix freely per axis:
   - **Two-point (recommended):** type the data values of two printed ticks into
     `X at pt 1 / X at pt 2` (and Y), press **Calibrate X**, then click those two
     ticks on the plot. Each click snaps to the **90-degree corner** where the
     tick stroke crosses the axis line (sub-pixel, grayscale-interpolated), so
     the anchor lands on the printed tick centre rather than on your cursor.
     The status line reports the snap distance or "no tick corner found - raw
     click used". **Esc** or pressing the same button again cancels; the
     *snap clicks to tick corners* checkbox turns snapping off. Only the
     coordinate along the calibrated axis matters, so you can click anywhere
     along the tick's alignment line. Click order is directional: pt 1 maps to
     the first value, pt 2 to the second, so right-to-left clicks simply
     reverse the axis. The accepted pair stays drawn in green until you press
     **Apply axes** (which clears the marks but keeps the calibration).
   - **Edge-based (original):** leave the anchors alone and enter the values at
     the crop edges, as before — the default anchors are exactly the crop edges.
     Re-calibrating on the two frame corners restores edge calibration.
   Then set the tick interval and `linear`/`log` per axis and press **Apply
   axes**. For a log axis the tick interval is the decade step; two points
   uniquely fix the decade scale and offset. Markers are pixel-anchored, so
   recalibrating re-values the table and CSV without moving any marker.
3. **Grid panel** — major/minor grid, line style, color, alpha, **per-axis** minor
   subdivisions (`X minor div` / `Y minor div`) and **per-axis** log minor mantissas
   (`X log minor` / `Y log minor`). The grid is derived from the axis transform, so it
   re-registers whenever the calibration changes; lining it up with the datasheet's own
   printed grid is the visual check that the calibration is right.
4. **Markers** — click the plot to add one, drag to move. Snapping projects onto the curve
   *centerline*; toggle **Snap to curve** off if datasheet coloring confuses detection, and
   use **show curve mask** to inspect what the detector sees. **ignore printed grid/axes**
   drops full-span thin dark rows/columns (plot frames, gridlines).
   - Hold **X** while dragging: horizontal-only motion (Y taken from the curve).
   - Hold **Y** while dragging: vertical-only motion.
   - Select a marker and press **Delete** to remove it.
5. **Table** `| # | X | Y | Condition |` — double-click a cell to edit. Editing X or Y moves
   the marker (no snapping: typed values are authoritative); Condition is free text such as
   `Vgs=4.5V`.
6. **Export CSV** — native Save As dialog, default folder `Projects/ImprovedVDMOS/data`,
   written as `index,x,y,condition` in datasheet units.

## Coordinate model

The image is always drawn in **crop-pixel space** (`imshow(extent=(0, W, H, 0))`); data
values appear only as tick labels, table entries and CSV output. A log axis printed on paper
already has compressed decade spacing in pixels, so warping the image onto a log-scaled
matplotlib axis would distort it — instead `transform.AxisTransform` maps pixels → data with
a linear or log law, and ticks/grid are placed at the pixel positions of nice data values.

## Snapping accuracy

Skeletonize the dark-pixel mask, then refine each ridge point by marching along the local
line normal (from the structure tensor) to both mask boundaries and taking the chord
midpoint. This is thickness-independent: measured centerline error is 0.000 px for
horizontal lines 1-9 px wide, and ≤0.5 px perpendicular distance on steep anti-aliased
curves.


## Calibration corner detection

`corners.py` finds the tick/axis vertex in a 28 px window around a calibration click:
the axis line is the row with the longest horizontal dark run; dark components outside
that band which touch it are candidate ticks; a component is accepted only if its
principal direction is within 12 degrees of perpendicular (rejecting shallow curve
crossings, arrowheads and text). Band centres come from grayscale threshold-crossing
interpolation, not from the binary mask, so anti-aliased printed lines place between
pixel centres. Positions follow the display convention of `imshow(extent=(0, W, H, 0))`:
pixel `j` spans `[j, j+1]`, so a band at columns 120-121 centres on 121.0. Fallbacks:
a lone perpendicular stroke (frame cropped out) snaps the along-axis coordinate only;
a click with no ink nearby stays raw. The corner mask is deliberately separate from
`CurveSnapper`'s — the snapper deletes frames and small components that corner
detection needs — but both share the dark-threshold slider.
