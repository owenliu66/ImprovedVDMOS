"""Curve detection and snapping for the digitizer.

A datasheet trace is typically 2-4 pixels wide.  Snapping to the nearest dark
pixel would bias every extracted point by up to half a line width, so the
snapper works on the trace *centerline*:

1. threshold the crop to a dark-pixel mask (printed gridlines/frames removed);
2. skeletonize the mask (scikit-image when available, distance-transform ridge
   otherwise) to get one point per line cross-section;
3. refine each skeleton point to sub-pixel accuracy by marching along the local
   line normal until the mask exits on both sides and taking the midpoint of
   that chord -- this is thickness independent and degenerates to the pixel
   itself for a 1-pixel line;
4. index the refined cloud with a cKDTree and snap cursor positions to it.
"""

from __future__ import annotations

import numpy as np
from scipy import ndimage
from scipy.spatial import cKDTree

try:  # optional; a distance-transform ridge is used when missing
    from skimage.morphology import skeletonize

    _HAVE_SKIMAGE = True
except ImportError:  # pragma: no cover
    _HAVE_SKIMAGE = False

LOCK_NONE, LOCK_X, LOCK_Y = None, "x", "y"


class CurveSnapper:
    def __init__(
        self,
        image,
        threshold: float = 0.55,
        suppress_printed_grid: bool = True,
        radius: float = 12.0,
        min_component: int = 6,
    ) -> None:
        self.gray = np.asarray(image.convert("L"), dtype=np.float32) / 255.0
        self.height, self.width = self.gray.shape
        self.threshold = float(threshold)
        self.suppress_printed_grid = suppress_printed_grid
        self.radius = float(radius)
        self.min_component = int(min_component)
        self._nearest_index_cache = None
        self.rebuild()

    # ------------------------------------------------------------------ mask
    def rebuild(
        self,
        threshold: float | None = None,
        suppress_printed_grid: bool | None = None,
    ) -> None:
        if threshold is not None:
            self.threshold = float(threshold)
        if suppress_printed_grid is not None:
            self.suppress_printed_grid = bool(suppress_printed_grid)

        mask = self.gray < self.threshold
        if self.suppress_printed_grid:
            mask = self._suppress_straight_lines(mask)
        mask = self._drop_small_components(mask)

        self.mask = mask
        self._nearest_index_cache = None
        # distance from every pixel to the nearest dark pixel (0 on the mask)
        self.dist_dark = ndimage.distance_transform_edt(~mask)
        # distance from every mask pixel to the nearest background pixel,
        # i.e. half the local line thickness on the centerline
        self.half_width = ndimage.distance_transform_edt(mask)
        self.normal_x, self.normal_y = self._normals(mask)

        cloud = self._centerline_cloud(mask)
        self.cloud = cloud
        self.tree = cKDTree(cloud) if len(cloud) else None

    def _suppress_straight_lines(self, mask: np.ndarray) -> np.ndarray:
        """Blank out full-span thin dark rows/columns (plot frames, gridlines).

        The test is deliberately strict: a row must be dark across essentially
        the whole crop (``span_frac``) *and* the dark band must be thin
        (``max_thickness``).  A saturation plateau of an I-V curve can easily
        span 80% of the crop width, so a looser rule would delete real data.
        """
        span_frac, max_thickness = 0.97, 3
        out = mask.copy()
        for axis in (0, 1):
            # fraction of dark pixels along each row (axis=1) / column (axis=0)
            profile = mask.mean(axis=axis)
            candidate = profile > span_frac
            for start, end in self._runs(candidate):
                if end - start + 1 <= max_thickness:
                    if axis == 1:
                        out[start : end + 1, :] = False
                    else:
                        out[:, start : end + 1] = False
        return out

    @staticmethod
    def _runs(flags: np.ndarray) -> list[tuple[int, int]]:
        """Index ranges of consecutive True entries in a 1-D boolean array."""
        idx = np.nonzero(flags)[0]
        if idx.size == 0:
            return []
        breaks = np.nonzero(np.diff(idx) > 1)[0]
        starts = np.concatenate(([idx[0]], idx[breaks + 1]))
        ends = np.concatenate((idx[breaks], [idx[-1]]))
        return list(zip(starts.tolist(), ends.tolist(), strict=True))

    def _drop_small_components(self, mask: np.ndarray) -> np.ndarray:
        """Remove specks without fragmenting steep or anti-aliased traces.

        Components are labelled on a 1-pixel *dilation* of the mask: a steep
        trace rasterizes with row gaps, so labelling the raw mask would split it
        into many tiny pieces that the size test would then delete.  Sizes are
        counted on the original mask and the returned mask is never dilated.
        """
        if not mask.any():
            return mask
        structure = np.ones((3, 3), dtype=int)
        labels, count = ndimage.label(ndimage.binary_dilation(mask, structure=structure),
                                      structure=structure)
        if count <= 1:
            return mask
        sizes = np.bincount(labels[mask], minlength=count + 1)
        keep = sizes >= self.min_component
        keep[0] = False
        return mask & keep[labels]

    # -------------------------------------------------------------- geometry
    def _normals(self, mask: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Unit normal field from the structure tensor of the smoothed mask.

        The dominant eigenvector of the structure tensor points along the
        direction of fastest intensity change, which is the line normal.
        """
        smooth = ndimage.gaussian_filter(mask.astype(np.float32), 1.2)
        gy, gx = np.gradient(smooth)
        jxx = ndimage.gaussian_filter(gx * gx, 2.0)
        jyy = ndimage.gaussian_filter(gy * gy, 2.0)
        jxy = ndimage.gaussian_filter(gx * gy, 2.0)

        trace = jxx + jyy
        det = jxx * jyy - jxy * jxy
        disc = np.maximum(trace * trace / 4.0 - det, 0.0)
        lambda_max = trace / 2.0 + np.sqrt(disc)

        nx = jxy
        ny = lambda_max - jxx
        norm = np.hypot(nx, ny)
        fallback_ok = norm > 1e-6
        # Degenerate (isotropic) regions: fall back to the intensity gradient.
        nx = np.where(fallback_ok, nx, gx)
        ny = np.where(fallback_ok, ny, gy)
        norm = np.hypot(nx, ny)
        np.divide(nx, norm, out=nx, where=norm > 1e-6)
        np.divide(ny, norm, out=ny, where=norm > 1e-6)
        return nx, ny

    def _normal_at(self, u: float, v: float) -> tuple[float, float] | None:
        i = int(round(v))
        j = int(round(u))
        if not (0 <= i < self.height and 0 <= j < self.width):
            return None
        nx = float(self.normal_x[i, j])
        ny = float(self.normal_y[i, j])
        if nx == 0.0 and ny == 0.0:
            return None
        return nx, ny

    def _inside(self, u: float, v: float) -> bool:
        i = int(round(v))
        j = int(round(u))
        return 0 <= i < self.height and 0 <= j < self.width and bool(self.mask[i, j])

    def _march(self, u: float, v: float, du: float, dv: float, limit: float) -> tuple[float, float]:
        """Walk along (du, dv) while inside the mask; return the last inside point."""
        last_u, last_v = u, v
        step = 0.5
        t = step
        while t <= limit:
            nu, nv = u + du * t, v + dv * t
            if not self._inside(nu, nv):
                break
            last_u, last_v = nu, nv
            t += step
        return last_u, last_v

    def refine(self, u: float, v: float, iterations: int = 2) -> tuple[float, float]:
        """Project a point onto the local line centerline (sub-pixel)."""
        for _ in range(iterations):
            normal = self._normal_at(u, v)
            if normal is None:
                break
            reach = 4.0 * max(self.half_width_at(u, v), 0.5) + 3.0
            b1 = self._march(u, v, normal[0], normal[1], reach)
            b2 = self._march(u, v, -normal[0], -normal[1], reach)
            nu, nv = 0.5 * (b1[0] + b2[0]), 0.5 * (b1[1] + b2[1])
            if not self._inside(nu, nv):
                break
            if abs(nu - u) < 0.2 and abs(nv - v) < 0.2:
                u, v = nu, nv
                break
            u, v = nu, nv
        return u, v

    def half_width_at(self, u: float, v: float) -> float:
        i = int(round(v))
        j = int(round(u))
        if not (0 <= i < self.height and 0 <= j < self.width):
            return 0.0
        return float(self.half_width[i, j])

    def _centerline_cloud(self, mask: np.ndarray) -> np.ndarray:
        if not mask.any():
            return np.empty((0, 2), dtype=np.float64)

        if _HAVE_SKIMAGE:
            ridge = skeletonize(mask)
        else:  # distance-transform ridge fallback
            ridge = mask & (self.half_width >= ndimage.maximum_filter(self.half_width, 3))
        rows, cols = np.nonzero(ridge)
        if rows.size == 0:
            return np.empty((0, 2), dtype=np.float64)

        refined = np.empty((rows.size, 2), dtype=np.float64)
        kept = []
        for k, (r, c) in enumerate(zip(rows, cols, strict=False)):
            u, v = self.refine(float(c), float(r))
            if abs(u - c) > 2.5 or abs(v - r) > 2.5:
                continue  # junctions/endpoints where the chord is meaningless
            refined[k] = (u, v)
            kept.append(k)
        return refined[kept] if kept else np.stack([cols, rows], axis=1).astype(np.float64)

    # ------------------------------------------------------------------ snap
    def _nearest_dark_pixel(self, u: float, v: float) -> tuple[float, float] | None:
        if self.mask.any():
            if self._nearest_index_cache is None:
                self._nearest_index_cache = ndimage.distance_transform_edt(
                    self.mask, return_indices=True
                )
            _, (ri, ci) = self._nearest_index_cache
            i = int(np.clip(round(v), 0, self.height - 1))
            j = int(np.clip(round(u), 0, self.width - 1))
            return float(ci[i, j]), float(ri[i, j])
        return None

    def snap(
        self,
        u: float,
        v: float,
        lock: str | None = None,
        radius: float | None = None,
    ) -> tuple[float, float] | None:
        """Return the snapped (u, v) for a cursor position, or None.

        ``lock='x'`` means the X key is held (horizontal-only motion): u is kept
        and v is taken from the curve at that column.  ``lock='y'`` is symmetric.
        """
        radius = self.radius if radius is None else float(radius)
        if self.tree is None:
            return None

        point = None
        idx = self.tree.query((u, v), distance_upper_bound=radius)
        if np.isfinite(idx[0]):
            point = self.refine(*self.cloud[idx[1]])
        else:
            dist = self.dist_dark_at(u, v)
            if dist <= radius:
                dark = self._nearest_dark_pixel(u, v)
                if dark is not None:
                    point = self.refine(*dark)
        if point is None:
            return None

        if lock == LOCK_X:
            candidates = self._cloud_at_column(u, v, radius)
            return (u, candidates) if candidates is not None else (u, point[1])
        if lock == LOCK_Y:
            candidates = self._cloud_at_row(v, u, radius)
            return (candidates, v) if candidates is not None else (point[0], v)
        return point

    def _cloud_at_column(self, u: float, v: float, radius: float) -> float | None:
        """Curve value v(u): average centerline points in a +-1.5 px column band."""
        if self.tree is None:
            return None
        hits = self.tree.query_ball_point((u, v), radius)
        if not hits:
            return None
        pts = self.cloud[hits]
        band = pts[np.abs(pts[:, 0] - u) <= 1.5]
        if band.size == 0:
            return None
        if len(band) > 1:
            band = band[np.argsort(np.abs(band[:, 1] - v))]
            near = band[np.abs(band[:, 1] - band[0, 1]) <= 3.0]
            return float(np.mean(near[:, 1]))
        return float(band[0, 1])

    def _cloud_at_row(self, v: float, u: float, radius: float) -> float | None:
        if self.tree is None:
            return None
        hits = self.tree.query_ball_point((u, v), radius)
        if not hits:
            return None
        pts = self.cloud[hits]
        band = pts[np.abs(pts[:, 1] - v) <= 1.5]
        if band.size == 0:
            return None
        if len(band) > 1:
            band = band[np.argsort(np.abs(band[:, 0] - u))]
            near = band[np.abs(band[:, 0] - band[0, 0]) <= 3.0]
            return float(np.mean(near[:, 0]))
        return float(band[0, 0])

    def dist_dark_at(self, u: float, v: float) -> float:
        i = int(np.clip(round(v), 0, self.height - 1))
        j = int(np.clip(round(u), 0, self.width - 1))
        return float(self.dist_dark[i, j])
