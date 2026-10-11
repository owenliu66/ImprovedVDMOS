"""PDF page and region rendering via PyMuPDF."""

from __future__ import annotations

import io
from pathlib import Path

from PIL import Image

try:  # PyMuPDF: ``pymupdf`` is the current name, ``fitz`` the legacy alias
    import pymupdf as fitz
except ImportError:  # pragma: no cover
    try:
        import fitz
    except ImportError:
        fitz = None

POINTS_PER_INCH = 72.0


def _dpi_int(dpi: float) -> int:
    """PyMuPDF's set_dpi() only accepts integers."""
    return max(1, int(round(float(dpi))))


class PdfSource:
    def __init__(self, path: str | Path) -> None:
        if fitz is None:
            raise RuntimeError(
                "PyMuPDF is required to read PDFs. Install it with: pip install pymupdf"
            )
        self.path = Path(path)
        self.doc = fitz.open(self.path)
        if self.doc.page_count == 0:
            raise RuntimeError(f"{self.path} contains no pages")

    @property
    def page_count(self) -> int:
        return self.doc.page_count

    def page_size_pt(self, page: int) -> tuple[float, float]:
        rect = self.doc[page].rect
        return float(rect.width), float(rect.height)

    def render_page(self, page: int, dpi: float) -> Image.Image:
        pix = self.doc[page].get_pixmap(dpi=_dpi_int(dpi))
        return Image.open(io.BytesIO(pix.tobytes("png"))).convert("RGB")

    def render_region(
        self, page: int, rect_pt: tuple[float, float, float, float], dpi: float
    ) -> Image.Image:
        """Render a sub-rectangle (in PDF points) at high resolution."""
        clip = fitz.Rect(*rect_pt)
        pix = self.doc[page].get_pixmap(dpi=_dpi_int(dpi), clip=clip)
        return Image.open(io.BytesIO(pix.tobytes("png"))).convert("RGB")

    def close(self) -> None:
        try:
            self.doc.close()
        except Exception:
            pass
