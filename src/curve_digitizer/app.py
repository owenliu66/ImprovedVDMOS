"""Main window: load a PDF datasheet and rubber-band a plot region."""

from __future__ import annotations

import tkinter as tk
from tkinter import filedialog, messagebox, ttk

from PIL import ImageTk

from .digitizer_window import DigitizerWindow
from .pdf_source import POINTS_PER_INCH, PdfSource

MIN_REGION_PX = 15


class DigitizerApp:
    def __init__(self, root: tk.Tk) -> None:
        self.root = root
        root.title("Curve Digitizer - load a PDF and drag a box over the plot")
        root.geometry("1100x820")

        self.source: PdfSource | None = None
        self.page = 0
        self.dpi = POINTS_PER_INCH
        self.offset_x = 0.0
        self.offset_y = 0.0
        self.tk_image: ImageTk.PhotoImage | None = None
        self.region_count = 0
        self.windows: list[DigitizerWindow] = []
        self.drag_origin: tuple[int, int] | None = None
        self.drag_item: int | None = None
        self._configure_job: str | None = None

        self._build_toolbar()
        self.canvas = tk.Canvas(root, bg="#2b2b2b", highlightthickness=0)
        self.canvas.pack(fill="both", expand=True)
        self.canvas.bind("<Configure>", self._on_configure)
        self.canvas.bind("<ButtonPress-1>", self._on_press)
        self.canvas.bind("<B1-Motion>", self._on_motion)
        self.canvas.bind("<ButtonRelease-1>", self._on_release)

        self.status = ttk.Label(root, text="No PDF loaded.", padding=(6, 3))
        self.status.pack(side="bottom", fill="x")

    # ------------------------------------------------------------- toolbar
    def _build_toolbar(self) -> None:
        bar = ttk.Frame(self.root, padding=(6, 4))
        bar.pack(side="top", fill="x")

        ttk.Button(bar, text="Load PDF...", command=self.load_pdf).pack(side="left")
        ttk.Button(bar, text="<", width=3, command=lambda: self.goto_page(self.page - 1)).pack(
            side="left", padx=(12, 0)
        )
        self.page_var = tk.IntVar(value=1)
        self.page_entry = ttk.Entry(bar, textvariable=self.page_var, width=5)
        self.page_entry.pack(side="left")
        self.page_entry.bind("<Return>", lambda _e: self.goto_page(self.page_var.get() - 1))
        ttk.Label(bar, text="/ 0").pack(side="left", padx=(2, 0))
        self.page_total = ttk.Label(bar, text="")
        self.page_total.pack(side="left")
        ttk.Button(bar, text=">", width=3, command=lambda: self.goto_page(self.page + 1)).pack(
            side="left", padx=(0, 12)
        )

        ttk.Label(bar, text="crop DPI").pack(side="left")
        self.crop_dpi = tk.IntVar(value=300)
        ttk.Spinbox(bar, from_=100, to=800, increment=50, textvariable=self.crop_dpi, width=6).pack(
            side="left", padx=(2, 12)
        )
        ttk.Label(
            bar, text="Drag a rectangle over the plot; a digitizer window opens."
        ).pack(side="left")

    # ---------------------------------------------------------------- load
    def load_pdf(self) -> None:
        path = filedialog.askopenfilename(
            title="Open datasheet PDF",
            filetypes=[("PDF files", "*.pdf"), ("All files", "*.*")],
        )
        if not path:
            return
        try:
            source = PdfSource(path)
        except Exception as exc:  # pragma: no cover - user-facing guard
            messagebox.showerror("Cannot open PDF", str(exc), parent=self.root)
            return
        if self.source is not None:
            self.source.close()
        self.source = source
        self.page = 0
        self.page_var.set(1)
        self.page_total.config(text=f"/ {source.page_count}")
        self._render_page()

    def goto_page(self, index: int) -> None:
        if self.source is None:
            return
        self.page = max(0, min(index, self.source.page_count - 1))
        self.page_var.set(self.page + 1)
        self._render_page()

    # ------------------------------------------------------------ rendering
    def _on_configure(self, _event) -> None:
        if self.source is None:
            return
        if self._configure_job:
            self.root.after_cancel(self._configure_job)
        self._configure_job = self.root.after(200, self._render_page)

    def _render_page(self) -> None:
        self._configure_job = None
        if self.source is None:
            return
        canvas_w = max(self.canvas.winfo_width(), 200)
        canvas_h = max(self.canvas.winfo_height(), 200)
        page_w_pt, page_h_pt = self.source.page_size_pt(self.page)
        fit = min(canvas_w / page_w_pt, canvas_h / page_h_pt)
        self.dpi = max(40.0, min(fit * POINTS_PER_INCH, 220.0))
        image = self.source.render_page(self.page, self.dpi)
        self.offset_x = max((canvas_w - image.width) / 2.0, 0.0)
        self.offset_y = max((canvas_h - image.height) / 2.0, 0.0)
        self.tk_image = ImageTk.PhotoImage(image)
        self.canvas.delete("all")
        self.canvas.create_image(self.offset_x, self.offset_y, anchor="nw", image=self.tk_image)
        self.status.config(
            text=(
                f"{self.source.path.name} - page {self.page + 1}/{self.source.page_count} "
                f"at {self.dpi:.0f} dpi"
            )
        )

    # ------------------------------------------------------------- selection
    def _canvas_point(self, event) -> tuple[float, float]:
        return event.x - self.offset_x, event.y - self.offset_y

    def _on_press(self, event) -> None:
        if self.source is None:
            return
        self.drag_origin = (event.x, event.y)
        self.drag_item = self.canvas.create_rectangle(
            event.x, event.y, event.x, event.y,
            outline="#ff4040", width=2, dash=(6, 3),
        )

    def _on_motion(self, event) -> None:
        if self.drag_item is not None:
            self.canvas.coords(self.drag_item, *self.drag_origin, event.x, event.y)

    def _on_release(self, event) -> None:
        if self.source is None or self.drag_origin is None:
            return
        x0, y0 = self.drag_origin
        x1, y1 = event.x, event.y
        self.drag_origin = None
        self.canvas.delete(self.drag_item)
        self.drag_item = None
        left, right = min(x0, x1), max(x0, x1)
        top, bottom = min(y0, y1), max(y0, y1)
        if right - left < MIN_REGION_PX or bottom - top < MIN_REGION_PX:
            self.status.config(text="Region too small - drag a box around the plot.")
            return

        scale = POINTS_PER_INCH / self.dpi  # canvas px -> PDF points
        rect_pt = (
            (left - self.offset_x) * scale,
            (top - self.offset_y) * scale,
            (right - self.offset_x) * scale,
            (bottom - self.offset_y) * scale,
        )
        page_w_pt, page_h_pt = self.source.page_size_pt(self.page)
        rect_pt = (
            max(rect_pt[0], 0.0),
            max(rect_pt[1], 0.0),
            min(rect_pt[2], page_w_pt),
            min(rect_pt[3], page_h_pt),
        )
        try:
            crop = self.source.render_region(self.page, rect_pt, float(self.crop_dpi.get()))
        except Exception as exc:  # pragma: no cover - user-facing guard
            messagebox.showerror("Cannot render region", str(exc), parent=self.root)
            return
        if crop.width < MIN_REGION_PX or crop.height < MIN_REGION_PX:
            self.status.config(text="Region too small at this crop DPI.")
            return

        self.region_count += 1
        window = DigitizerWindow(self.root, crop, self.page, self.region_count)
        self.windows.append(window)
        window.bind(
            "<Destroy>",
            lambda _e, w=window: self.windows.remove(w) if w in self.windows else None,
        )
        self.status.config(
            text=f"Opened region {self.region_count} ({crop.width}x{crop.height} px)."
        )


def main() -> None:
    root = tk.Tk()
    DigitizerApp(root)
    root.mainloop()


if __name__ == "__main__":  # pragma: no cover
    main()
