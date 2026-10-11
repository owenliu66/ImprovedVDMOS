"""Editable marker table: | Index | X | Y | Condition |."""

from __future__ import annotations

import tkinter as tk
from tkinter import font as tkfont
from tkinter import ttk

COLUMNS = ("index", "x", "y", "gap", "condition")
HEADINGS = {"index": "#", "x": "X", "y": "Y", "gap": "", "condition": "Condition"}
COLUMN_IDS = {"#1": "index", "#2": "x", "#3": "y", "#4": "gap", "#5": "condition"}
NON_EDITABLE = ("index", "gap")
NUMERIC = ("x", "y")


class MarkerTable:
    """Treeview with in-place cell editing, kept in sync with the plot markers.

    ``host`` must expose ``row_text(index) -> tuple[str, str, str, str]``,
    ``on_cell_edited(index, column, text)``, ``on_row_selected(index)`` and
    ``on_delete_requested(index)``.
    """

    def __init__(self, parent: tk.Widget, host) -> None:
        self.host = host
        self._editor: ttk.Entry | None = None
        self._updating = False

        frame = ttk.Frame(parent)
        self.frame = frame
        # The frame itself must be packed, otherwise the whole table is invisible.
        frame.pack(fill="both", expand=True)
        self.tree = ttk.Treeview(
            frame, columns=COLUMNS, show="headings", height=8, selectmode="browse"
        )
        self._apply_row_metrics()
        widths = {"index": 46, "x": 130, "y": 130, "gap": 16, "condition": 200}
        for column in COLUMNS:
            self.tree.heading(column, text=HEADINGS[column])
            anchor = "e" if column in NUMERIC or column == "index" else "w"
            stretch = column != "gap"  # the gap column is a fixed visual separator
            self.tree.column(column, width=widths[column], anchor=anchor, stretch=stretch)

        bar = ttk.Scrollbar(frame, orient="vertical", command=self.tree.yview)
        self.tree.configure(yscrollcommand=bar.set)
        self.tree.pack(side="left", fill="both", expand=True)
        bar.pack(side="right", fill="y")

        self.tree.bind("<Double-1>", self._begin_edit)
        self.tree.bind("<<TreeviewSelect>>", self._on_select)
        self.tree.bind("<Delete>", self._on_delete)
        self.tree.bind("<Escape>", lambda _e: self._cancel_edit())

    def _apply_row_metrics(self) -> None:
        """Give rows enough height for the font.

        ``height=`` only sets the number of visible rows; the pixel row height
        comes from the style's ``rowheight``, whose default is smaller than the
        line height of the UI font at scaled DPI -- which clips cell text.

        The theme's own ``Treeview`` style is configured rather than a derived
        custom name: a derived style keeps the layout but loses the theme's
        per-state text colours, which renders white text on white rows.
        """
        base = tkfont.nametofont("TkDefaultFont")
        family, size = base.actual("family"), base.actual("size")
        linespace = base.metrics("linespace")
        rowheight = max(24, int(round(linespace * 1.45)))

        style = ttk.Style(self.tree)
        style.configure(
            "Treeview",
            rowheight=rowheight,
            font=(family, size),
            background="#ffffff",
            foreground="#101010",
            fieldbackground="#ffffff",
            borderwidth=1,
        )
        style.configure(
            "Treeview.Heading", font=(family, size, "bold"), rowheight=rowheight
        )
        style.map(
            "Treeview",
            background=[("selected", "#2a7fd4")],
            foreground=[("selected", "#ffffff")],
        )
        self.rowheight = rowheight

    def _values(self, index: int) -> tuple[str, str, str, str, str]:
        """Row text with the empty separator column inserted before Condition."""
        index_text, x_text, y_text, condition = self.host.row_text(index)
        return (index_text, x_text, y_text, "", condition)

    # ------------------------------------------------------------- refresh
    def refresh(self) -> None:
        self._updating = True
        try:
            count = self.host.marker_count()
            existing = self.tree.get_children()
            for item in existing:
                self.tree.delete(item)
            for index in range(count):
                self.tree.insert("", "end", iid=str(index), values=self._values(index))
        finally:
            self._updating = False

    def update_row(self, index: int) -> None:
        item = str(index)
        if self.tree.exists(item):
            self.tree.item(item, values=self._values(index))

    def select_row(self, index: int | None) -> None:
        if self._updating:
            return
        self.tree.selection_set()
        if index is None:
            return
        item = str(index)
        if self.tree.exists(item):
            self.tree.selection_set(item)
            self.tree.see(item)

    # -------------------------------------------------------------- events
    def _on_select(self, _event) -> None:
        if self._updating:
            return
        selection = self.tree.selection()
        self.host.on_row_selected(int(selection[0]) if selection else None)

    def _on_delete(self, _event) -> None:
        selection = self.tree.selection()
        if selection:
            self.host.on_delete_requested(int(selection[0]))

    def _begin_edit(self, event) -> None:
        if self._editor is not None:
            self._commit_edit()
        item = self.tree.identify_row(event.y)
        column_id = self.tree.identify_column(event.x)
        if not item or column_id not in COLUMN_IDS:
            return
        column = COLUMN_IDS[column_id]
        if column in NON_EDITABLE:
            return
        bbox = self.tree.bbox(item, column)
        if not bbox:
            return
        x, y, width, height = bbox
        entry = ttk.Entry(self.tree, font=("TkDefaultFont",))
        entry.place(x=x, y=y, width=width, height=max(height, self.rowheight))
        entry.insert(0, str(self.tree.set(item, column)))
        entry.select_range(0, tk.END)
        entry.focus_set()
        entry.bind("<Return>", lambda _e: self._commit_edit())
        entry.bind("<FocusOut>", lambda _e: self._commit_edit())
        entry.bind("<Escape>", lambda _e: self._cancel_edit())
        self._editor = entry
        self._editor_target = (int(item), column)

    def _cancel_edit(self) -> None:
        if self._editor is not None:
            self._editor.destroy()
            self._editor = None

    def _commit_edit(self) -> None:
        if self._editor is None:
            return
        entry, (index, column) = self._editor, self._editor_target
        text = entry.get()
        entry.destroy()
        self._editor = None
        self.host.on_cell_edited(index, column, text)
