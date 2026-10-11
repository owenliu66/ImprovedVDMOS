"""Interactive curve digitizer for MOSFET datasheet PDFs.

Run with::

    python -m curve_digitizer        # with src/ on PYTHONPATH

Workflow: load a PDF, drag a rectangle over a plot, then calibrate the pop-up
axes, place/drag markers and export them as ``index,x,y,condition`` CSV.
"""

from .app import DigitizerApp, main

__all__ = ["DigitizerApp", "main"]
