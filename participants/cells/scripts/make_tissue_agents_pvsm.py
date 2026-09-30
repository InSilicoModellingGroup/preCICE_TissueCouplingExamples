#!/usr/bin/env python3
"""Build a ParaView state for tissue_41 agent VTUs (run with pvbatch/pvpython).

Usage:
  pvbatch scripts/make_tissue_agents_pvsm.py
"""
from __future__ import annotations

import glob
import os
import re
import sys

from paraview.simple import (
    XMLPartitionedUnstructuredGridReader,
    Glyph,
    GetActiveViewOrCreate,
    GetColorTransferFunction,
    GetOpacityTransferFunction,
    Show,
    ColorBy,
    Hide,
    Render,
    ResetCamera,
    SaveState,
    Delete,
)


def main() -> None:
    root = os.path.abspath(
        os.path.join(os.path.dirname(__file__), "..", "output", "tissue_41_diffusion", "cells")
    )
    files = sorted(
        glob.glob(os.path.join(root, "MyCell-*.pvtu")),
        key=lambda p: int(re.search(r"MyCell-(\d+)\.pvtu$", os.path.basename(p)).group(1)),
    )
    if not files:
        sys.exit(f"No MyCell-*.pvtu under {root}")

    reader = XMLPartitionedUnstructuredGridReader(FileName=files)
    reader.PointArrayStatus = ["cell_color_", "diameter_", "h2o2_", "no2_", "track_id_"]

    view = GetActiveViewOrCreate("RenderView")
    view.Background = [0.08, 0.08, 0.1]

    # Points alone are invisible at mm-scale — glyph as spheres.
    glyph = Glyph(Input=reader, GlyphType="Sphere")
    glyph.OrientationArray = ["POINTS", "No orientation array"]
    glyph.ScaleArray = ["POINTS", "diameter_"]
    glyph.ScaleFactor = 1.0
    glyph.GlyphMode = "All Points"
    try:
        glyph.GlyphType.Radius = 0.5  # unit sphere; scaled by diameter_
    except Exception:
        pass

    display = Show(glyph, view)
    ColorBy(display, ("POINTS", "h2o2_"))
    display.RescaleTransferFunctionToDataRange(True, False)
    lut = GetColorTransferFunction("h2o2_")
    lut.RescaleTransferFunction(0.0, 100.0)
    GetOpacityTransferFunction("h2o2_").RescaleTransferFunction(0.0, 100.0)
    display.SetScalarBarVisibility(view, True)

    Hide(reader, view)
    ResetCamera(view)
    # Zoom to ~spheroid near top of 2 mm cylinder
    cam = view.GetActiveCamera()
    cam.SetFocalPoint(0.0, 0.0, 0.00177)
    cam.SetPosition(0.0, -0.0015, 0.00177)
    cam.SetViewUp(0.0, 0.0, 1.0)
    view.CameraParallelScale = 0.00035
    Render()

    out = os.path.join(root, "cells.pvsm")
    SaveState(out)
    print(f"Wrote {out}  ({len(files)} timesteps)")
    print(f"Also open: {os.path.join(root, 'MyCell.pvd')}")


if __name__ == "__main__":
    main()
