#!/usr/bin/env python3
"""Rebuild MyCell-{iter}.pvtu + MyCell.pvd from BioDynaMo piece VTUs.

Needed when a coupled run is killed before BioDynaMo finalizes ParaView metadata
(no .pvtu / .pvsm were written).
"""
from __future__ import annotations

import argparse
import re
from collections import defaultdict
from pathlib import Path


def rebuild(cells_dir: Path, dt: float = 0.001) -> None:
    pieces = sorted(cells_dir.glob("MyCell-*_*.vtu"))
    by_iter: dict[int, list[tuple[int, str]]] = defaultdict(list)
    for p in pieces:
        m = re.match(r"MyCell-(\d+)_(\d+)\.vtu$", p.name)
        if not m:
            continue
        by_iter[int(m.group(1))].append((int(m.group(2)), p.name))
    if not by_iter:
        raise SystemExit(f"No MyCell-*_*.vtu pieces in {cells_dir}")

    arr_lines = [
        '      <PDataArray type="Float64" Name="cell_color_" NumberOfComponents="3"/>',
        '      <PDataArray type="Float64" Name="diameter_"/>',
        '      <PDataArray type="Float64" Name="h2o2_"/>',
        '      <PDataArray type="Float64" Name="no2_"/>',
        '      <PDataArray type="UInt64" Name="track_id_"/>',
    ]
    pvd = [
        '<?xml version="1.0"?>',
        '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">',
        "  <Collection>",
    ]
    for it in sorted(by_iter):
        ranks = sorted(by_iter[it], key=lambda x: x[0])
        body = [
            '<VTKFile type="PUnstructuredGrid" version="0.1" byte_order="LittleEndian" '
            'header_type="UInt32" compressor="vtkZLibDataCompressor">',
            '  <PUnstructuredGrid GhostLevel="0">',
            "    <PPointData>",
            *arr_lines,
            "    </PPointData>",
            "    <PPoints>",
            '      <PDataArray type="Float64" Name="position_" NumberOfComponents="3"/>',
            "    </PPoints>",
        ]
        for _, name in ranks:
            body.append(f'    <Piece Source="{name}"/>')
        body += ["  </PUnstructuredGrid>", "</VTKFile>", ""]
        (cells_dir / f"MyCell-{it}.pvtu").write_text("\n".join(body))
        pvd.append(
            f'    <DataSet timestep="{it * dt:.6g}" group="" part="0" '
            f'file="MyCell-{it}.pvtu"/>'
        )
    pvd += ["  </Collection>", "</VTKFile>", ""]
    text = "\n".join(pvd)
    (cells_dir / "MyCell.pvd").write_text(text)
    (cells_dir / "cells.pvd").write_text(text)
    print(
        f"Wrote {len(by_iter)} pvtu + MyCell.pvd under {cells_dir}\n"
        f"  t = 0 .. {max(by_iter) * dt:.3f} s"
    )


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "cells_dir",
        nargs="?",
        type=Path,
        default=Path(__file__).resolve().parents[1]
        / "output"
        / "tissue_41_diffusion"
        / "cells",
    )
    p.add_argument("--dt", type=float, default=0.001)
    args = p.parse_args()
    rebuild(args.cells_dir, args.dt)


if __name__ == "__main__":
    main()
