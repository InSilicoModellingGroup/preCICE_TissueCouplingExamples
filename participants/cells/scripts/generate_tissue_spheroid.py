#!/usr/bin/env python3
"""Compact HCP spheroid agent positions for the tissue cylinder cases.

Default: 5000 cells, diameter 20 µm, centered in the cylinder.
Domain: cylinder R=0.5 mm, H=1 mm, axis +z, top at z=H.
"""
from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def generate(
    n: int,
    diameter: float,
    out: Path,
    r_cyl: float = 0.0005,
    height: float = 0.001,
    spacing_factor: float = 1.02,
    packing_phi: float = 0.74,
) -> None:
    r = diameter / 2.0
    spacing = diameter * spacing_factor
    margin = diameter * 1.5
    r_sph = r * (n / packing_phi) ** (1.0 / 3.0) * 1.05

    cx = cy = 0.0
    cz = height / 2.0

    if r_sph > r_cyl - margin or r_sph > height / 2.0 - margin:
        raise SystemExit(
            "Requested spheroid does not fit inside the cylinder with the "
            f"{margin*1e6:.1f} µm boundary margin"
        )

    a = spacing
    v0 = (a, 0.0, 0.0)
    v1 = (a * 0.5, a * math.sqrt(3) / 2.0, 0.0)
    v2 = (a * 0.5, a * math.sqrt(3) / 6.0, a * math.sqrt(6) / 3.0)

    n_lat = int(math.ceil(2.0 * r_sph / a)) + 3
    candidates: list[tuple[float, float, float, float]] = []
    for i in range(-n_lat, n_lat + 1):
        for j in range(-n_lat, n_lat + 1):
            for k in range(-n_lat, n_lat + 1):
                x = cx + i * v0[0] + j * v1[0] + k * v2[0]
                y = cy + i * v0[1] + j * v1[1] + k * v2[1]
                z = cz + i * v0[2] + j * v1[2] + k * v2[2]
                dist2 = (x - cx) ** 2 + (y - cy) ** 2 + (z - cz) ** 2
                if dist2 > r_sph ** 2:
                    continue
                if x * x + y * y > (r_cyl - margin) ** 2:
                    continue
                if z < margin or z > height - margin:
                    continue
                candidates.append((dist2, x, y, z))

    candidates.sort(key=lambda t: t[0])
    if len(candidates) < n:
        raise SystemExit(
            f"Only {len(candidates)} lattice sites inside spheroid+cylinder; need {n}"
        )

    pts = [(x, y, z) for _, x, y, z in candidates[:n]]
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["x", "y", "z"])
        for x, y, z in pts:
            w.writerow([f"{x:.8f}", f"{y:.8f}", f"{z:.8f}"])

    zs = [p[2] for p in pts]
    print(
        f"wrote {len(pts)} agents -> {out}\n"
        f"  d={diameter*1e6:.1f} µm  R_sph≈{r_sph*1e6:.1f} µm  "
        f"center_z={cz*1e3:.3f} mm\n"
        f"  z=[{min(zs)*1e3:.3f}, {max(zs)*1e3:.3f}] mm"
    )


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("-n", type=int, default=5000)
    p.add_argument("--diameter", type=float, default=20e-6)
    p.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path(__file__).resolve().parents[1]
        / "agents_positions-tissue_41_diffusion.csv",
    )
    args = p.parse_args()
    generate(args.n, args.diameter, args.output)


if __name__ == "__main__":
    main()
