#!/usr/bin/env python3
"""Plots for the translated Chapter 6 outputs."""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "out"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def plot_sphere_rcs() -> None:
    path = OUT_DIR / "sphere_rcs_convergence.csv"
    if not path.exists():
        return

    rows = read_csv(path)
    ka = [float(row["ka"]) for row in rows]
    rcs_1 = [float(row["rcs_n1_over_pi_a2"]) for row in rows]
    rcs_5 = [float(row["rcs_n5_over_pi_a2"]) for row in rows]
    rcs_10 = [float(row["rcs_n10_over_pi_a2"]) for row in rows]
    rcs_50 = [float(row["rcs_n50_over_pi_a2"]) for row in rows]

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.semilogy(ka, rcs_1, "k:", label="N=1")
    ax.semilogy(ka, rcs_5, "k--", label="N=5")
    ax.semilogy(ka, rcs_10, "k-.", label="N=10")
    ax.semilogy(ka, rcs_50, "k-", label="N=50")
    ax.set_xlabel("ka")
    ax.set_ylabel("sigma / (pi a^2)")
    ax.set_title("Convergencia da RCS de esfera PEC")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT_DIR / "sphere_rcs_convergence.png", dpi=150)
    plt.close(fig)


def plot_mom3d_cuts() -> None:
    path = OUT_DIR / "mom3d_cuts.csv"
    if not path.exists():
        return

    grouped: dict[str, list[tuple[float, float]]] = defaultdict(list)
    for row in read_csv(path):
        grouped[row["series"]].append((float(row["coordinate_m"]), float(row["value"])))

    fig, ax = plt.subplots(figsize=(7, 4))

    if grouped["AA_computed"]:
        xy = sorted(grouped["AA_computed"])
        ax.plot([p[0] for p in xy], [p[1] for p in xy], "xb", label="Cut AA: this code")
    if grouped["AA_reference"]:
        xy = sorted(grouped["AA_reference"])
        ax.plot([p[0] for p in xy], [p[1] for p in xy], "-ob", label="Cut AA [Lit]")
    if grouped["BB_computed"]:
        xy = sorted(grouped["BB_computed"])
        ax.plot([p[0] for p in xy], [p[1] for p in xy], "+r", label="Cut BB: this code")
    if grouped["BB_reference"]:
        xy = sorted(grouped["BB_reference"])
        ax.plot([p[0] for p in xy], [p[1] for p in xy], "--sr", label="Cut BB [Lit]")

    ax.set_xlabel("X / lambda, Y / lambda")
    ax.set_ylabel("|J_x / H_inc|")
    ax.set_title("RWG plate current cuts")
    ax.grid(True, alpha=0.3)
    if ax.has_data():
        ax.legend()
    fig.tight_layout()
    fig.savefig(OUT_DIR / "mom3d_cuts.png", dpi=150)
    plt.close(fig)


def main() -> None:
    plot_sphere_rcs()
    plot_mom3d_cuts()
    print("Figuras geradas em Cap_06/out/")


if __name__ == "__main__":
    main()
