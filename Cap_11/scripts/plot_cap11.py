#!/usr/bin/env python3
"""Plots for the translated Chapter 11 outputs."""

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


def to_float_rows(rows: list[dict[str, str]]) -> list[dict[str, float]]:
    converted: list[dict[str, float]] = []
    for row in rows:
        converted.append({key: float(value) for key, value in row.items()})
    return converted


def column_as_float(rows: list[dict[str, str]], name: str) -> list[float]:
    return [float(row[name]) for row in rows]


def plot_fetd_fdtd() -> None:
    path = OUT_DIR / "fetd_fdtd_ctlN_funcs.csv"
    if not path.exists():
        return

    rows = to_float_rows(read_csv(path))
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot([row["y"] for row in rows], [row["we1"] for row in rows], "--", label="w_e1")
    ax.plot([row["y"] for row in rows], [row["we2"] for row in rows], "-.", label="w_e2")
    ax.plot([row["y"] for row in rows], [row["we1_dot_we2"] for row in rows], "-", label="w_e1 . w_e2")
    quad_rows = [row for row in rows if int(row["is_quad_point"]) == 1]
    ax.plot([row["y"] for row in quad_rows], [row["quad_value"] for row in quad_rows], "ko", label="Quadrature points")
    ax.set_xlabel("y")
    ax.set_ylabel("f(y)")
    ax.set_title("CTLN functions for FETD-FDTD note")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fetd_fdtd_ctlN_funcs.png", dpi=150)
    plt.close(fig)


def plot_test_tet_quad() -> None:
    path = OUT_DIR / "test_tet_quad_summary.csv"
    if not path.exists():
        return

    rows = read_csv(path)
    labels = [row["test"] for row in rows]
    errors = column_as_float(rows, "relative_error")

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.bar(labels, errors, color="white", edgecolor="black")
    ax.set_ylabel("Relative error")
    ax.set_title("Tetrahedral quadrature rule check")
    ax.set_yscale("log")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "test_tet_quad_errors.png", dpi=150)
    plt.close(fig)


def plot_modes(csv_name: str, png_name: str, title: str) -> None:
    path = OUT_DIR / csv_name
    if not path.exists():
        return

    rows = to_float_rows(read_csv(path))
    ranks = [row["rank"] for row in rows]
    kc = [row["kc"] for row in rows]
    rel_err = [row["relative_error"] for row in rows]

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(ranks, kc, marker="o")
    ax.set_xlabel("Mode rank")
    ax.set_ylabel("k_c")
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT_DIR / png_name, dpi=150)
    plt.close(fig)

    physical_rows = [(row["rank"], row["relative_error"]) for row in rows if row["relative_error"] >= 0.0]
    if physical_rows:
        fig, ax = plt.subplots(figsize=(7, 4))
        ax.semilogy([row[0] for row in physical_rows], [row[1] for row in physical_rows], marker="o")
        ax.set_xlabel("Mode rank")
        ax.set_ylabel("Relative error")
        ax.set_title(f"{title} - relative error")
        ax.grid(True, which="both", alpha=0.3)
        fig.tight_layout()
        fig.savefig(OUT_DIR / png_name.replace(".png", "_relerr.png"), dpi=150)
        plt.close(fig)


def main() -> None:
    plot_fetd_fdtd()
    plot_test_tet_quad()
    plot_modes("eigdata3D_v0_modes.csv", "eigen3d_v0_modes.png", "Eigen3D v0")
    plot_modes("eigdata3D_ctlN_modes.csv", "eigen3d_ctln_modes.png", "Eigen3D CTLN")
    plot_modes("eigdata3D_ltqn_modes.csv", "eigen3d_ltqn_modes.png", "Eigen3D LTQN")
    print("Figuras geradas em Cap_11/out/")


if __name__ == "__main__":
    main()
