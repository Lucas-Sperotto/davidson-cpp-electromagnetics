#!/usr/bin/env python3
"""Plots for the translated Chapter 4 outputs."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "out"


def read_csv(path: Path) -> list[dict[str, float]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        rows: list[dict[str, float]] = []
        for row in reader:
            rows.append({key: float(value) for key, value in row.items()})
    return rows


def plot_static_mom() -> None:
    path = OUT_DIR / "static_mom_charge.csv"
    if not path.exists():
        return

    rows = read_csv(path)
    z = [row["z_m"] for row in rows]
    charge = [row["line_charge_pC_per_m"] for row in rows]

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.bar(z, charge, width=(z[1] - z[0]) if len(z) > 1 else 0.2, color="white", edgecolor="black")
    ax.set_xlabel("Distance along wire [m]")
    ax.set_ylabel("Line charge density [pC/m]")
    ax.set_title("Static MoM wire charge distribution")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "static_mom_charge.png", dpi=150)
    plt.close(fig)


def plot_thin_dipole() -> None:
    path = OUT_DIR / "thin_dipole_current.csv"
    if not path.exists():
        return

    rows = read_csv(path)
    z = [row["z_over_lambda"] for row in rows]
    i_delta = [row["abs_i_delta_gap"] for row in rows]
    i_frill = [row["abs_i_mag_frill"] for row in rows]

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(z, i_delta, label="Delta gap")
    ax.plot(z, i_frill, "-.", label="Magnetic frill")
    ax.set_xlabel("z / lambda")
    ax.set_ylabel("|I| [A/m]")
    ax.set_title("Thin dipole current distribution")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "thin_dipole_current.png", dpi=150)
    plt.close(fig)


def plot_mom_current() -> None:
    path = OUT_DIR / "mom_tm_current.csv"
    if not path.exists():
        return

    rows = read_csv(path)
    phi = [row["phi_deg"] for row in rows]
    single = [row["abs_js_single_over_h_inc"] for row in rows]
    quad = [row["abs_js_quad_over_h_inc"] for row in rows]

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(phi, single, "+", label="Single point")
    ax.plot(phi, quad, "^", label="Quadrature")
    ax.set_xlabel("Phi [deg]")
    ax.set_ylabel("|J_s / H_inc|")
    ax.set_title("TM current on PEC cylinder")
    ax.set_xlim(0.0, 180.0)
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "mom_tm_current.png", dpi=150)
    plt.close(fig)


def plot_mom_rcs() -> None:
    path = OUT_DIR / "mom_tm_rcs.csv"
    if not path.exists():
        return

    rows = read_csv(path)
    ka = [row["ka"] for row in rows]
    a_over_lambda = [row["a_over_lambda"] for row in rows]
    quad_pi_a = [row["rcs_quad_over_pi_a"] for row in rows]
    single_pi_a = [row["rcs_single_over_pi_a"] for row in rows]
    exact_lambda = [row["exact_over_lambda"] for row in rows]
    quad_lambda = [row["rcs_quad_over_lambda"] for row in rows]
    single_lambda = [row["rcs_single_over_lambda"] for row in rows]

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.loglog(ka, quad_pi_a, "-", label="Quadrature")
    ax.loglog(ka, single_pi_a, "-.", label="Single point")
    ax.set_xlabel("ka")
    ax.set_ylabel("sigma_TM / (pi a)")
    ax.set_title("TM echo width normalized by pi a")
    ax.legend()
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "mom_tm_rcs_pi_a.png", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(a_over_lambda, exact_lambda, "k-", label="Eigenfunction")
    ax.plot(a_over_lambda, quad_lambda, "k:", label="Quadrature")
    ax.plot(a_over_lambda, single_lambda, "k-.", label="Single point")
    ax.set_xlabel("a / lambda")
    ax.set_ylabel("sigma_TM / lambda")
    ax.set_title("TM echo width normalized by lambda")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "mom_tm_rcs_lambda.png", dpi=150)
    plt.close(fig)


def main() -> None:
    plot_static_mom()
    plot_thin_dipole()
    plot_mom_current()
    plot_mom_rcs()
    print("Figuras geradas em Cap_04/out/")


if __name__ == "__main__":
    main()
