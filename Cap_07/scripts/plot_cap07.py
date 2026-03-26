#!/usr/bin/env python3
"""Plots for the translated Chapter 7 outputs."""

from __future__ import annotations

import csv
from collections import defaultdict
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


def group_by(rows: list[dict[str, float]], key: str) -> dict[float, list[dict[str, float]]]:
    grouped: dict[float, list[dict[str, float]]] = defaultdict(list)
    for row in rows:
        grouped[row[key]].append(row)
    return dict(grouped)


def plot_scalar_pot() -> None:
    dtm = read_csv(OUT_DIR / "scalar_pot_dtm.csv")
    dte = read_csv(OUT_DIR / "scalar_pot_dte.csv")
    integrand = read_csv(OUT_DIR / "scalar_pot_integrand.csv")
    region1 = read_csv(OUT_DIR / "scalar_pot_region1.csv")
    region2_lambda = read_csv(OUT_DIR / "scalar_pot_region2_lambda.csv")
    region2_t = read_csv(OUT_DIR / "scalar_pot_region2_t.csv")
    region3 = read_csv(OUT_DIR / "scalar_pot_region3.csv")

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot([row["lambda_over_k0"] for row in dtm], [row["re_d_tm"] for row in dtm], label="Real")
    ax.plot([row["lambda_over_k0"] for row in dtm], [row["im_d_tm"] for row in dtm], ":", label="Imag")
    ax.set_xlabel("lambda / k0")
    ax.set_ylabel("D_TM")
    ax.set_title("TM denominator")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "scalar_pot_dtm.png", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot([row["lambda_over_k0"] for row in dte], [row["abs_d_te"] for row in dte], color="black")
    ax.set_xlabel("lambda / k0")
    ax.set_ylabel("|D_TE|")
    ax.set_title("TE denominator magnitude")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "scalar_pot_dte.png", dpi=150)
    plt.close(fig)

    x = [row["lambda_over_k0"] for row in integrand]
    y_real = [row["re_integrand"] for row in integrand]
    y_imag = [row["im_integrand"] for row in integrand]

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(x, y_real, label="Real")
    ax.plot(x, y_imag, "-.", label="Imag")
    ax.set_xlabel("lambda / k0")
    ax.set_ylabel("F")
    ax.set_title("Scalar-potential integrand")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT_DIR / "scalar_pot_integrand.png", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(x, y_real, label="Real")
    ax.plot(x, y_imag, "-.", label="Imag")
    ax.set_xlim(0.9, 1.4)
    ax.set_ylim(-0.6, 0.2)
    ax.set_xlabel("lambda / k0")
    ax.set_ylabel("F")
    ax.set_title("Scalar-potential integrand (zoom)")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT_DIR / "scalar_pot_integrand_zoom.png", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot([row["t"] for row in region1], [row["re_integrand"] for row in region1], label="Real")
    ax.plot([row["t"] for row in region1], [row["im_integrand"] for row in region1], "-.", label="Imag")
    ax.set_xlabel("t")
    ax.set_ylabel("Region 1 integrand")
    ax.set_title("Region 1")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT_DIR / "scalar_pot_region1.png", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot([row["lambda_over_k0"] for row in region2_lambda],
            [row["re_full"] for row in region2_lambda], label="F")
    ax.plot([row["lambda_over_k0"] for row in region2_lambda],
            [row["re_subtracted"] for row in region2_lambda], "-.", label="F - F_sing")
    ax.set_xlabel("lambda / k0")
    ax.set_ylabel("Real part")
    ax.set_title("Region 2 subtraction")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT_DIR / "scalar_pot_region2_lambda.png", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot([row["t"] for row in region2_t], [row["abs_f"] for row in region2_t], "-.", label="|F|")
    ax.plot([row["t"] for row in region2_t], [row["abs_f_sing"] for row in region2_t], ":", label="|F_sing|")
    ax.plot([row["t"] for row in region2_t],
            [row["re_integrand_normalized"] for row in region2_t], label="Normalized difference")
    ax.set_xlabel("t")
    ax.set_ylabel("Magnitude / normalized value")
    ax.set_title("Region 2")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT_DIR / "scalar_pot_region2_t.png", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot([row["lambda_over_k0"] for row in region3],
            [row["re_full"] for row in region3], label="F")
    ax.plot([row["lambda_over_k0"] for row in region3],
            [row["re_minus_static"] for row in region3], ":", label="F - F_static")
    ax.set_xlabel("lambda / k0")
    ax.set_ylabel("Real part")
    ax.set_title("Region 3")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT_DIR / "scalar_pot_region3.png", dpi=150)
    plt.close(fig)


def plot_v_pot_eps() -> None:
    rows = read_csv(OUT_DIR / "v_pot_eps.csv")
    grouped = group_by(rows, "eps_r_prime")

    fig, ax = plt.subplots(figsize=(7, 4))
    styles = ["k-", "k:", "k-.", "k--"]
    for style, (eps_r_prime, values) in zip(styles, sorted(grouped.items())):
        ax.loglog([row["k0_rho"] for row in values],
                  [row["normalized_abs_v"] for row in values],
                  style, label=f"eps_r={eps_r_prime:g}")
    ax.set_xlabel("k0 rho")
    ax.set_ylabel("|V_n|")
    ax.set_title("Scalar potential vs relative permittivity")
    ax.legend()
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "v_pot_eps.png", dpi=150)
    plt.close(fig)


def plot_v_pot_height() -> None:
    mag_rows = read_csv(OUT_DIR / "v_pot_height_mag.csv")
    phase_rows = read_csv(OUT_DIR / "v_pot_height_phase.csv")
    grouped_mag = group_by(mag_rows, "b")
    grouped_phase = group_by(phase_rows, "b")

    fig, ax = plt.subplots(figsize=(7, 4))
    styles = ["k-", "k:", "k-."]
    for style, (b_value, values) in zip(styles, sorted(grouped_mag.items(), reverse=True)):
        ax.loglog([row["k_rho"] for row in values],
                  [row["normalized_abs_v"] for row in values],
                  style, label=f"b={b_value:g}")
    ax.set_xlabel("k rho")
    ax.set_ylabel("|V_n|")
    ax.set_title("Scalar potential vs height")
    ax.legend()
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "v_pot_height_mag.png", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4))
    for style, (b_value, values) in zip(styles, sorted(grouped_phase.items(), reverse=True)):
        ax.semilogx([row["k_rho"] for row in values],
                    [row["phase_deg"] for row in values],
                    style, label=f"b={b_value:g}")
    ax.set_xlabel("k rho")
    ax.set_ylabel("Phase angle V_n [deg]")
    ax.set_title("Scalar potential phase vs height")
    ax.legend()
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "v_pot_height_phase.png", dpi=150)
    plt.close(fig)


def plot_mom_som() -> None:
    impedance_rows = read_csv(OUT_DIR / "mom_som_impedance.csv")
    current_rows = read_csv(OUT_DIR / "mom_som_currents.csv")
    grouped_currents = group_by(current_rows, "freq_hz")

    fig, ax = plt.subplots(figsize=(7, 4))
    freq_ghz = [row["freq_ghz"] for row in impedance_rows]
    ax.plot(freq_ghz, [row["re_zin"] for row in impedance_rows], label="Real")
    ax.plot(freq_ghz, [row["im_zin"] for row in impedance_rows], "--", label="Imag")
    ax.set_xlabel("Freq [GHz]")
    ax.set_ylabel("Impedance [ohm]")
    ax.set_title("Printed dipole input impedance")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT_DIR / "mom_som_impedance.png", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4))
    for freq_hz, values in sorted(grouped_currents.items()):
        ax.plot([row["x_over_L"] for row in values],
                [row["abs_i"] for row in values],
                label=f"{freq_hz/1e9:.2f} GHz")
    ax.set_xlabel("x / L")
    ax.set_ylabel("|I|")
    ax.set_title("Printed dipole current profiles")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT_DIR / "mom_som_currents.png", dpi=150)
    plt.close(fig)


def main() -> None:
    plot_scalar_pot()
    plot_v_pot_eps()
    plot_v_pot_height()
    plot_mom_som()
    print("Figuras geradas em Cap_07/out/")


if __name__ == "__main__":
    main()
