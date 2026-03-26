from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


SCRIPT_DIR = Path(__file__).resolve().parent
OUT_DIR = SCRIPT_DIR.parent / "out"


def load_signal(path: Path) -> tuple[np.ndarray, np.ndarray]:
    data = np.genfromtxt(path, delimiter=",", names=True, dtype=float, encoding="utf-8")
    names = set(data.dtype.names or ())
    required = {"tempo_ns", "Ey_Vpm"}
    missing = required.difference(names)
    if missing:
        raise ValueError(f"{path} nao possui colunas obrigatorias: {sorted(missing)}")
    return np.atleast_1d(data["tempo_ns"]), np.atleast_1d(data["Ey_Vpm"])


def load_meta(path: Path | None) -> dict[str, str]:
    if path is None or not path.exists():
        return {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        return {row["key"]: row["value"] for row in reader if "key" in row and "value" in row}


def make_label(default_label: str, meta: dict[str, str]) -> str:
    d_cell = meta.get("d_cell")
    poly_m = meta.get("poly_m")
    if d_cell is not None and poly_m is not None:
        return f"d={d_cell}; m={poly_m}"
    return default_label


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Traducao do PMLperformance.m para comparar desempenho de diferentes configuracoes de PML."
    )
    parser.add_argument("--thin", default=str(OUT_DIR / "ey_point1_pml.csv"))
    parser.add_argument("--thick", required=True, help="CSV da segunda configuracao PML.")
    parser.add_argument("--reference", required=True, help="CSV de referencia.")
    parser.add_argument("--thin-meta", default=str(OUT_DIR / "ey_point1_pml_meta.csv"))
    parser.add_argument("--thick-meta", default=None)
    parser.add_argument("--out-prefix", default="pml_performance")
    args = parser.parse_args()

    thin_path = Path(args.thin)
    thick_path = Path(args.thick)
    reference_path = Path(args.reference)
    thin_meta_path = Path(args.thin_meta) if args.thin_meta else None
    thick_meta_path = Path(args.thick_meta) if args.thick_meta else None

    thin_time, thin_signal = load_signal(thin_path)
    thick_time, thick_signal = load_signal(thick_path)
    reference_time, ref = load_signal(reference_path)

    thin_meta = load_meta(thin_meta_path)
    thick_meta = load_meta(thick_meta_path)

    length = min(len(thin_time), len(thick_time), len(reference_time))
    if length < 2:
        raise ValueError("Series muito curtas para comparacao.")

    time_ns = reference_time[:length]
    ref = ref[:length]
    thin_signal = thin_signal[:length]
    thick_signal = thick_signal[:length]

    refl_thin = thin_signal - ref
    refl_thick = thick_signal - ref

    gate = int(round(0.75 * length))
    scale = np.max(np.abs(ref))
    if scale == 0.0:
        scale = 1.0

    refl_thin_gated = np.zeros_like(refl_thin)
    refl_thick_gated = np.zeros_like(refl_thick)
    refl_thin_gated[:gate] = refl_thin[:gate]
    refl_thick_gated[:gate] = refl_thick[:gate]

    refl_thin_norm = refl_thin_gated / scale
    refl_thick_norm = refl_thick_gated / scale

    thin_label = make_label(thin_path.stem, thin_meta)
    thick_label = make_label(thick_path.stem, thick_meta)

    gamma_png = OUT_DIR / f"{args.out_prefix}_gamma.png"
    gamma_db_png = OUT_DIR / f"{args.out_prefix}_gamma_db.png"

    plt.figure(figsize=(8, 4.5))
    plt.plot(time_ns[:gate], refl_thin_norm[:gate], "k-", label=thin_label)
    plt.plot(time_ns[:gate], refl_thick_norm[:gate], "k-.", label=thick_label)
    plt.xlabel("t [ns]")
    plt.ylabel(r"$\Gamma_{norm}$")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(gamma_png, dpi=300)

    eps = 1e-30
    plt.figure(figsize=(8, 4.5))
    plt.plot(time_ns[:gate], 20.0 * np.log10(np.maximum(np.abs(refl_thin_norm[:gate]), eps)), "k-", label=thin_label)
    plt.plot(time_ns[:gate], 20.0 * np.log10(np.maximum(np.abs(refl_thick_norm[:gate]), eps)), "k-.", label=thick_label)
    plt.xlabel("t [ns]")
    plt.ylabel(r"$|\Gamma_{norm}|$ [dB]")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(gamma_db_png, dpi=300)

    print(f"Graficos salvos em: {gamma_png} e {gamma_db_png}")


if __name__ == "__main__":
    main()
