from __future__ import annotations

import argparse
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


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Traducao do cyl_test.m para comparar diferentes rodadas do problema do cilindro PEC."
    )
    parser.add_argument("--old", required=True, help="CSV do caso antigo.")
    parser.add_argument("--new", required=True, help="CSV do caso novo.")
    parser.add_argument("--new-longer", required=True, help="CSV do caso novo com dominio/tempo maior.")
    parser.add_argument("--output", default=str(OUT_DIR / "cyl_test_comparison.png"))
    args = parser.parse_args()

    old_path = Path(args.old)
    new_path = Path(args.new)
    new_longer_path = Path(args.new_longer)
    output_path = Path(args.output)

    old_time, old_signal = load_signal(old_path)
    new_time, new_signal = load_signal(new_path)
    new_longer_time, new_longer_signal = load_signal(new_longer_path)

    plt.figure(figsize=(8, 4.5))
    plt.plot(old_time, old_signal, label=old_path.stem)
    plt.plot(new_time, new_signal, ":", label=new_path.stem)
    plt.plot(new_longer_time, new_longer_signal, "--", label=new_longer_path.stem)
    plt.xlabel("t [ns]")
    plt.ylabel(r"$E_y$ [V/m]")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()

    print(f"Grafico salvo em: {output_path}")


if __name__ == "__main__":
    main()
