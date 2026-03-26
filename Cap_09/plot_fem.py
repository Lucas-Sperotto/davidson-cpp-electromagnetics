from pathlib import Path
import csv
import matplotlib.pyplot as plt

BASE_DIR = Path(__file__).resolve().parent
OUT_DIR = BASE_DIR / "out"

def read_csv(filename):
    with open(filename, 'r', newline='') as f:
        reader = csv.reader(f)
        rows = list(reader)
    header = rows[0]
    data = list(zip(*rows[1:]))  # transpose
    cols = {h: [float(x) for x in col] for h, col in zip(header, data)}
    return cols

# Grafico de convergencia
conv = read_csv(OUT_DIR / "fem_convergence.csv")
plt.figure()
plt.loglog(conv["h_over_lambda"], conv["rms_err"], marker='o')
plt.xlabel("h/λ")
plt.ylabel("RMS error")
plt.title("Convergence (FEM vs Analítico)")
plt.grid(True, which='both')
plt.savefig(OUT_DIR / "convergence.png", dpi=200)
plt.close()

# Perfis FEM x Analitico para cada estagio
stage_files = sorted(OUT_DIR.glob("fem_profile_stage_*.csv"))
for fname in stage_files:
    cols = read_csv(fname)
    plt.figure()
    plt.plot(cols["z"], cols["V_fem"], label="FEM")
    plt.plot(cols["z"], cols["V_analitico"], linestyle='--', label="Analítico")
    plt.xlabel("z")
    plt.ylabel("V(z)")
    plt.title(fname.stem)
    plt.legend()
    plt.grid(True)
    out = fname.with_suffix(".png")
    plt.savefig(out, dpi=200)
    plt.close()

print(f"Graficos gerados em {OUT_DIR}: convergence.png e fem_profile_stage_#.png")
