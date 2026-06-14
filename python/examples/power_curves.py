"""Power-vs-N curve for the two-group Gaussian design (d = 0.5). Writes a CSV for the
figure and prints the curve, including the Type M exaggeration ratio (which shrinks toward
1 as power grows -- a design-analysis point worth teaching)."""
import os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
SPEC = os.path.join(HERE, "..", "..", "spec", "examples")
BUILD = os.path.join(HERE, "..", "..", "build")
os.makedirs(BUILD, exist_ok=True)
from pilotr import load_spec, power_curve


def main():
    spec = load_spec(os.path.join(SPEC, "between_2group_gaussian.json"))
    grid = [20, 40, 60, 80, 100, 120, 140, 160]   # total N (per group = N / 2)
    curve = power_curve(spec, grid, n_sims=2000)

    with open(os.path.join(BUILD, "power_curve_gaussian.csv"), "w", newline="") as f:
        f.write("n_subject,power,type_m\n")
        for r in curve:
            f.write(f"{r['n_subject']},{r['power']},{r['type_m']}\n")

    print("=== power curve: two-group Gaussian, d = 0.5, alpha = .05, 2000 sims/point ===")
    print("  total N   per group   power   Type M")
    for r in curve:
        tm = "  n/a" if r["type_m"] != r["type_m"] else f"{r['type_m']:.2f}"
        print(f"  {r['n_subject']:5d}   {r['n_subject'] // 2:7d}    {r['power']:.3f}   {tm}")
    hit = next((r["n_subject"] for r in curve if r["power"] >= 0.80), None)
    print(f"  ~80% power reached by total N = {hit} ({hit // 2}/group)" if hit
          else "  (extend the grid to reach 80% power)")


if __name__ == "__main__":
    main()
