"""Power-vs-N curve for the two-group Gaussian design (d = 0.5). Writes a CSV for the
figure and prints the curve, including the Type M exaggeration ratio (which shrinks towards
1 as power grows, a design-analysis point worth teaching) and the sample size the curve
implies for 80% power."""
import os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
SPEC = os.path.join(HERE, "..", "..", "spec", "examples")
BUILD = os.path.join(HERE, "..", "..", "build")
os.makedirs(BUILD, exist_ok=True)
from pilotr import load_spec, power_curve, target_n


def main():
    spec = load_spec(os.path.join(SPEC, "between_2group_gaussian.json"))
    grid = [20, 40, 60, 80, 100, 120, 140, 160]   # total N (per group = N / 2)
    curve = power_curve(spec, grid, n_sims=2000)

    with open(os.path.join(BUILD, "power_curve_gaussian.csv"), "w", newline="") as f:
        # The replicate counts go into the file too, so the saved curve can be solved later
        # without rerunning the simulation.
        f.write("n_subject,power,type_m,n_sims,n_significant\n")
        for r in curve:
            f.write(f"{r['n_subject']},{r['power']},{r['type_m']},"
                    f"{r['n_sims']},{r['n_significant']}\n")

    print("=== power curve: two-group Gaussian, d = 0.5, alpha = .05, 2000 sims/point ===")
    print("  total N   per group   power   Type M")
    for r in curve:
        tm = "  n/a" if r["type_m"] != r["type_m"] else f"{r['type_m']:.2f}"
        print(f"  {r['n_subject']:5d}   {r['n_subject'] // 2:7d}    {r['power']:.3f}   {tm}")
    # The first grid point at or above 0.80 is a property of where the grid falls, not of the
    # design. target_n() fits the curve and inverts the fit, and refuses rather than
    # extrapolate when the grid does not reach the target.
    try:
        solved = target_n(curve, target=0.80)
        print(f"  80% power at total N = {solved['n']} ({solved['n'] // 2}/group), "
              f"95% interval {solved['n_lo']} to {solved['n_hi']}")
    except ValueError as e:
        print(f"  {e}")


if __name__ == "__main__":
    main()
