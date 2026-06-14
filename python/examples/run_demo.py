"""Demo: simulate both example designs, recover parameters, run a power/design analysis."""
import os, sys, statistics

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))          # import pilotr without installing
SPEC = os.path.join(HERE, "..", "..", "spec", "examples")
BUILD = os.path.join(HERE, "..", "..", "build")
os.makedirs(BUILD, exist_ok=True)

from pilotr import simulate, load_spec, power


def main():
    # ---- 1. Between-subjects Gaussian design ----
    d = simulate(os.path.join(SPEC, "between_2group_gaussian.json"))
    d.to_csv(os.path.join(BUILD, "py_between.csv"))
    g0 = [r["score"] for r in d.rows if r["group"] == "control"]
    g1 = [r["score"] for r in d.rows if r["group"] == "treatment"]
    print("=== between_2group_gaussian ===")
    print(f"  N = {len(d)} ; control mean = {statistics.mean(g0):.3f}, "
          f"treatment mean = {statistics.mean(g1):.3f}")
    print(f"  observed difference = {statistics.mean(g1) - statistics.mean(g0):.3f} "
          f"(one noisy n=64 sample; true effect = 5.0. The Monte Carlo mean below recovers it.)")

    # ---- 2. Crossed mixed-effects reaction-time design ----
    d2 = simulate(os.path.join(SPEC, "crossed_mixed_rt.json"))
    d2.to_csv(os.path.join(BUILD, "py_crossed.csv"))
    rel = [r["RT"] for r in d2.rows if r["condition"] == "related"]
    unr = [r["RT"] for r in d2.rows if r["condition"] == "unrelated"]
    print("\n=== crossed_mixed_rt (subjects x items x condition) ===")
    print(f"  N = {len(d2)} rows ; mean RT related = {statistics.mean(rel):.1f} ms, "
          f"unrelated = {statistics.mean(unr):.1f} ms")
    print(f"  priming effect = {statistics.mean(unr) - statistics.mean(rel):.1f} ms "
          f"(positive => unrelated slower, as specified by cond beta = 0.05 on log scale)")
    print("  first 3 rows:")
    for r in d2.head(3):
        print("   ", r)

    # ---- continuous predictors + interactions + continuous random slopes ----
    d3 = simulate(os.path.join(SPEC, "reading_time_continuous.json"))
    d3.to_csv(os.path.join(BUILD, "py_continuous.csv"))
    print("\n=== reading_time_continuous (continuous predictors + interactions) ===")
    print(f"  N = {len(d3)} rows; columns = {d3.columns}")
    print("  predictors vary by item/subject; 3 interactions; by-subject random slopes on")
    print("  continuous predictors -> recover with precision_design_analysis.R")

    # ---- additional grouping factor: subjects nested in higher-level clusters ----
    d4 = simulate(os.path.join(SPEC, "nested_clusters.json"))
    d4.to_csv(os.path.join(BUILD, "py_nested.csv"))
    print("\n=== nested_clusters (subjects nested in clusters) ===")
    print(f"  N = {len(d4)} rows; columns = {d4.columns}")

    # ---- Beta family (proportions) + partial crossing (item subset per subject) ----
    for nm, fn in (("beta", "beta_proportion.json"), ("partial", "partial_crossing.json")):
        dd = simulate(os.path.join(SPEC, fn))
        dd.to_csv(os.path.join(BUILD, f"py_{nm}.csv"))
        print(f"=== {fn} === N = {len(dd)} rows; columns = {dd.columns}")

    # ---- 3. Simulation-based power + design analysis (Type S / Type M) ----
    spec = load_spec(os.path.join(SPEC, "between_2group_gaussian.json"))
    res = power(spec, n_sims=2000)
    print("\n=== simulation-based power + design analysis ===")
    print(f"  design: n=64 (32/group), d=0.5, alpha=.05, {res['n_sims']} simulations")
    print(f"  power            = {res['power']:.3f}")
    print(f"  Type S error     = {res['type_s']:.4f}  (sign errors among significant)")
    print(f"  Type M (exagg.)  = {res['type_m']:.3f}  (mean |estimate|/true | significant)")
    print(f"  true effect = {res['true_effect']:.2f}, mean estimate = {res['mean_estimate']:.3f}")


if __name__ == "__main__":
    main()
