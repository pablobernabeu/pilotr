"""Exercise the realistic response families: ordinal (Likert) and Poisson (counts)."""
import os, sys, statistics
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
SPEC = os.path.join(HERE, "..", "..", "spec", "examples")
from simdgp import simulate


def main():
    # ---- Ordinal / Likert via cumulative-logit ----
    d = simulate(os.path.join(SPEC, "ordinal_likert_between.json"))
    print("=== ordinal_likert_between (5-point Likert, cumulative-logit) ===")
    for grp in ("control", "treatment"):
        ratings = [r["rating"] for r in d.rows if r["group"] == grp]
        dist = Counter(ratings)
        props = {k: round(dist[k] / len(ratings), 3) for k in sorted(dist)}
        print(f"  {grp:9s} mean = {statistics.mean(ratings):.2f}  category proportions = {props}")
    print("  (treatment shifted to higher categories, as set by grp beta = 1.0 on the logit scale)")

    # ---- Poisson counts via log link ----
    d2 = simulate(os.path.join(SPEC, "poisson_counts_between.json"))
    print("\n=== poisson_counts_between (counts, log link) ===")
    for grp, expected in (("control", 4.09), ("treatment", 6.11)):
        counts = [r["count"] for r in d2.rows if r["group"] == grp]
        print(f"  {grp:9s} mean count = {statistics.mean(counts):.2f}  (expected ~{expected})")
    print("  (intercept ln(5)=1.609; grp beta 0.4 on log scale => rate ratio exp(0.4)=1.49)")


if __name__ == "__main__":
    main()
