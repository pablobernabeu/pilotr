"""Verify that the R and Python implementations produce identical data from one spec.

Reads the CSVs written by both run_demo scripts (build/) and compares them
cell-by-cell: numeric columns within tolerance, string columns exactly. Run order:
  python run_demo.py  ->  Rscript run_demo.R  ->  python parity_check.py

Passing --csv PATH also writes the per-design summary as a CSV, which is how the R
package's cross-language vignette gets its table: the vignette cannot run the check at
build time (that needs both languages installed), so it ships a fixture regenerated from
this script rather than a hand-transcribed one.
"""
import csv, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
BUILD = os.path.join(HERE, "..", "..", "build")


def load(path):
    with open(path) as f:
        rd = csv.reader(f)
        return next(rd), [r for r in rd]


def _num(x):
    try:
        return float(x)
    except ValueError:
        return None


def compare(name, summary=None):
    py, r = os.path.join(BUILD, f"py_{name}.csv"), os.path.join(BUILD, f"r_{name}.csv")
    for p in (py, r):
        if not os.path.exists(p):
            print(f"  [{name}] MISSING {os.path.basename(p)}. Run both demos first.")
            return False
    h1, r1 = load(py)
    h2, r2 = load(r)
    ok = h1 == h2 and len(r1) == len(r2)
    if not ok:
        print(f"  [{name}] structure mismatch: headers {h1==h2}, rows {len(r1)} vs {len(r2)}")
        return False
    max_diff, str_mismatch = 0.0, 0
    for a, b in zip(r1, r2):
        for av, bv in zip(a, b):
            fa, fb = _num(av), _num(bv)
            if fa is not None and fb is not None:
                max_diff = max(max_diff, abs(fa - fb))
            elif av != bv:
                str_mismatch += 1
    verdict = "IDENTICAL" if (max_diff < 1e-6 and str_mismatch == 0) else "DIFFERS"
    print(f"  [{name}] rows={len(r1)} cols={len(h1)}  max abs numeric diff = {max_diff:.3e}  "
          f"string mismatches = {str_mismatch}  -> {verdict}")
    if summary is not None:
        summary.append({"design": name, "rows": len(r1), "columns": len(h1),
                        "max_abs_diff": "0" if max_diff == 0 else f"{max_diff:.3e}",
                        "categorical_mismatches": str_mismatch})
    return max_diff < 1e-6 and str_mismatch == 0


def main():
    out = None
    if "--csv" in sys.argv:
        out = sys.argv[sys.argv.index("--csv") + 1]
    summary = [] if out else None
    print("=== R vs Python parity check (same spec + seed) ===")
    results = [compare(n, summary)
               for n in ("between", "crossed", "continuous", "nested", "beta", "partial")]
    if out:
        with open(out, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=list(summary[0]))
            w.writeheader()
            w.writerows(summary)
        print(f"\nSummary written to {out}")
    print("\nRESULT:", "ALL IDENTICAL [PASS]" if all(results) else "PARITY FAILURE [FAIL]")
    sys.exit(0 if all(results) else 1)


if __name__ == "__main__":
    main()
