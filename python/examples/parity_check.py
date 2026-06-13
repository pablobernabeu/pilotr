"""Verify that the R and Python implementations produce identical data from one spec.

Reads the CSVs written by both run_demo scripts (toolkit/build/) and compares them
row-by-row. Run order:  python run_demo.py  ->  Rscript run_demo.R  ->  python parity_check.py
"""
import csv, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
BUILD = os.path.join(HERE, "..", "..", "build")


def load(path):
    with open(path) as f:
        rd = csv.reader(f)
        header = next(rd)
        rows = [r for r in rd]
    return header, rows


def compare(name, response_col):
    py = os.path.join(BUILD, f"py_{name}.csv")
    r = os.path.join(BUILD, f"r_{name}.csv")
    for p in (py, r):
        if not os.path.exists(p):
            print(f"  [{name}] MISSING {os.path.basename(p)} -- run both demos first.")
            return False
    h1, r1 = load(py)
    h2, r2 = load(r)
    ok = True
    if h1 != h2:
        print(f"  [{name}] header mismatch: {h1} vs {h2}"); ok = False
    if len(r1) != len(r2):
        print(f"  [{name}] row count mismatch: {len(r1)} vs {len(r2)}"); ok = False
    ci = h1.index(response_col)
    label_cols = [i for i, c in enumerate(h1) if c not in (response_col,)]
    max_diff = 0.0
    label_mismatch = 0
    for a, b in zip(r1, r2):
        max_diff = max(max_diff, abs(float(a[ci]) - float(b[ci])))
        for i in label_cols:
            if a[i] != b[i]:
                label_mismatch += 1
    verdict = "IDENTICAL" if (ok and max_diff < 1e-6 and label_mismatch == 0) else "DIFFERS"
    print(f"  [{name}] rows={len(r1)}  max abs diff ({response_col}) = {max_diff:.3e}  "
          f"label mismatches = {label_mismatch}  -> {verdict}")
    return ok and max_diff < 1e-6 and label_mismatch == 0


def main():
    print("=== R vs Python parity check (same spec + seed) ===")
    a = compare("between", "score")
    b = compare("crossed", "RT")
    print("\nRESULT:", "ALL IDENTICAL [PASS]" if (a and b) else "PARITY FAILURE [FAIL]")
    sys.exit(0 if (a and b) else 1)


if __name__ == "__main__":
    main()
