"""Cross-language parity harness, comparison and golden-file check.

Diffs the dumps written by run_r.R and run_py.py byte for byte, and checks the R dumps
against recorded SHA-256 hashes so that an unintended change to the generative core is caught
even when both languages change together.

For every differing cell it reports the distance in units in the last place, which separates
a genuine logic error (large ulp gap, or a differing string) from a floating-point
accumulation-order artefact (one or two ulps).

Usage:
    python tools/parity/compare.py              # compare, and check against golden.json
    python tools/parity/compare.py --update     # rewrite golden.json from the current R dumps

Exit status is 0 when the two languages agree and the hashes match, 1 otherwise.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import struct
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
OUT = os.path.join(ROOT, "tools", "parity", "out")
GOLDEN = os.path.join(ROOT, "tools", "parity", "golden.json")
TOLERANCE = os.path.join(ROOT, "tools", "parity", "tolerance.json")

MAX_REPORTED = 5      # differing cells shown per case, so a systematic fault stays readable


def _tolerances() -> tuple[int, dict[str, int], set[str]]:
    """Per-case ulp allowance, and the transcendental classification that decides which
    cases the golden anchor may pin. See tolerance.json for why each entry is what it is.
    """
    if not os.path.exists(TOLERANCE):
        return 0, {}, set()
    with open(TOLERANCE) as f:
        cfg = json.load(f)
    default_ulp = int(cfg.get("default_max_ulp", 0))
    case_ulp = dict(cfg.get("cases", {}))
    transcendental = set(cfg.get("transcendental", []))
    # A case that needs an ulp allowance needs it because a libm function shaped its bytes,
    # so leaving it out of the transcendental list would re-admit it to the golden anchor.
    unclassified = [c for c, ulp in case_ulp.items() if ulp and c not in transcendental]
    if unclassified:
        raise SystemExit("tolerance.json is inconsistent: %s carr%s an ulp allowance but "
                         "%s not classified transcendental"
                         % (", ".join(sorted(unclassified)),
                            "ies" if len(unclassified) == 1 else "y",
                            "is" if len(unclassified) == 1 else "are"))
    return default_ulp, case_ulp, transcendental


def _ulps(a: float, b: float) -> float:
    """Signed-magnitude ordinal distance between two doubles, in units in the last place."""
    if a == b:
        return 0.0
    if math.isnan(a) or math.isnan(b) or math.isinf(a) or math.isinf(b):
        return math.inf

    def ordinal(x: float) -> int:
        (n,) = struct.unpack("<q", struct.pack("<d", x))
        return n if n >= 0 else -(n & 0x7FFFFFFFFFFFFFFF) - (1 << 63)

    return float(abs(ordinal(a) - ordinal(b)))


def _sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()


def _describe(r_cell: str, p_cell: str) -> str:
    try:
        gap = _ulps(float(r_cell), float(p_cell))
    except ValueError:
        return "R=%s Py=%s (non-numeric)" % (r_cell, p_cell)
    return "R=%s Py=%s (%s ulp)" % (r_cell, p_cell, "inf" if math.isinf(gap) else "%d" % gap)


def _compare_case(r_path: str, p_path: str, max_ulp: int) -> tuple[int, int, int, list[str]]:
    """Return (failing cells, tolerated cells, total cells, sample of readable differences).

    A cell fails when it differs by more than `max_ulp`, or differs in a way no ulp count can
    describe (a non-numeric cell, an infinity, a NaN). Cells within the allowance are counted
    separately rather than ignored, so that tolerated drift stays visible.
    """
    with open(r_path, "rb") as f:
        r_lines = f.read().decode().splitlines()
    with open(p_path, "rb") as f:
        p_lines = f.read().decode().splitlines()

    if len(r_lines) != len(p_lines):
        return 1, 0, max(len(r_lines), len(p_lines)), [
            "row count differs: R has %d lines, Python %d" % (len(r_lines), len(p_lines))]
    if r_lines[0] != p_lines[0]:
        return 1, 0, len(r_lines), ["header differs: R=%s Py=%s" % (r_lines[0], p_lines[0])]

    header = r_lines[0].split(",")
    n_fail = n_tol = n_cells = 0
    report: list[str] = []
    for i, (rl, pl) in enumerate(zip(r_lines[1:], p_lines[1:]), start=1):
        if rl == pl:
            n_cells += len(header)
            continue
        rc, pc = rl.split(","), pl.split(",")
        n_cells += len(rc)
        for j, (a, b) in enumerate(zip(rc, pc)):
            if a == b:
                continue
            try:
                gap = _ulps(float(a), float(b))
            except ValueError:
                gap = math.inf
            if gap <= max_ulp:
                n_tol += 1
                continue
            n_fail += 1
            if len(report) < MAX_REPORTED:
                col = header[j] if j < len(header) else "col%d" % j
                report.append("row %d, %s: %s" % (i, col, _describe(a, b)))
    return n_fail, n_tol, n_cells, report


def main() -> int:
    update = "--update" in sys.argv
    r_dir, p_dir = os.path.join(OUT, "r"), os.path.join(OUT, "py")
    if not os.path.isdir(r_dir) or not os.path.isdir(p_dir):
        raise SystemExit("run tools/parity/run_r.R and tools/parity/run_py.py first")

    cases = sorted(f for f in os.listdir(r_dir) if f.endswith(".txt"))
    if not cases:
        raise SystemExit("no dumps found in " + r_dir)

    golden = {}
    if os.path.exists(GOLDEN) and not update:
        with open(GOLDEN) as f:
            golden = json.load(f)["hashes"]

    default_ulp, case_ulp, transcendental = _tolerances()
    failures, new_hashes = 0, {}
    print("%-42s %-10s %s" % ("case", "parity", "detail"))
    print("-" * 96)
    for case in cases:
        r_path, p_path = os.path.join(r_dir, case), os.path.join(p_dir, case)
        new_hashes[case] = _sha256(r_path)
        if not os.path.exists(p_path):
            print("%-42s %-10s %s" % (case, "MISSING", "no Python dump"))
            failures += 1
            continue
        allowance = case_ulp.get(case, default_ulp)
        n_fail, n_tol, n_cells, report = _compare_case(r_path, p_path, allowance)
        if n_fail:
            status, detail = "DIFFERS", "%d/%d cells (%.2f%%)" % (
                n_fail, n_cells, 100.0 * n_fail / max(n_cells, 1))
        elif n_tol:
            status, detail = "ok", "%d/%d cells within the %d-ulp allowance" % (
                n_tol, n_cells, allowance)
        else:
            status, detail = "ok", "bit-identical"
        print("%-42s %-10s %s" % (case, status, detail))
        for line in report:
            print("%-42s %-10s   %s" % ("", "", line))
        if n_fail:
            failures += 1

    # The anchor covers exactly the cases whose bytes are IEEE-754-exact
    # arithmetic: the Gaussian cases, which apply no transcendental to anything
    # that reaches the dump. A case classified transcendental in tolerance.json
    # is excluded because it cannot be anchored at all: its values pass through
    # exp(), log() or pow(), whose rounding IEEE-754 leaves to the maths
    # library, and Windows (MinGW/UCRT) and Linux (glibc) do not share a libm,
    # so the same golden.json could fail on whichever platform did not record
    # it. That some of those cases (beta, Poisson, ordinal, and the rounded
    # dumps of the exp families) have hashed identically on every platform
    # measured so far is an accident of the platforms measured, not a
    # guarantee, which is why the classification is by construction rather
    # than by observed agreement. The transcendental cases are gated by the
    # R-versus-Python comparison above instead, which is the contract, at an
    # allowance of zero ulp unless tolerance.json grants one.
    anchored = {c: h for c, h in new_hashes.items() if c not in transcendental}

    if update:
        with open(GOLDEN, "w") as f:
            json.dump({"note": "SHA-256 of the IEEE-exact (non-transcendental) R dumps in "
                               "tools/parity/out/r; regenerate with "
                               "`python tools/parity/compare.py --update`",
                       "hashes": anchored}, f, indent=2, sort_keys=True)
            f.write("\n")
        print("\ngolden.json updated (%d anchored cases; %d transcendental cases are "
              "gated by the cross-language comparison alone)"
              % (len(anchored), len(new_hashes) - len(anchored)))
        return 1 if failures else 0

    if golden:
        drifted = [c for c, h in anchored.items() if golden.get(c) != h]
        unanchorable = sorted(c for c in golden if c in transcendental)
        if unanchorable:
            # A stale golden.json from before the classification, or a case reclassified
            # without regenerating: either way the file pins bytes the rule says it must not.
            print("\ngolden.json anchors transcendental case(s): %s; regenerate it with "
                  "`python tools/parity/compare.py --update`" % ", ".join(unanchorable))
            failures += 1
        missing = [c for c in golden if c not in new_hashes and c not in transcendental]
        if drifted or missing:
            print("\ngolden-file drift:")
            for c in drifted:
                print("  %s: expected %s, got %s"
                      % (c, (golden.get(c) or "absent")[:16], new_hashes[c][:16]))
            for c in missing:
                print("  %s: case disappeared" % c)
            failures += len(drifted) + len(missing)
        else:
            print("\ngolden-file hashes match for all %d anchored cases"
                  " (%d cases are classified transcendental and are gated by the"
                  " cross-language comparison alone)"
                  % (len(anchored), len(new_hashes) - len(anchored)))
    else:
        print("\nno golden.json yet; create one with --update")

    print("\n%d of %d cases failed" % (failures, len(cases)))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
