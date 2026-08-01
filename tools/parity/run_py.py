"""Cross-language parity harness, Python side.

Simulates every specification in spec/examples/ from the *local* package sources (not an
installed copy) and writes one canonical dump per case, in the same format as run_r.R.
compare.py then diffs the two sets byte for byte.

Each specification is run twice: "asis", and "noround" with ``response.round`` deleted. The
second case matters because six of the eight shipped examples set ``round``, which quantises
away exactly the last-ulp divergences this harness exists to detect.

Usage: python tools/parity/run_py.py [output-dir]
"""

from __future__ import annotations

import json
import os
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "python"))

from pilotr.simulate import simulate  # noqa: E402  (needs the sys.path line above)


def _check_formatter() -> None:
    """printf must round-trip a double exactly, or the whole comparison is meaningless."""
    probes = [1 / 3, 0.1, 3.141592653589793, 1e-300, sys.float_info.max, -2.5e-17,
              0.0, 1.0, -1.0, sys.float_info.epsilon, 1 - sys.float_info.epsilon / 2]
    for v in probes:
        if float("%.17g" % v) != v:
            raise SystemExit("'%%.17g' does not round-trip %r" % v)


def _cell(v) -> str:
    """Language-neutral cell rendering: strings verbatim, everything else at 17 significant
    digits as a double. Integers therefore print as "1", matching R."""
    if isinstance(v, str):
        return v
    return "%.17g" % float(v)


def _dump(ds, path: str) -> None:
    with open(path, "wb") as f:                     # binary: LF endings on every platform
        f.write((",".join(ds.columns) + "\n").encode())
        for r in ds.rows:
            f.write((",".join(_cell(r[c]) for c in ds.columns) + "\n").encode())


def main() -> None:
    _check_formatter()
    out_dir = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "tools", "parity", "out", "py")
    os.makedirs(out_dir, exist_ok=True)

    # Two sources. spec/examples/ is what ships, and must never change silently.
    # tools/parity/cases/ holds adversarial specifications that deliberately exercise what the
    # shipped examples do not: long random-effect vectors (so the Cholesky and matrix-vector
    # inner products run to three terms or more, where the two languages' built-in reductions
    # diverge), long linear predictors, interaction random slopes, and a Gaussian response, so
    # that no libm transcendental masks or manufactures a difference.
    spec_dirs = [os.path.join(ROOT, "spec", "examples"),
                 os.path.join(ROOT, "tools", "parity", "cases")]
    specs = [(d, n) for d in spec_dirs if os.path.isdir(d)
             for n in sorted(os.listdir(d)) if n.endswith(".json")]
    if not specs:
        raise SystemExit("no specifications found in " + " or ".join(spec_dirs))

    for spec_dir, name in specs:
        base = name[:-len(".json")]
        for variant in ("asis", "noround"):
            with open(os.path.join(spec_dir, name)) as f:
                spec = json.load(f)
            if variant == "noround":
                spec["response"].pop("round", None)
            _dump(simulate(spec), os.path.join(out_dir, "%s.%s.txt" % (base, variant)))
        print("simulated:", base)

    print("Python dumps written to", out_dir)


if __name__ == "__main__":
    main()
