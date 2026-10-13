"""Rounding a response: the Python twin rounds by R's algorithm, so the twins agree at decimal ties.

The specification rounds by the algorithm of R's round(x, digits) (src/nmath/fround.c, R 4.0.0 and
later), and `pilotr.core.fround` ports it. CPython's built-in round() rounds the exact binary
value, so near a decimal tie it can take the other candidate. round(1.05, 1) is 1 in R and 1.1 in
CPython.

The fixture under tests/fixtures was written by R itself:

    Rscript tools/parity/round_probes.R

Each row holds a double, a number of digits and R's round() of that double, with the doubles in
C's hexadecimal notation so that they read back exactly. The file holds 5,000 doubles within three
units in the last place of a decimal tie, at 1 to 4 digits, and 1,000 random doubles. It also holds
57 edge cases covering each branch of fround.c and 310 doubles at 31 numbers of digits from 16 to
324.
"""
import copy, csv, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from pilotr import simulate
from pilotr.core import fround

FIXTURE = os.path.join(os.path.dirname(__file__), "fixtures", "r_round_probes.csv")

# NEAR_TIE is the design of tools/parity/cases/round_near_tie.json. Its sigma is far below half a
# unit in the last place of 1.05, so every unrounded response is the double nearest 1.05.
NEAR_TIE = {"spec_version": "0.3", "name": "round_near_tie", "seed": 1,
            "units": {"subject": {"n": 10}},
            "fixed": {"intercept": 1.05, "coefficients": {}},
            "response": {"family": "gaussian", "name": "y", "sigma": 1e-17, "round": 1}}


def _probes():
    with open(FIXTURE, newline="", encoding="ascii") as f:
        return [(r["kind"], float.fromhex(r["x"]), float(r["digits"]), float.fromhex(r["r"]))
                for r in csv.DictReader(f)]


def test_fround_matches_r_on_every_probe():
    # Compared through float.hex(), which tells -0.0 from 0.0 and reads NaN as equal to itself.
    probes = _probes()
    assert len(probes) == 6367
    wrong = [(kind, x.hex(), d, want.hex(), fround(x, d).hex())
             for kind, x, d, want in probes if fround(x, d).hex() != want.hex()]
    assert not wrong, "%d of %d probes differ from R, the first being %s" % (
        len(wrong), len(probes), wrong[:5])


def test_the_probes_separate_r_from_cpython():
    # Without this, the test above could pass with CPython's round() itself.
    ties = [(x, int(d), want) for kind, x, d, want in _probes() if kind == "tie"]
    differ = sum(round(x, d).hex() != want.hex() for x, d, want in ties)
    assert differ >= 100  # 166 of the 5,000 when the fixture was written


def test_a_decimal_tie_goes_to_the_candidate_nearer_in_double_precision():
    # In double precision, 1.1 - 1.05 and 1.05 - 1 are equal and 10 is even, so R gives 1. The
    # double nearest 1.05 lies just above it, so CPython gives 1.1.
    assert fround(1.05, 1) == 1.0 and round(1.05, 1) == 1.1
    assert fround(4.5499999999999998, 1) == 4.6 and round(4.5499999999999998, 1) == 4.5


def test_a_response_at_a_decimal_tie_rounds_as_in_r():
    # R writes 1 in every row, where Python wrote 1.1.
    assert simulate(NEAR_TIE).column("y") == [1.0] * 10


def test_a_round_that_is_not_whole_counts_as_the_nearest_whole_number_of_digits():
    # Validation refuses such a value. Without validation, R rounds to floor(digits + 0.5)
    # digits, where Python raised TypeError.
    spec = copy.deepcopy(NEAR_TIE)
    spec["fixed"]["intercept"] = 1.23456
    three = copy.deepcopy(spec)
    spec["response"]["round"] = 2.5
    three["response"]["round"] = 3
    assert simulate(spec, validate=False).column("y") == simulate(three).column("y")
    assert simulate(three).column("y") == [1.235] * 10
