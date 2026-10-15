"""`power_report`: a sample-size justification in Markdown, from a power result.

The report is rendered from a record of counts and strings, and the R twin renders the same record
into the same bytes. fixtures/power_report_record.json is one such record, the case "mixed_full"
of tools/parity/report_cases.json, and fixtures/power_report_expected.md is its text. The R
package's tests/testthat/fixtures holds the same two files, and tools/parity/report_cross.py
renders every case of report_cases.json in both twins and requires identical text.
"""
import copy, hashlib, json, os, platform, sys, warnings
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import pytest

import pilotr
from pilotr import power_report, solve_curve, target_n
from pilotr.mcse import _rate_with_error, _type_s_interval
from pilotr.report import _power_report_record, _render

FIXTURES = os.path.join(os.path.dirname(__file__), "fixtures")


def _fixture(name, mode="r"):
    path = os.path.join(FIXTURES, name)
    if mode == "rb":
        with open(path, "rb") as f:
            return f.read()
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def _expected():
    # A checkout on Windows can give the file CRLF line endings, which are read back as the
    # newlines the report writes.
    return _fixture("power_report_expected.md", "rb").decode("utf-8").replace("\r\n", "\n")


def _fake_mixed(backend="statsmodels MixedLM (crossed variance components, REML; Wald z tests)"):
    """A power_mixed() result in the shape with one entry per focal effect under `effects`.

    Built by hand, which needs neither statsmodels nor any fitting. The second effect returned an
    estimate in 198 of the 200 fits.
    """
    return {
        "backend": backend, "n_sims": 200, "alpha": 0.05,
        "n_attempted": 200, "n_returned": 200, "n_converged": 151, "n_singular": 46,
        "n_warning": 3,
        "effects": {
            "effect": {"power": 141 / 200, "power_mcse": 0.0322, "power_lo": 0.638,
                       "power_hi": 0.763, "n_significant": 141, "true_effect": 0.06,
                       "mean_estimate": 0.061, "type_s": 0.0, "type_s_lo": 0.0,
                       "type_s_hi": 0.0265, "type_m": 1.2734, "type_m_mcse": 0.0213},
            "age": {"power": 37 / 198, "power_mcse": 0.0277, "power_lo": 0.139,
                    "power_hi": 0.247, "n_significant": 37, "true_effect": -0.02,
                    "mean_estimate": -0.021, "type_s": 2 / 37, "type_s_lo": 0.015,
                    "type_s_hi": 0.177, "type_m": 2.1, "type_m_mcse": 0.19},
        },
    }


def _fake_null():
    out = _fake_mixed()
    out["n_returned"] = 199
    for name, count in (("effect", 11), ("age", 9)):
        e = out["effects"][name]
        e.update(power=count / 199, n_significant=count, true_effect=0.0,
                 type_s=float("nan"), type_m=float("nan"))
    return out


def _fake_solved():
    return {"value": 51.37236, "lo": 40.12, "hi": 67.4, "level": 0.95, "target": 0.8, "se": 0.4,
            "dispersion": 1.3071, "x": "n_subject", "y": "power", "transform": "sqrt",
            "intercept": -2.0, "slope": 0.4, "n_points": 7, "x_min": 8.0, "x_max": 56.0,
            "n": 52, "n_lo": 41, "n_hi": 68}


def _rate_line(count, n):
    """The report's line for a rate, built here from the counts with the package's own helpers."""
    r = _rate_with_error(count, n, "p")
    return ("%.3f (Monte Carlo standard error %.3f; 95%% Wilson interval %.3f to %.3f)"
            % (r["p"], r["p_mcse"], r["p_lo"], r["p_hi"]))


def _two_group_spec(n=64, effect=5):
    return {
        "name": "d", "seed": 1,
        "units": {"subject": {"n": n}},
        "factors": [{"name": "group", "levels": ["a", "b"],
                     "contrasts": {"effect": [-0.5, 0.5]}, "between": "subject"}],
        "fixed": {"intercept": 100, "coefficients": {"effect": effect}},
        "response": {"family": "gaussian", "name": "score", "sigma": 10},
    }


def test_stored_record_renders_as_the_expected_text():
    assert _render(_fixture("power_report_record.json")) == _expected()


def test_the_fixtures_are_the_parity_case_and_the_r_twins_text():
    # The fixture is the case "mixed_full" of the parity cases, and the expected text is the R
    # twin's file, byte for byte, wherever the repository is present.
    root = os.path.join(os.path.dirname(__file__), "..", "..")
    cases = os.path.join(root, "tools", "parity", "report_cases.json")
    r_text = os.path.join(root, "r", "pilotr", "tests", "testthat", "fixtures",
                          "power-report-expected.md")
    if not (os.path.exists(cases) and os.path.exists(r_text)):
        pytest.skip("the repository's parity cases and R fixtures are not here")
    with open(cases, encoding="utf-8") as f:
        assert json.load(f)["mixed_full"] == _fixture("power_report_record.json")
    with open(r_text, "rb") as f:
        assert f.read().replace(b"\r\n", b"\n") == _expected().encode("utf-8")


def test_a_power_mixed_result_becomes_the_record_the_parity_case_holds():
    stored = _fixture("power_report_record.json")
    rec = _power_report_record(_fake_mixed(), solved=_fake_solved(), spec_file=None,
                               justification={"effect": stored["effects"][0]["justification"]},
                               null_result=_fake_null())
    # The R twin fits this design with lmerTest and a ladder, which the Python result cannot
    # carry, so the analysis and the software line come from the stored record.
    rec["analysis"] = stored["analysis"]
    rec["reproducibility"] = stored["reproducibility"]
    assert _render(rec) == _expected()


def test_each_rate_is_computed_in_the_report_from_its_counts():
    rep = power_report(_fake_mixed())
    assert isinstance(rep, str)
    assert "- Power: %s, from 141 of 200 replicates" % _rate_line(141, 200) in rep
    assert "- Power: %s, from 37 of 198 replicates" % _rate_line(37, 198) in rep
    lo, hi = _type_s_interval(2, 37, -0.02)
    assert "- Type S: 0.054 (95%% Wilson interval %.3f to %.3f)" % (lo, hi) in rep
    assert "- Type M: 2.100 (Monte Carlo standard error 0.190)" in rep


def test_the_primary_analysis_names_the_backend_alpha_and_sidedness():
    rep = power_report(_fake_mixed())
    assert ("- Model and test: statsmodels MixedLM (crossed variance components, REML; Wald z "
            "tests)\n") in rep
    assert "- Significance level: alpha = 0.05, two-sided\n" in rep
    one = power_report(_fake_mixed(
        "statsmodels MixedLM (crossed variance components, REML; one-sided Wald z tests, less)"))
    assert "- Significance level: alpha = 0.05, one-sided (less)\n" in one


def test_the_scalar_keys_of_an_older_power_mixed_result_are_read():
    # The shape power_mixed() returned before `effects`: one effect, the within contrast.
    x = {"backend": "statsmodels MixedLM (crossed variance components, REML; Wald z tests)",
         "n_sims": 50, "alpha": 0.05, "n_attempted": 50, "n_returned": 48, "n_converged": 40,
         "n_singular": 5, "n_warning": 3, "power": 30 / 48, "n_significant": 30,
         "true_effect": 0.05, "type_s": 0.0, "type_m": 1.4, "type_m_mcse": 0.05}
    rep = power_report(x, justification="A pilot estimate.")
    assert "### The within-factor contrast" in rep
    assert "- Power: %s, from 30 of 48 replicates" % _rate_line(30, 48) in rep


def test_a_power_result_with_effects_is_read_without_its_deprecated_keys():
    # A result that also carries the deprecated scalar keys, in a dict that warns when one is
    # read, is reported from `effects` alone.
    class Warns(dict):
        def __getitem__(self, key):
            if key in ("power", "n_significant", "true_effect", "type_s", "type_m"):
                warnings.warn("deprecated key %s" % key, FutureWarning, stacklevel=2)
            return dict.__getitem__(self, key)

        def get(self, key, default=None):
            if key in ("power", "n_significant", "true_effect", "type_s", "type_m"):
                warnings.warn("deprecated key %s" % key, FutureWarning, stacklevel=2)
            return dict.get(self, key, default)

    x = Warns(_fake_mixed())
    dict.update(x, x["effects"]["effect"])
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        rep = power_report(x)
    assert "### `age`" in rep


def test_a_live_power_mixed_result_and_its_null_are_reported_from_effects():
    pytest.importorskip("statsmodels")
    pytest.importorskip("pandas")
    from pilotr import power_mixed
    spec = {
        "name": "small", "seed": 5,
        "units": {"subject": {"n": 10}, "item": {"n": 6}},
        "factors": [{"name": "condition", "levels": ["a", "b"],
                     "contrasts": {"cond": [-0.5, 0.5]}, "vary_within": ["subject", "item"]}],
        "fixed": {"intercept": 0, "coefficients": {"cond": 0.5}},
        "random": {"subject": {"intercept_sd": 0.5, "slopes": {"cond": 0.2}},
                   "item": {"intercept_sd": 0.3, "slopes": {"cond": 0.1}}},
        "response": {"family": "gaussian", "name": "y", "sigma": 1},
    }
    res = power_mixed(spec, n_sims=2)
    null_spec = copy.deepcopy(spec)
    null_spec["fixed"]["coefficients"]["cond"] = 0.0
    null = power_mixed(null_spec, n_sims=2)
    with warnings.catch_warnings():
        warnings.simplefilter("error", FutureWarning)
        rep = power_report(res, null_result=null, justification={"cond": "Half a unit."})
        k, n = res["effects"]["cond"]["n_significant"], res["n_returned"]
    assert "### `cond`\n\n- True value: 0.5\n- Justification: Half a unit." in rep
    assert "from %d of %d replicates" % (k, n) in rep
    assert "## Power when the null hypothesis is true" in rep
    assert "statsmodels " in rep


def test_a_power_result_is_reported_over_all_its_replicates():
    pytest.importorskip("scipy")
    from pilotr import power
    res = power(_two_group_spec(), n_sims=40, alternative="greater")
    rep = power_report(res, justification="The smallest difference worth the cost of the course.")
    assert "### The difference between the two groups" in rep
    assert ("- Power: %.3f (Monte Carlo standard error %.3f; 95%% Wilson interval %.3f to %.3f), "
            "from %d of 40 replicates" % (res["power"], res["power_mcse"], res["power_lo"],
                                          res["power_hi"], res["n_significant"])) in rep
    assert "- Model and test: scipy.stats.ttest_ind, two-sample t with pooled variance" in rep
    assert "alpha = 0.05, one-sided (greater)" in rep
    assert "- Replicates: 40 simulated\n" in rep
    assert "- Justification: The smallest difference worth the cost of the course." in rep
    assert "- Fits:" not in rep
    assert "SciPy " in rep


def test_a_solved_curve_gives_the_sample_size_with_its_interval():
    rep = power_report(_fake_mixed(), solved=_fake_solved())
    assert "## Desired power and sample size" in rep
    assert "- Desired power: 0.8\n" in rep
    assert ("- Sample size: 52 (95% confidence interval 41 to 68), the exact solution 51.37 "
            "rounded up") in rep
    assert "the square root of `n_subject`" in rep
    assert "- Heterogeneity factor: 1.31" in rep

    curve = [{"n_subject": n, "power": p, "n_returned": m} for n, p, m in (
        (20, 0.21, 200), (40, 0.38, 200), (60, 0.52, 190), (80, 0.63, 200), (120, 0.78, 180),
        (160, 0.87, 200))]
    solved = target_n(curve, target=0.8)
    rep = power_report(_fake_mixed(), solved=solved)
    assert ("- Sample size: %d (95%% confidence interval %d to %d), the exact solution %.2f "
            "rounded up" % (solved["n"], solved["n_lo"], solved["n_hi"], solved["value"])) in rep

    effects = [{"cond": c, "power": p, "n_sims": 400} for c, p in (
        (0.1, 0.11), (0.2, 0.31), (0.3, 0.58), (0.4, 0.79), (0.5, 0.92))]
    mde = solve_curve(effects, target=0.8, transform="identity")
    rep = power_report(_fake_mixed(), solved=mde)
    assert "## Desired power and the solved value" in rep
    assert ("- Solved value of `cond`: %g (95%% confidence interval %g to %g)"
            % (mde["value"], mde["lo"], mde["hi"])) in rep
    with pytest.raises(ValueError) as caught:
        power_report(_fake_mixed(), solved={"value": 3})
    assert str(caught.value) == "`solved` must be the result of target_n() or solve_curve()."


def test_the_null_result_gives_the_rejection_rate():
    rep = power_report(_fake_mixed(), null_result=_fake_null())
    assert "## Power when the null hypothesis is true" in rep
    assert "- `effect`: %s, from 11 of 199 replicates" % _rate_line(11, 199) in rep
    assert "- Null condition: 200 simulated, 199 returned an estimate" in rep
    with pytest.raises(ValueError) as caught:
        power_report(_fake_mixed(), null_result=_fake_mixed())
    assert str(caught.value) == (
        "`null_result` must come from a design in which every focal effect is zero. Its effect "
        "'effect' has a true value of 0.06.")


def test_justifications_attach_to_their_effects():
    rep = power_report(_fake_mixed(), justification={"age": "Ageing studies report this."})
    assert ("### `age`\n\n- True value: -0.02\n- Justification: Ageing studies report this."
            in rep)
    assert "### `effect`\n\n- True value: 0.06\n- Justification: not given" in rep
    both = power_report(_fake_mixed(), justification="One reason for both.")
    assert both.count("One reason for both.") == 2
    with pytest.raises(ValueError) as caught:
        power_report(_fake_mixed(), justification={"cond": "x"})
    assert str(caught.value) == (
        "`justification` names 'cond', which is not a focal effect of `x`. Its focal effects are "
        "'effect' and 'age'.")
    form = ("`justification` must be one string, which applies to every focal effect, or strings "
            "named by focal effect.")
    for bad in (["a", "b"], 1, {"effect": 2}):
        with pytest.raises(ValueError) as caught:
            power_report(_fake_mixed(), justification=bad)
        assert str(caught.value) == form
    scalar = {"n_sims": 10, "alpha": 0.05, "alternative": "two.sided", "power": 0.5,
              "n_significant": 5, "true_effect": 5.0, "type_s": 0.0, "type_m": 1.5,
              "type_m_mcse": 0.1}
    with pytest.raises(ValueError) as caught:
        power_report(scalar, justification={"group": "x"})
    assert str(caught.value) == (
        "`x` has one effect, without a name, so `justification` must be a single string.")


def test_the_specification_file_gives_its_checksum_seed_and_units(tmp_path):
    path = tmp_path / "two group.json"
    path.write_text(json.dumps(_two_group_spec(), indent=2), encoding="utf-8")
    x = _fake_mixed()
    rep = power_report(x, spec_file=path)
    digest = hashlib.md5(path.read_bytes()).hexdigest()
    assert "- Specification: two group.json, MD5 checksum %s\n" % digest in rep
    assert "- Seed: 1\n" in rep
    assert "- Simulated units: subject 64\n" in rep
    assert ("- Software: pilotr %s, Python %s" % (pilotr.__version__, platform.python_version())
            in rep)
    none = power_report(x)
    assert "- Specification: not given, so neither its seed nor its checksum is recorded" in none
    missing = str(tmp_path / "no-such-spec.json")
    with pytest.raises(ValueError) as caught:
        power_report(x, spec_file=missing)
    assert str(caught.value) == (
        "`spec_file` must be the path of the specification file the analysis read. No file is "
        "found at '%s'." % missing)


def test_a_curve_or_anything_else_is_refused():
    msg = ("`x` must be the result of power() or power_mixed() at one design. A curve goes to "
           "target_n(), and its result to `solved`.")
    for bad in ([{"n_subject": 10, "power": 0.2, "n_sims": 50}], {"a": 1}, "power"):
        with pytest.raises(ValueError) as caught:
            power_report(bad)
        assert str(caught.value) == msg


def test_an_effect_without_a_true_value_or_estimates_is_reported_as_such():
    x = _fake_mixed()
    for e in x["effects"].values():
        e.update(true_effect=float("nan"), type_s=float("nan"), type_m=float("nan"),
                 type_m_mcse=float("nan"))
    x["effects"]["age"].update(power=float("nan"), n_significant=0)
    rep = power_report(x)
    assert "- True value: not given, so Type S and Type M are not defined" in rep
    assert ("- Power: not available, since no replicate returned an estimate for this effect"
            in rep)
    y = _fake_mixed()
    y["effects"]["age"]["type_m_mcse"] = float("nan")
    assert ("- Type M: 2.100 (too few significant replicates for a Monte Carlo standard error)"
            in power_report(y))


def test_the_report_is_markdown_and_the_same_on_every_call():
    rep = power_report(_fake_mixed(), solved=_fake_solved(), null_result=_fake_null())
    assert rep.startswith("# Sample-size justification\n")
    assert rep.endswith("\n") and not rep.endswith("\n\n")
    assert rep == power_report(copy.deepcopy(_fake_mixed()), solved=_fake_solved(),
                               null_result=_fake_null())
