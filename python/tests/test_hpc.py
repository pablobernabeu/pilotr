"""The SLURM array job in the repository's hpc/ folder runs on any cluster.

The job ran only under the author's account and directory layout: it charged a fixed account,
wrote its logs and results under the author's project paths, ran a copy of its runner kept outside
the repository and needed a bootstrap script the repository never held. The job is R, and
r/pilotr/tests/testthat/test-hpc.R also runs it against precision_curve(). That test skips under
R CMD check, which builds the package without the folders around it, so this one holds the same
guard where CI runs the suite on the whole checkout. Twinned with test-hpc.R.
"""

import os
import re

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.normpath(os.path.join(HERE, "..", ".."))
BANNED = ("educ-intract", "educ1242", "/data/educ", "pilotr_toolkit", "bootstrap_rlib")


def _hpc_text(name):
    path = os.path.join(REPO, "hpc", name)
    if not os.path.exists(path):
        pytest.skip(f"hpc/{name} not available outside the repository")
    with open(path, encoding="utf-8") as fh:
        return fh.read()


@pytest.mark.parametrize("name", ["precision_array.slurm", "precision_array.R"])
def test_the_hpc_job_names_nothing_of_the_authors_cluster(name):
    text = _hpc_text(name)
    for banned in BANNED:
        assert banned not in text, f"{name} holds {banned}"


def test_the_slurm_wrapper_carries_placeholders_and_runs_the_runner_beside_it():
    text = _hpc_text("precision_array.slurm")
    lines = text.splitlines()
    assert "EDIT THESE" in text
    assert "#SBATCH --account=EDIT_ME_ACCOUNT" in lines
    assert any(re.match(r"PROJECT_DIR=.*EDIT THIS", line) for line in lines)
    # A log path is resolved where the job runs, so an absolute one names somebody's directory.
    assert not any(re.match(r"#SBATCH +--(output|error)=/", line) for line in lines)
    assert not any(line.startswith("module load R/") for line in lines)
    assert any(line.startswith("SCRIPT_DIR=${SLURM_SUBMIT_DIR") for line in lines)
    assert 'Rscript "$SCRIPT_DIR/precision_array.R"' in text


def test_the_runner_calls_precision_design_from_the_installed_package():
    text = _hpc_text("precision_array.R")
    assert re.search(r"^[^#\n]*library\(pilotr\)", text, re.MULTILINE)
    assert not re.search(r"^[^#\n]*\bsource\(", text, re.MULTILINE)
    assert re.search(r"precision_design\(spec, focal = focal, rope = rope, n_sims = n_sims,", text)
