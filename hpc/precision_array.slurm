#!/bin/bash
# =============================================================================
# SLURM array job: simulation-based precision and ROPE design analysis with
# pilotr. One array task per sample size. Each task calls precision_design()
# for its size, with the replicates spread over its cores, and writes one CSV to
# $PROJECT_DIR/results. Bound together with rbind(), the CSVs form the curve that
# precision_curve() returns, ready for solve_curve().
#
# The R library at $PROJECT_DIR/Rlib has to hold pilotr, lme4 and lmerTest,
# installed once on a login node (see "Running at scale" in the README). Install
# pilotr from this clone, so that the package and this job come from one commit.
#
# Submit from the hpc/ directory of the clone, so that this wrapper finds
# precision_array.R beside itself. SLURM writes each task's log there too.
#   sbatch precision_array.slurm                              # full sweep, N_SIMS=200
#   sbatch --export=ALL,N_SIMS=4 --array=0 precision_array.slurm   # smoke test
#   ROPE=-0.02,0.05 sbatch precision_array.slurm              # a region given as an interval
# =============================================================================
# ---- EDIT THESE two values before submitting --------------------------------
# 1. The account (allocation) the scheduler charges. The scheduler reads #SBATCH
#    directives before the shell runs, so this one cannot come from a variable.
#SBATCH --account=EDIT_ME_ACCOUNT
# 2. PROJECT_DIR, just after the directives below: a directory you can write to,
#    holding the R library at Rlib/.
# Add a --partition directive if your cluster has no default partition, or if
# its default does not allow the time requested below.
# -----------------------------------------------------------------------------
#SBATCH --job-name=pilotr_prec
#SBATCH --array=0-7                 # 8 sample sizes: 20 30 40 50 60 70 80 100
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8           # one worker per core, forked from the task's R process
#SBATCH --mem=16G
#SBATCH --time=02:00:00

set -uo pipefail
conda deactivate 2>/dev/null || true
module purge 2>/dev/null || true
module load R    # or the versioned R module your cluster names

PROJECT_DIR=$HOME/pilotr_runs    # EDIT THIS: writable, holds Rlib/ (see above)

# The runner is the one saved next to this wrapper. sbatch runs a spooled copy of
# the wrapper, so $0 does not name the saved location. SLURM_SUBMIT_DIR is the
# directory sbatch was invoked from, which the submit instructions above pin to
# hpc/. $0 covers running the wrapper directly.
SCRIPT_DIR=${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}
if [ ! -f "$SCRIPT_DIR/precision_array.R" ]; then
  echo "precision_array.R is not in $SCRIPT_DIR. Submit from the hpc/ directory of the clone." >&2
  exit 1
fi

export R_LIBS=$PROJECT_DIR/Rlib:${R_LIBS:-}   # pilotr, lme4 and lmerTest
export OUTDIR=$PROJECT_DIR/results
export N_SIMS=${N_SIMS:-200}
# ROPE is a half-width or an interval written lower,upper. sbatch's --export separates variables
# with commas and so cannot carry an interval, which is set in the submitting shell, as above.
export ROPE=${ROPE:-0.05}
export TASK_ID=$SLURM_ARRAY_TASK_ID
export OMP_NUM_THREADS=1                      # one thread per worker, the workers parallelise
mkdir -p "$OUTDIR" "$PROJECT_DIR/tmp"
export TMPDIR=$PROJECT_DIR/tmp

echo "Host $(hostname) | task $SLURM_ARRAY_TASK_ID | cpus ${SLURM_CPUS_PER_TASK} | N_SIMS=$N_SIMS | $(date)"
Rscript "$SCRIPT_DIR/precision_array.R"
status=$?
echo "Exit $status | $(date)"
exit $status    # so that sacct marks a task whose R run failed as FAILED
