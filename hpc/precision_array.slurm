#!/bin/bash
# =============================================================================
# Oxford ARC SLURM array job: simulation-based precision/ROPE design analysis
# for the pilotr toolkit. One array task per sample size; each task runs N_SIMS
# replicate fits in parallel across its cores and writes one CSV to the project
# results dir. Modelled on the project's existing run_R.slurm / bayesian_model.slurm.
#
# Submit (after bootstrapping the R library once via scripts/bootstrap_rlib.sh):
#   sbatch precision_array.slurm                       # full sweep, N_SIMS=200
#   sbatch --export=ALL,N_SIMS=4 --array=0 --partition=devel precision_array.slurm   # smoke test
# =============================================================================
#SBATCH --account=educ-intract
#SBATCH --job-name=pilotr_prec
#SBATCH --array=0-7                 # 8 sample sizes: 20 30 40 50 60 70 80 100
#SBATCH --partition=short           # short <=12h
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8           # replicates run in parallel via mclapply
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/data/educ-intract/educ1242/pilotr_toolkit/logs/%x_%A_%a.out
#SBATCH --error=/data/educ-intract/educ1242/pilotr_toolkit/logs/%x_%A_%a.err

set -uo pipefail
conda deactivate 2>/dev/null || true
module purge 2>/dev/null || true
module load R/4.5.2-gfbf-2025b

PD=/data/educ-intract/educ1242/pilotr_toolkit
export R_LIBS=$PD/Rlib:${R_LIBS:-}        # lme4/lmerTest live here; jsonlite from the module
export TOOLKIT=$HOME/pilotr_toolkit/toolkit
export OUTDIR=$PD/results
export N_SIMS=${N_SIMS:-200}
export ROPE=${ROPE:-0.05}
export TASK_ID=$SLURM_ARRAY_TASK_ID
export OMP_NUM_THREADS=1                  # mclapply provides the parallelism
mkdir -p "$OUTDIR" "$PD/logs" "$PD/tmp"
export TMPDIR=$PD/tmp

echo "Host $(hostname) | task $SLURM_ARRAY_TASK_ID | cpus ${SLURM_CPUS_PER_TASK} | N_SIMS=$N_SIMS | $(date)"
Rscript "$HOME/pilotr_toolkit/scripts/precision_array.R"
echo "Exit $? | $(date)"
