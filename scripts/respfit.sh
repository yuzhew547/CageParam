#!/usr/bin/env bash
#SBATCH --job-name=respfit
#SBATCH --output=respfit.out
#SBATCH --time=96:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --export=ALL

# Generic respfit SLURM wrapper. Submit from a cage directory containing
# input.molden:
#   cd /path/to/some/cage/
#   sbatch /path/to/cagepipe/respfit.sh
# Multiwfn is OpenMP, so use --cpus-per-task (not ntasks). respfit.py reads
# $SLURM_CPUS_PER_TASK and writes a local settings.ini with that nthreads.

module load multiwfn/noGUI
cd "${SLURM_SUBMIT_DIR:-$PWD}"
# Requires the `cagepipe` conda env (provides cagepipe-respfit on PATH).
cagepipe-respfit input.molden -o input.chg --keep-aux
