#!/usr/bin/env bash
#SBATCH --job-name=respfit
#SBATCH --output=respfit.out
#SBATCH --time=96:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --export=ALL

# Per-cage respfit wrapper. Submit with sbatch --chdir=<cage_dir>:
#   sbatch --chdir=/path/to/cage  /path/to/cagepipe/respfit_cage.sh
# Reads scr_1/inputfile.molden and writes inputfile.chg in the cage dir.

module load multiwfn/noGUI
# SLURM's --chdir sets cwd to the cage dir; do not cd elsewhere.
# Requires the `cagepipe` conda env (provides respfit on PATH).
respfit scr_1/inputfile.molden -o inputfile.chg --keep-aux
