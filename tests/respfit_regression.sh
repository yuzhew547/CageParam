#!/usr/bin/env bash
#SBATCH --job-name=respfit_regr_1013
#SBATCH --output=respfit_test.log
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --export=ALL
#
# One-off regression test for the new cagepipe-respfit console script.
# Reads scr_1/inputfile.molden and writes inputfile.chg.respfit_test_<ts>
# in the cage dir. Does NOT touch the existing inputfile.chg.
#
# Usage:
#   sbatch --chdir=/home/gridsan/ywang6/PFASuptake/Cages/1013 \
#          /home/gridsan/ywang6/PFASuptake/Cages/cagepipe/tests/respfit_regression.sh
set -eo pipefail

TS=$(date +%Y%m%d_%H%M%S)
OUT_CHG="inputfile.chg.respfit_test_${TS}"

# Activate the AmberTools25 env (currently the only env where pip-installed
# cagepipe lives on PATH via ~/.local/bin); also brings python on PATH.
source /state/partition1/llgrid/pkg/anaconda/python-LLM-2023b/etc/profile.d/conda.sh
conda activate AmberTools25

# The cagepipe-respfit script is in ~/.local/bin (from `pip install -e .` we ran).
export PATH="$HOME/.local/bin:$PATH"

module load multiwfn/noGUI || true   # tolerate missing module on dev nodes

echo "[regression] Job started $(date)"
echo "[regression] PWD=$PWD"
echo "[regression] OUT_CHG=$OUT_CHG"
echo "[regression] which cagepipe-respfit: $(command -v cagepipe-respfit || echo MISSING)"
echo "[regression] which Multiwfn_noGUI: $(command -v Multiwfn_noGUI || echo MISSING)"
echo "[regression] SLURM_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK"

cagepipe-respfit scr_1/inputfile.molden -o "$OUT_CHG" --keep-aux

echo "[regression] Job finished $(date)"
echo "[regression] Output: $PWD/$OUT_CHG"

# Compare against the existing chg if present.
if [[ -f inputfile.chg ]]; then
    python - <<PY
import sys
def load(p):
    rows = []
    with open(p) as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 5:
                try: rows.append((parts[0], float(parts[1]), float(parts[2]), float(parts[3]), float(parts[4])))
                except ValueError: pass
    return rows
old = load("inputfile.chg")
new = load("${OUT_CHG}")
if len(old) != len(new):
    print(f"[regression] DIFFERS: atom count old={len(old)} new={len(new)}")
    sys.exit(0)
import math
max_dq = max(abs(o[4] - n[4]) for o,n in zip(old, new))
rms_dq = math.sqrt(sum((o[4] - n[4])**2 for o,n in zip(old, new))/len(old))
print(f"[regression] {len(old)} atoms, max|dq|={max_dq:.4f} e, RMS dq={rms_dq:.4f} e")
print(f"[regression] {'PASS' if max_dq < 1e-3 else ('CLOSE' if max_dq < 0.05 else 'DIFF')}")
PY
else
    echo "[regression] No baseline inputfile.chg to diff against."
fi
