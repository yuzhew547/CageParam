#!/usr/bin/env bash
# Run the cagepipe dry-parametrization pipeline in the current cage directory.
# Requires: inputfile.xyz + inputfile.chg in cwd; the `cagepipe` conda env active
# (which provides tleap + the cagepipe-* console scripts via `pip install -e .`).
#
# Two ways to use this script:
#   1. Direct CLI chain (this script, original behavior preserved):
#        cd <cage_dir> && /path/to/cagepipe/scripts/run_pipeline.sh
#   2. Snakemake (recommended once you have multiple cages):
#        snakemake --snakefile /path/to/cagepipe/workflow/Snakefile \
#                  --cores 1 -d <cage_dir>
set -eo pipefail

# Allow opt-in: source conda only if requested. Otherwise the caller is
# expected to have already activated the cagepipe env.
if [[ "${CAGEPIPE_AUTO_ACTIVATE:-0}" == "1" ]]; then
  source /state/partition1/llgrid/pkg/anaconda/python-LLM-2023b/etc/profile.d/conda.sh
  conda activate cagepipe
fi

for bin in cagepipe-pdb4munro cagepipe-munro cagepipe-tleapgen tleap; do
  command -v "$bin" >/dev/null 2>&1 || { echo "ERROR: '$bin' not on PATH. Did you 'conda activate cagepipe' and 'pip install -e <cagepipe>'?" >&2; exit 1; }
done

if [[ ! -f inputfile.xyz ]]; then
  echo "ERROR: inputfile.xyz not found in $PWD" >&2; exit 1
fi
if [[ ! -f inputfile.chg ]]; then
  echo "ERROR: inputfile.chg not found in $PWD" >&2; exit 1
fi

echo "===== [1/4] cagepipe-pdb4munro ====="
cagepipe-pdb4munro inputfile.xyz --chg inputfile.chg

echo "===== [2/4] cagepipe-munro --auto-from-pdb ====="
cagepipe-munro -p bone.pdb --auto-from-pdb       # GAFF auto-resolves to bundled gaff2.dat

echo "===== [3/4] cagepipe-tleapgen ====="
cagepipe-tleapgen -p bone.pdb -o tleap.in

echo "===== [4/4] tleap -s -f tleap.in ====="
tleap -s -f tleap.in

if [[ -f ori_dry.prmtop && -f ori_dry.inpcrd && -f ori_dry.pdb ]]; then
  echo "===== DONE: ori_dry.{pdb,prmtop,inpcrd} produced ====="
else
  echo "===== FAIL: missing ori_dry.* outputs ====="; exit 2
fi
