# cagepipe

AMBER GAFF2 parametrization pipeline for metal-organic cages (Pd/Pt M2L4, M6L12, M12L24).
Produces tleap-ready force-field files (`ori_dry.prmtop`, `ori_dry.inpcrd`) from a cage
geometry (`inputfile.xyz`) and RESP charges (`inputfile.chg`).

## What it does

Per-cage chain (also expressed as a Snakemake DAG in `workflow/Snakefile`):

```
scr_1/inputfile.molden     # optional QM molden (Gaussian/Psi4/ORCA)
        |
        | cagepipe-respfit  (Multiwfn 2-stage RESP)
        v
inputfile.xyz + inputfile.chg
        |
        | cagepipe-pdb4munro  (OpenBabel for XYZ->PDB; antechamber for Sybyl types)
        v
bone.pdb + LA*.mol2, P*.mol2, Ltemp*_template.{pdb,mol2}
        |
        | cagepipe-munro --auto-from-pdb  (uses bundled gaff2.dat)
        v
munro.frcmod (+ ligand .frcmod files)
        |
        | cagepipe-tleapgen
        v
tleap.in
        |
        | tleap -s -f tleap.in
        v
ori_dry.{pdb,prmtop,inpcrd}
```

## Install

```bash
# 1. Create the conda env (brings in AmberTools, OpenBabel, Snakemake)
conda env create -f envs/cagepipe.yaml
conda activate cagepipe

# 2. The env file already does `pip install -e ..` for cagepipe itself, but
#    if you skipped that, do it manually:
pip install -e .
```

## Multiwfn (closed-source — install separately)

The RESP step uses [Multiwfn](https://sobereva.com/multiwfn/). Download the noGUI
Linux build and point cagepipe at it via either:

```bash
export CAGEPIPE_MULTIWFN=/path/to/Multiwfn_noGUI
# or
cagepipe-respfit ... --multiwfn /path/to/Multiwfn_noGUI
# or, in Snakemake:
snakemake --config multiwfn=/path/to/Multiwfn_noGUI ...
```

## Usage

### Option 1: direct CLI chain

```bash
cd <cage_dir>     # must contain inputfile.xyz and inputfile.chg
# Either run the four stages individually:
cagepipe-pdb4munro inputfile.xyz --chg inputfile.chg
cagepipe-munro     -p bone.pdb --auto-from-pdb
cagepipe-tleapgen  -p bone.pdb -o tleap.in
tleap -s -f tleap.in
# Or use the bundled wrapper:
bash /path/to/cagepipe/scripts/run_pipeline.sh
```

### Option 2: Snakemake (recommended for batch work)

```bash
# Single cage
snakemake --snakefile /path/to/cagepipe/workflow/Snakefile --cores 1 -d <cage_dir>

# Many cages, in parallel (limit to N at a time)
for d in 1013 1036 1052 ...; do
  snakemake --snakefile /path/to/cagepipe/workflow/Snakefile --cores 1 -d "$d" &
done; wait

# Dry-run the DAG without executing anything
snakemake --snakefile /path/to/cagepipe/workflow/Snakefile --cores 1 -d <cage_dir> -n
```

## Console scripts

| Name | Purpose |
| --- | --- |
| `cagepipe-pdb4munro`      | XYZ/PDB → `bone.pdb` + per-residue `*.mol2` + ligand templates |
| `cagepipe-munro`          | `bone.pdb` + GAFF → `munro.frcmod` (auto-runs antechamber+parmchk2) |
| `cagepipe-tleapgen`       | `bone.pdb` → `tleap.in` |
| `cagepipe-respfit`        | molden → `inputfile.chg` (whole-cage or differential RESP) |
| `cagepipe-chgass`         | per-residue MOL2s from pre-charged templates (M12L24 cleanup) |
| `cagepipe-seasoning`      | place N anions around the cage (random shells outside cavity) |
| `cagepipe-seasoning-inside` | place 1 anion (e.g. PFOA) inside the cavity |
| `cagepipe-gaff-typing`    | run antechamber + parmchk2 on template PDB(s) |

## Package layout

```
cagepipe/
├── pyproject.toml           # pip-installable package definition
├── envs/cagepipe.yaml       # conda environment (AmberTools, OpenBabel, snakemake, cagepipe)
├── workflow/                # Snakemake workflow
│   ├── Snakefile
│   └── rules/{resp,dry,solvate}.smk
├── config/config.yaml       # workflow defaults; override via --config k=v
├── src/cagepipe/            # Python package
│   ├── __init__.py
│   ├── *.py                 # pdb4munro, munro, tleapgen, respfit, chgass, ...
│   └── data/                # gaff2.dat, gaff.dat, leaprc.gaff2 (bundled)
└── scripts/                 # legacy bash wrappers (kept for back-compat)
    ├── run_pipeline.sh
    ├── respfit.sh
    └── respfit_cage.sh
```

## Notes / quirks

- The pipeline expects a Pd-N or Pt-N metal-coordinating cage. Bond detection uses a
  1.90 Å covalent cutoff, which naturally excludes metal–ligand bonds.
- `munro --auto-from-pdb` will run `antechamber` and `parmchk2` on every
  `*_template.pdb` it finds — make sure your `cagepipe` env is active so they're on PATH.
- The bundled `gaff2.dat` ships at `src/cagepipe/data/gaff2.dat`. Override with
  `cagepipe-munro -g /custom/gaff2.dat`.
- Multiwfn redistribution is restricted by its license — that's why it's not in the
  conda env file.
