# cagepipe

This package provides an automated workflow for parameterizing metal–organic cages (MOCs) in AMBER, including Pd/Pt-based assemblies such as M₂L₄, M₆L₁₂, and M₁₂L₂₄ cages (M=Pd/Pt). Starting from an optimized cage geometry (cage.xyz) and DFT-generated electronic structure data, the toolkit transfers existing force-field parameters template to the target cage, derives RESP atomic charges from Molden files, and assigns the charges to MOL2 files.

The package prepares AMBER-compatible simulation inputs by automatically generating tleap.in files for system construction and setup. It also supports addition of counterions, placement of guest molecules inside the cage cavity, and generation of cntl files for [alchemical free-energy calculations](https://github.com/Gallicchio-Lab/AToM-OpenMM).


## Installation

### Install via mamba (recommended):

```bash
mamba env create -f envs/cagepipe.yaml       
conda activate cagepipe
```
### Install via conda (slower):

```bash
conda env create -f envs/cagepipe.yaml       
conda activate cagepipe
```
### Verify the install:

```bash
which pdb4munro munro tleapgen respfit obabel tleap snakemake
# all paths should be inside $CONDA_PREFIX/bin
```

### Multiwfn (optional— install separately)

If you wish to use this package to fit the resp charge, you need to install [Multiwfn](https://sobereva.com/multiwfn/) separately. Download the noGUI Linux build and point cagepipe at it via either:

```bash
# export the path and add it to .bashrc
export CAGEPIPE_MULTIWFN=/path/to/Multiwfn_noGUI

# or when you need to specify the path of multiwfn everytime when you run respfit command
respfit ... --multiwfn /path/to/Multiwfn_noGUI

# or, in Snakemake:
snakemake --config multiwfn=/path/to/Multiwfn_noGUI ...
```

## Basic Usage

Use -h/--help to check how to use these commands. Here the most important functions are shown here.

### respfit: RESP charge fitting

### pdb4munro: PDB file preparation and mol2 file generation for tleap file

```bash
pdb4munro input [-o OUTPUT] [--chg CHG]
# Required arguments: input: input cage coordinates (pdb format or xyz file)
# Optional arguments: [-o OUTPUT] Output pdb file name of the cage system; [--chg CHG] resp charge file
# Example: pdb4munro input.xyz -o bone.pdb --chg input.chg
```

### seasoning: place one guest molecule inside the cage or multiple counterions around the cage

### munro: generate frcmod file for MOCs

```bash
munro -p PDB 
```

### tleapgen: generate tleap file to system

```bash
tleapgen -p PDB [-o OUTPUT] [--solvent SOLVENT]
# Required arguments: -p PDB: input pdb file
```

### cntlgen: generate cntl file for alchemical transfer method

```bash
cntlgen [--pdb PDB] [--resname RESNAME] [--out OUT]
# Optional arguments: [--pdb PDB] input system coordinates; [--resname RESNAME] resp charge file; [-o OUTPUT] Output cntl file name.
# Example: cntlgen cage_solv.pdb -out bone.pdb -resname GS1
```

```bash
# Either run the four stages individually:
pdb4munro inputfile.xyz --chg inputfile.chg
munro     -p bone.pdb --auto-from-pdb
tleapgen  -p bone.pdb -o tleap.in
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

### Option 3: ATM cntl generation

After solvation (`cage_solv.pdb` exists with cage + guest + counterions + water),
build the ATM control file:

```bash
cd <cage_dir>
cntlgen                          # PDB=gus_solv.pdb, guest resname=GS1 (filling default)
cntlgen --resname BFA            # legacy systems where BFA is the guest
cntlgen --template my.cntl       # use a custom template (only atom-list fields are rewritten)
```

`cntlgen` writes 0-based atom indices for `LIGAND_ATOMS`, `LIGAND_CM_ATOMS`,
`RCPT_CM_ATOMS`, and `POS_RESTRAINED_ATOMS`; all other ATM settings come from the
template. Receptor atoms are anything that is not the guest, not water (HOH,
WAT, T3P, T4P, T5P, OPC, OPC3, SOL, TIP*), not a monatomic ion (Na+, Cl-, K+,
...), and not a `seasoning` counterion (BFA..BFZ, BGA..). `POS_RESTRAINED_ATOMS`
picks up receptor residues whose names start with `P` or `M` (the metal centers
— override with `--metal-prefix`).

## Console scripts

| Name | Purpose |
| --- | --- |
| `pdb4munro`      | XYZ/PDB → `bone.pdb` + per-residue `*.mol2` + ligand templates |
| `munro`          | `bone.pdb` + GAFF → `munro.frcmod` (auto-runs antechamber+parmchk2) |
| `tleapgen`       | `bone.pdb` → `tleap.in` |
| `respfit`        | molden → `inputfile.chg` (whole-cage or differential RESP) |
| `chgass`         | per-residue MOL2s from pre-charged templates (M12L24 cleanup) |
| `seasoning`      | place N anions around the cage (random shells outside cavity) |
| `filling`        | parametrize a guest (.pdb or .xyz; `--autoparam` runs antechamber + parmchk2), place 1 copy (default resname `GS1`) inside the cavity, and optionally also place `--counterions N` BFA/BFB/... outside (subsumes `seasoning`) |
| `gaff-typing`    | run antechamber + parmchk2 on template PDB(s) |
| `cntlgen`        | generate `gus_solv.cntl` for ATM from a solvated PDB; auto-fills `LIGAND_ATOMS`, `RCPT_CM_ATOMS`, and `POS_RESTRAINED_ATOMS` (0-based) by residue-name classification (guest defaults to `GS1`) |

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
  `munro -g /custom/gaff2.dat`.
- Multiwfn redistribution is restricted by its license — that's why it's not in the
  conda env file.
