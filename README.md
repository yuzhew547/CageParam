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

respfit command can autodetect the integraty of the cage molecule. For small cage, it will average the charge of same ligand and keep the charge difference in different ligand (heteroleptic cages). For large system, it will take forvever to get molden file/fit the charge. User can use ML4 fragment to parametrize the system. The respfit will compute the charge difference of the ligand before and after bind one metal and apply the difference to both side of the ligand. The charge of the metal will be applied to all metal centers.

```bash
respfit [-o OUTPUT] molden
# Required arguments: molden molden file for the cage obtained from DFT calculation
# Optional arguments: [-o OUTPUT] Output chg file name of the cage system;
# Example: prespfit input.molden -o input.chg
```
Note: some cages have very large size. You may need to submit a jobscript to do this. The script support

### pdb4munro: PDB file preparation and mol2 file generation for tleap file

```bash
pdb4munro input [-o OUTPUT] [--chg CHG]
# Required arguments: input: input cage coordinates (pdb format or xyz file)
# Optional arguments: [-o OUTPUT] Output pdb file name of the cage system; [--chg CHG] resp charge file
# Example: pdb4munro input.xyz -o bone.pdb --chg input.chg
```

### seasoning: place one guest molecule inside the cage or multiple counterions around the cage

```bash
seasoning
# seasoning is an interactive script.
```

### munro: generate frcmod file for MOCs

```bash
munro -p PDB [-o OUTPUT]
# Required arguments: input: input cage pdb file
# Optional arguments:[-o OUTPUT]: output frcmod file name;
```

### tleapgen: generate tleap file to system

```bash
tleapgen -p PDB [-o OUTPUT] [--solvent SOLVENT] [-f FRCMOD]
# Required arguments: -p PDB: input pdb file 
# Optional arguments:[-o OUTPUT]: output file name; [-f FRCMOD] frcmod file name; [--solvent SOLVENT] solvent for the system
# (default: water or you can put the lib file name here for custom solvent box. None means no solvent)
# Example: tleapgen -p tastybone.pdb -o tleap.in -f munro.frcmod --solvent None
```

### cntlgen: generate cntl file for alchemical transfer method

```bash
cntlgen [--pdb PDB] [--resname RESNAME] [--out OUT]
# Optional arguments: [--pdb PDB] input system pdb file; [--resname RESNAME] residue name of the guest molecule; [-o OUTPUT] Output cntl file name.
# Example: cntlgen cage_solv.pdb -out bone.pdb -resname GS1
```

```bash
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


## Notes / quirks

- The pipeline expects a Pd-N or Pt-N metal-coordinating cage. Bond detection uses a
  1.90 Å covalent cutoff, which naturally excludes metal–ligand bonds.
- `munro --auto-from-pdb` will run `antechamber` and `parmchk2` on every
  `*_template.pdb` it finds — make sure your `cagepipe` env is active so they're on PATH.
- The bundled `gaff2.dat` ships at `src/cagepipe/data/gaff2.dat`. Override with
  `munro -g /custom/gaff2.dat`.
- Multiwfn redistribution is restricted by its license — that's why it's not in the
  conda env file.
