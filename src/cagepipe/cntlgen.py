#!/usr/bin/env python3
"""cntlgen - generate an ATM ``gus_solv.cntl`` from a solvated cage+guest PDB.

Console-script name: ``cntlgen``.

Reads a solvated PDB (default ``gus_solv.pdb``), classifies atoms, and writes
a cntl file with these fields populated (0-based indices, matching ATM):

  LIGAND_ATOMS / LIGAND_CM_ATOMS   atoms whose residue name == --resname
                                   (default ``GS1``, the guest name placed by
                                   ``filling`` / ``seasoning_inside``)
  RCPT_CM_ATOMS                    cage atoms: everything that is NOT the guest,
                                   not solvent, not a monatomic ion, and not a
                                   ``seasoning`` counterion (BFA..BFZ, BGA..)
  POS_RESTRAINED_ATOMS             receptor atoms whose residue name starts
                                   with a metal prefix (default ``P,M``; e.g.
                                   ``P1``/``P2`` for Pd, ``M1``..``M12``)

All other settings (TEMPERATURES, LAMBDAS, WALL_TIME, restraint constants, ...)
come from a template cntl. A sensible default template is built in; pass
``--template path/to/your.cntl`` to use a custom one.

Usage:
  cntlgen                                   # PDB=gus_solv.pdb, resname=GS1
  cntlgen --pdb gus_solv.pdb --resname GS1
  cntlgen --template my_template.cntl --out gus_solv.cntl
"""
import argparse
import os
import re
import sys


# Residue names to exclude from the receptor.
_WATER = {"HOH", "WAT", "SOL", "TIP", "TIP3", "TIP4", "TIP5",
          "T3P", "T4P", "T5P", "OPC", "OPC3"}
_IONS = {"NA", "CL", "K", "MG", "CA", "ZN",
         "Na+", "Cl-", "K+", "NA+", "CL-"}
# Counterions placed by ``seasoning`` cycle BFA..BFZ, BGA..BGZ, ... (B + [F-Z] + [A-Z])
_COUNTERION_RX = re.compile(r"^B[F-Z][A-Z]$")


def _is_solvent_or_ion(resn):
    return resn in _WATER or resn in _IONS or bool(_COUNTERION_RX.match(resn))


def parse_pdb_atoms(pdb_path):
    """Return list of (index_0based, resname, atom_name) in file order."""
    atoms = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                resn = line[17:20].strip()
                name = line[12:16].strip()
                atoms.append((len(atoms), resn, name))
    return atoms


def _fmt_atom_list(atoms):
    return ", ".join(str(i) for i in atoms)


_FIELD_RX = {
    "LIGAND_ATOMS":         re.compile(r"^(\s*LIGAND_ATOMS\s*=\s*).*$",         re.M),
    "LIGAND_CM_ATOMS":      re.compile(r"^(\s*LIGAND_CM_ATOMS\s*=\s*).*$",      re.M),
    "RCPT_CM_ATOMS":        re.compile(r"^(\s*RCPT_CM_ATOMS\s*=\s*).*$",        re.M),
    "POS_RESTRAINED_ATOMS": re.compile(r"^(\s*POS_RESTRAINED_ATOMS\s*=\s*).*$", re.M),
}


def _substitute(template, mapping):
    """Update atom-list fields in the template. Supports both placeholder
    tokens (``__LIGAND_ATOMS__``) and inline replacement of existing values."""
    out = template
    for key, value in mapping.items():
        token = f"__{key}__"
        if token in out:
            out = out.replace(token, value)
            continue
        rx = _FIELD_RX[key]
        if rx.search(out):
            out = rx.sub(lambda m, v=value: m.group(1) + v, out)
        else:
            sys.stderr.write(f"warning: template has no {key} line; appending\n")
            out += f"\n{key} = {value}\n"
    return out


_DEFAULT_TEMPLATE = """\

#The job transport is the mean in which replicas are executed on GPU devices
#LOCAL_OPENMM is the only job transport system currently supported. Each local GPU is
#managed by a different process using the python multiprocessing module
JOB_TRANSPORT = 'LOCAL_OPENMM'

#The basename of the job. Input amber files are expected to be called <jobname>.prmtop and <jobname>.inpcrd
#The checkpoint file is expected to be called <jobname>_0.xml
BASENAME = 'gus_solv'

#Arrays of thermodynamic states in temperature and alchemical space.
TEMPERATURES = '300'
LAMBDAS =    '0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00'
DIRECTION=   '   1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1'
INTERMEDIATE='   0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    1,    1,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0'
LAMBDA1 =    '0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.10, 0.20, 0.30, 0.40, 0.50, 0.50, 0.40, 0.30, 0.20, 0.10, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00'
LAMBDA2 =    '0.00, 0.10, 0.20, 0.30, 0.40, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.40, 0.30, 0.20, 0.10, 0.00'
ALPHA =      '0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10'
U0 =         '110., 110., 110., 110., 110., 110., 110., 110., 110., 110., 110., 110., 110., 110., 110., 110., 110., 110., 110., 110., 110., 110.'
W0COEFF =    '0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00'

#The displacement vector that brings the ligand from the binding site to a position in the bulk
DISPLACEMENT = '22, 22, 22'

#The position of the binding site restraint. In leg1 it is zero to center it on the receptor
LIGOFFSET = '0., 0., 0.'

#Execution time in minutes
WALL_TIME = 900

#Frequency of replica exchange attempts in seconds
CYCLE_TIME = 10

#Frequency of saving checkpoint files in minutes
CHECKPOINT_TIME = 600

#The nodefile. Each line corresponds to a GPU device.
NODEFILE = 'nodefile'

#Number of replicas to keep in a fast execution queue. It is expressed as a fraction of the number of compute devices.
#With one device a value of 1 here keeps one replica in the queue.
SUBJOBS_BUFFER_SIZE = '1.0'

#MD steps per replica
PRODUCTION_STEPS = '2500'

#frequency of printing information after a replica run. Must be a multiple of PRODUCTION_STEPS
PRNT_FREQUENCY = '2500'

#frequency of saving trajectory frames. Must be a multiple of PRODUCTION_STEPS
TRJ_FREQUENCY = '2500'

#list of ligand atoms.
LIGAND_ATOMS = __LIGAND_ATOMS__

#list of atoms of the ligand that define the centroid of the ligand.
LIGAND_CM_ATOMS = __LIGAND_CM_ATOMS__

#list of atoms of the ligand that define the centroid of the binding site.
RCPT_CM_ATOMS = __RCPT_CM_ATOMS__

#force constant (in kcal/(mol A^2)) and tolerance (in A) of the binding site restraint potential
CM_KF = 15.00
CM_TOL = 8.00

#list of atoms that are restrained followed by the corresponding force constant and tolerance,
#in kcal/(mol A^2) and angstroms, respectively
POS_RESTRAINED_ATOMS = __POS_RESTRAINED_ATOMS__
POSRE_FORCE_CONSTANT = 25.0
POSRE_TOLERANCE = 1.5

#softcore parameters in kcal/mol, acore is dimensionless
UMAX = 200.00
ACORE = 0.062500
UBCORE = 100.0

#thermostat friction coefficient in 1/ps
FRICTION_COEFF = 0.500000

#MD time step in ps
TIME_STEP = 0.002

#GPU platform
OPENMM_PLATFORM = OpenCL

#set to 'yes' to turn on verbose logging
VERBOSE = 'no'
"""


def build_parser():
    ap = argparse.ArgumentParser(
        prog="cntlgen",
        description="Generate an ATM gus_solv.cntl from a solvated PDB.",
    )
    ap.add_argument("--pdb", default="gus_solv.pdb",
                    help="solvated PDB to read (default: gus_solv.pdb)")
    ap.add_argument("--resname", default="GS1",
                    help="guest residue name (default: GS1 — matches the "
                         "default from  `seasoning`)")
    ap.add_argument("--template", default=None,
                    help="path to a cntl template; defaults to the built-in template")
    ap.add_argument("--out", default="gus_solv.cntl",
                    help="output cntl path (default: gus_solv.cntl)")
    ap.add_argument("--metal-prefix", default="P,M",
                    help="comma-separated residue-name prefixes treated as metal "
                         "atoms for POS_RESTRAINED_ATOMS (default: 'P,M')")
    return ap


def main(argv=None):
    args = build_parser().parse_args(argv)

    if not os.path.isfile(args.pdb):
        sys.stderr.write(f"cntlgen: PDB not found: {args.pdb}\n")
        sys.exit(2)

    atoms = parse_pdb_atoms(args.pdb)
    if not atoms:
        sys.stderr.write(f"cntlgen: no ATOM/HETATM records in {args.pdb}\n")
        sys.exit(2)

    metal_prefixes = tuple(p.strip() for p in args.metal_prefix.split(",") if p.strip())

    guest = [i for i, r, _ in atoms if r == args.resname]
    if not guest:
        unique = sorted({r for _, r, _ in atoms})
        sys.stderr.write(
            f"cntlgen: no atoms with residue name '{args.resname}' in {args.pdb}\n"
            f"  unique resnames found: {unique}\n"
        )
        sys.exit(2)

    receptor = [i for i, r, _ in atoms
                if r != args.resname and not _is_solvent_or_ion(r)]
    metals = [i for i, r, _ in atoms
              if r != args.resname
              and not _is_solvent_or_ion(r)
              and r.startswith(metal_prefixes)]

    if args.template:
        with open(args.template) as f:
            template = f.read()
    else:
        template = _DEFAULT_TEMPLATE

    cntl = _substitute(template, {
        "LIGAND_ATOMS":         _fmt_atom_list(guest),
        "LIGAND_CM_ATOMS":      _fmt_atom_list(guest),
        "RCPT_CM_ATOMS":        _fmt_atom_list(receptor),
        "POS_RESTRAINED_ATOMS": _fmt_atom_list(metals),
    })

    with open(args.out, "w") as f:
        f.write(cntl)

    print(f"wrote {args.out}")
    print(f"  guest  ({args.resname}):           {len(guest):>4} atoms  "
          f"[{guest[0]}..{guest[-1]}]")
    if receptor:
        print(f"  receptor (cage):           {len(receptor):>4} atoms  "
              f"[{receptor[0]}..{receptor[-1]}]")
    else:
        print("  receptor (cage):           EMPTY  (check residue names)")
    if metals:
        print(f"  metals (POS_RESTRAINED):   {len(metals):>4} atoms  {metals}")
    else:
        print("  metals (POS_RESTRAINED):   none detected "
              f"(prefixes={list(metal_prefixes)})")


if __name__ == "__main__":
    main()
