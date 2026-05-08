#!/usr/bin/env python3
"""
gaff_typing.py - Run antechamber + parmchk2 on ligand template PDBs.

Wraps the two AmberTools commands the user used to call by hand:
  antechamber -fi pdb -fo mol2 -i <pdb> -o <mol2> -pf y -nc 0
  parmchk2    -i <mol2> -o <frcmod> -f mol2

Antechamber is run inside a tmp cwd so its side files (ANTECHAMBER_AC.AC,
NEWPDB.PDB, ATOMTYPE.INF, sqm.in, sqm.out, ...) don't pollute the project.

Usage as a library:
    from gaff_typing import generate_for_templates
    results = generate_for_templates(["LAtemp1_template.pdb", "LBtemp2_template.pdb"])
    # results = [(mol2_path, frcmod_path, prefix), ...]

Usage as a script:
    python gaff_typing.py LAtemp1_template.pdb LBtemp2_template.pdb -nc 0
"""
import argparse
import glob
import os
import re
import shutil
import subprocess
import sys
import tempfile

# Search order: PATH first (so the active conda env wins), then known fallbacks.
ANTECHAMBER_CANDIDATES = [
    "antechamber",
    "/home/gridsan/ywang6/.conda/envs/metallicious/bin/antechamber",
    "/home/gridsan/ywang6/.conda/envs/AmberTools25/bin/antechamber",
]
PARMCHK2_CANDIDATES = [
    "parmchk2",
    "/home/gridsan/ywang6/.conda/envs/metallicious/bin/parmchk2",
    "/home/gridsan/ywang6/.conda/envs/AmberTools25/bin/parmchk2",
]


def _which(candidates):
    for p in candidates:
        if os.path.isfile(p) and os.access(p, os.X_OK):
            return p
        located = shutil.which(p)
        if located:
            return located
    return None


def find_antechamber():
    return _which(ANTECHAMBER_CANDIDATES)


def find_parmchk2():
    return _which(PARMCHK2_CANDIDATES)


def derive_prefix(pdb_path):
    """
    Heuristic prefix from a template filename.
    Matches what pdb4munro.py / munro.py expect (LA, LB, LC, LD, ... or just L).

    Examples:
        LAtemp1_template.pdb -> LA
        Ltemp1_template.pdb  -> L
        L1_template.pdb      -> L
        cage_part.pdb        -> cage_part
    """
    base = os.path.splitext(os.path.basename(pdb_path))[0]
    m = re.match(r"^(L[A-Za-z]?)(?:\d|temp|_).*$", base)
    if m:
        return m.group(1)
    return base


def run_antechamber(pdb_path, mol2_path, net_charge=0, atom_type="gaff",
                    extra_args=None, verbose=False):
    """Run antechamber pdb -> mol2."""
    exe = find_antechamber()
    if exe is None:
        raise RuntimeError(
            "antechamber not found on PATH. Activate an AmberTools conda env "
            "(e.g. `conda activate metallicious`) or update ANTECHAMBER_CANDIDATES."
        )

    pdb_abs = os.path.abspath(pdb_path)
    mol2_abs = os.path.abspath(mol2_path)
    cmd = [
        exe,
        "-i", pdb_abs, "-fi", "pdb",
        "-o", mol2_abs, "-fo", "mol2",
        "-pf", "y",
        "-nc", str(net_charge),
        "-at", atom_type,
    ]
    if extra_args:
        cmd.extend(extra_args)

    work = tempfile.mkdtemp(prefix="ante_")
    try:
        if verbose:
            print(f"    [antechamber] {' '.join(cmd)}")
        result = subprocess.run(cmd, cwd=work, capture_output=True, text=True)
        if result.returncode != 0 or not os.path.isfile(mol2_abs):
            raise RuntimeError(
                f"antechamber failed (rc={result.returncode}) on {pdb_path}\n"
                f"--- stdout ---\n{result.stdout}\n"
                f"--- stderr ---\n{result.stderr}"
            )
    finally:
        shutil.rmtree(work, ignore_errors=True)
    return mol2_abs


def run_parmchk2(mol2_path, frcmod_path, atom_type="gaff", verbose=False):
    """Run parmchk2 mol2 -> frcmod."""
    exe = find_parmchk2()
    if exe is None:
        raise RuntimeError(
            "parmchk2 not found on PATH. Activate an AmberTools conda env."
        )
    mol2_abs = os.path.abspath(mol2_path)
    frcmod_abs = os.path.abspath(frcmod_path)
    cmd = [
        exe,
        "-i", mol2_abs, "-o", frcmod_abs,
        "-f", "mol2",
        "-s", atom_type,
    ]
    if verbose:
        print(f"    [parmchk2]    {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0 or not os.path.isfile(frcmod_abs):
        raise RuntimeError(
            f"parmchk2 failed (rc={result.returncode}) on {mol2_path}\n"
            f"--- stdout ---\n{result.stdout}\n"
            f"--- stderr ---\n{result.stderr}"
        )
    return frcmod_abs


def _is_cache_fresh(src, *outputs):
    if not all(os.path.isfile(o) for o in outputs):
        return False
    src_mtime = os.path.getmtime(src)
    return all(os.path.getmtime(o) >= src_mtime for o in outputs)


def generate(template_pdb, mol2_out=None, frcmod_out=None,
             net_charge=0, atom_type="gaff", force=False, verbose=False):
    """
    Run antechamber + parmchk2 on one template PDB. If both outputs are
    newer than the input PDB and `force` is False, this is a no-op.
    Returns (mol2_path, frcmod_path).
    """
    if not os.path.isfile(template_pdb):
        raise FileNotFoundError(template_pdb)
    base = os.path.splitext(template_pdb)[0]
    if mol2_out is None:
        mol2_out = base + ".mol2"
    if frcmod_out is None:
        frcmod_out = base + ".frcmod"

    if not force and _is_cache_fresh(template_pdb, mol2_out, frcmod_out):
        if verbose:
            print(f"  [cache] {template_pdb} -> {mol2_out}, {frcmod_out}")
        return mol2_out, frcmod_out

    print(f"  antechamber: {template_pdb} -> {mol2_out}")
    run_antechamber(template_pdb, mol2_out,
                    net_charge=net_charge, atom_type=atom_type, verbose=verbose)
    print(f"  parmchk2:    {mol2_out} -> {frcmod_out}")
    run_parmchk2(mol2_out, frcmod_out, atom_type=atom_type, verbose=verbose)
    return mol2_out, frcmod_out


def generate_for_templates(template_pdbs, output_dir=".",
                           net_charge=0, charges_by_prefix=None,
                           atom_type="gaff", force=False, verbose=False):
    """
    Process a list of template PDB files. For each, produce <prefix>.mol2 and
    <prefix>.frcmod under `output_dir`, where prefix is derived from the
    filename (see derive_prefix).

    `charges_by_prefix` lets you override net charge per template, e.g.
        {"LA": 0, "LB": -1}

    Returns: list of (mol2_path, frcmod_path, prefix).
    """
    os.makedirs(output_dir, exist_ok=True)
    charges_by_prefix = charges_by_prefix or {}
    out = []
    for pdb in template_pdbs:
        prefix = derive_prefix(pdb)
        nc = charges_by_prefix.get(prefix, net_charge)
        mol2 = os.path.join(output_dir, f"{prefix}.mol2")
        frcmod = os.path.join(output_dir, f"{prefix}.frcmod")
        m, f = generate(pdb, mol2, frcmod,
                        net_charge=nc, atom_type=atom_type,
                        force=force, verbose=verbose)
        out.append((m, f, prefix))
    return out


def discover_templates(pattern="*template*.pdb", search_dir="."):
    """Find candidate template PDB files by glob. Sorted for determinism."""
    return sorted(glob.glob(os.path.join(search_dir, pattern)))


def main():
    parser = argparse.ArgumentParser(
        description="Wrap antechamber + parmchk2 over one or more ligand template PDBs.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process every *_template.pdb in current dir, net charge 0:
  python gaff_typing.py --auto

  # Explicit templates:
  python gaff_typing.py LAtemp1_template.pdb LBtemp2_template.pdb -nc 0

  # Per-template charges:
  python gaff_typing.py LA*.pdb LB*.pdb --charge LA:0 LB:-1
        """,
    )
    parser.add_argument("templates", nargs="*",
                        help="Template PDB files (or use --auto)")
    parser.add_argument("--auto", action="store_true",
                        help="Glob *template*.pdb in current directory")
    parser.add_argument("-nc", "--net-charge", type=int, default=0,
                        help="Default net charge (default: 0)")
    parser.add_argument("--charge", nargs="+", default=[],
                        metavar="PREFIX:CHARGE",
                        help="Per-prefix charges, e.g. LA:0 LB:-1")
    parser.add_argument("-at", "--atom-type", default="gaff",
                        choices=["gaff", "gaff2", "amber", "bcc", "sybyl", "amber14sb"],
                        help="Atom-type style for antechamber/parmchk2 (default: gaff)")
    parser.add_argument("-d", "--output-dir", default=".",
                        help="Where to write .mol2 / .frcmod (default: cwd)")
    parser.add_argument("--force", action="store_true",
                        help="Re-run even if cached outputs are newer than the PDB")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()

    templates = list(args.templates)
    if args.auto or not templates:
        if not args.auto and not templates:
            print("No templates given; use --auto to glob *template*.pdb",
                  file=sys.stderr)
            sys.exit(2)
        templates.extend(discover_templates())
        # de-dup, preserving order
        seen = set()
        templates = [t for t in templates if not (t in seen or seen.add(t))]
    if not templates:
        print("No template PDB files found.", file=sys.stderr)
        sys.exit(1)

    charges_by_prefix = {}
    for spec in args.charge:
        if ":" not in spec:
            print(f"--charge expects PREFIX:VALUE, got {spec!r}", file=sys.stderr)
            sys.exit(2)
        k, v = spec.split(":", 1)
        charges_by_prefix[k] = int(v)

    print("=" * 70)
    print(f"GAFF typing for {len(templates)} template(s)")
    print(f"  antechamber: {find_antechamber()}")
    print(f"  parmchk2:    {find_parmchk2()}")
    print("=" * 70)

    results = generate_for_templates(
        templates,
        output_dir=args.output_dir,
        net_charge=args.net_charge,
        charges_by_prefix=charges_by_prefix,
        atom_type=args.atom_type,
        force=args.force,
        verbose=args.verbose,
    )

    print("\nDone. Generated:")
    for mol2, frcmod, prefix in results:
        print(f"  [{prefix}] {mol2}, {frcmod}")


if __name__ == "__main__":
    main()
