#!/usr/bin/env python3
"""
chgass.py - Per-residue MOL2 generation from pre-charged ligand & metal templates.

Used downstream of `respfit.py --free ... --cluster ...` (differential mode)
for M12L24-style cages where whole-cage RESP is intractable. Inputs are:
  - a cage PDB whose residues have already been split/named by pdb4munro.py
  - one or more ligand-template mol2 files carrying RESP (or RESP-corrected)
    charges, Sybyl types, and bonds — typically `L_corr.mol2` from
    `respfit.py` differential mode
  - one or more metal-template mol2 files carrying the metal RESP charge

For every residue in the PDB, this script:
  1. picks a matching ligand template by graph isomorphism (or a metal template
     by element);
  2. emits <resname>.mol2 with template-derived atom names, types, bonds, and
     charges, and PDB-derived coordinates;
  3. once all residues are processed, applies MUNRO atom types: M1/M2/... on
     metals and Y1/Y2/... on coordinating nitrogens (within METAL_BOND_CUTOFF
     of any metal), so the resulting mol2 set is consistent with tleapgen.py's
     addAtomTypes section.

Usage:
    python chgass.py cage.pdb \
        --ligand-template L_corr.mol2 [--ligand-template L2_corr.mol2 ...] \
        --metal-template  P_corr.mol2 [--metal-template Pt_corr.mol2 ...]

Writes one MOL2 per residue in the cwd (LA1.mol2, LA2.mol2, P1.mol2, ...).
"""
import argparse
import os
import sys
from collections import defaultdict

try:
    from . import mol2gen_helper as mgh          # package mode
except ImportError:
    import mol2gen_helper as mgh                  # script-mode fallback

METALS = mgh.METALS


# --------------------- Template loading ---------------------

def load_ligand_template(mol2_path):
    """Parse a fully-charged ligand template mol2.

    Returns dict with: atoms (list of mgh.Atom with .charge/.sybyl_type set),
    adj (build_adjacency), bonds (list of (i, j, btype)), and the template
    file path for diagnostics."""
    atoms = []
    bonds = []
    id_map = {}
    section = None
    with open(mol2_path) as f:
        for line in f:
            s = line.strip()
            if s.startswith("@<TRIPOS>"):
                section = s
                continue
            if not s:
                continue
            p = s.split()
            if section == "@<TRIPOS>ATOM" and len(p) >= 6:
                file_id = p[0]
                name = p[1]
                x, y, z = float(p[2]), float(p[3]), float(p[4])
                sybyl = p[5]
                charge = float(p[8]) if len(p) >= 9 else 0.0
                idx = len(atoms)
                id_map[file_id] = idx
                a = mgh.Atom(name, x, y, z, "TMP", 1)
                a.sybyl_type = sybyl
                a.charge = charge
                atoms.append(a)
            elif section == "@<TRIPOS>BOND" and len(p) >= 4:
                try:
                    bonds.append((id_map[p[1]], id_map[p[2]], p[3]))
                except KeyError:
                    continue
    if not atoms:
        raise RuntimeError(f"No atoms parsed from {mol2_path}")
    adj = mgh.build_adjacency(atoms, use_covalent_radii=True)
    return {
        'path': mol2_path,
        'atoms': atoms,
        'adj': adj,
        'bonds': bonds,
    }


def load_metal_template(mol2_path):
    """Parse a metal template mol2 (one heavy atom + a charge).

    Returns (element, charge, sybyl_type). The sybyl_type from the template
    is recorded but typically overwritten by apply_munro_types later."""
    section = None
    with open(mol2_path) as f:
        for line in f:
            s = line.strip()
            if s.startswith("@<TRIPOS>"):
                section = s
                continue
            if section == "@<TRIPOS>ATOM" and s:
                p = s.split()
                if len(p) >= 6:
                    name = p[1]
                    sybyl = p[5]
                    charge = float(p[8]) if len(p) >= 9 else 0.0
                    a = mgh.Atom(name, p[2], p[3], p[4], "TMP", 1)
                    if a.element not in METALS:
                        continue
                    return a.element, charge, sybyl
    raise RuntimeError(f"No metal atom found in {mol2_path}")


# --------------------- Per-residue assignment ---------------------

def assign_ligand(resid, pdb_atoms, templates, debug=False):
    """Match a ligand residue against the list of templates by isomorphism.

    Returns (template, mapping) where mapping is template_idx -> pdb_idx,
    or (None, None) if no template fits."""
    pdb_adj = mgh.build_adjacency(pdb_atoms, use_covalent_radii=True)
    for tpl in templates:
        if len(tpl['atoms']) != len(pdb_atoms):
            continue
        m = mgh.solve_isomorphism(tpl['atoms'], tpl['adj'],
                                  pdb_atoms, pdb_adj, debug=debug)
        if m is not None:
            return tpl, m
    return None, None


def build_residue_atoms(pdb_atoms, template, mapping, resname, resid):
    """Reorder PDB atoms into template order; copy template name/type/charge,
    keep PDB coordinates."""
    out = []
    for t_idx in range(len(template['atoms'])):
        p_idx = mapping[t_idx]
        pdb_atom = pdb_atoms[p_idx]
        t_atom = template['atoms'][t_idx]
        new_atom = mgh.Atom(t_atom.name, pdb_atom.x, pdb_atom.y, pdb_atom.z,
                            resname, resid, element=pdb_atom.element)
        new_atom.sybyl_type = t_atom.sybyl_type
        new_atom.charge = t_atom.charge
        out.append(new_atom)
    return out


# --------------------- Main ---------------------

def main():
    ap = argparse.ArgumentParser(
        description="Assign per-residue MOL2 files from a cage PDB and "
                    "pre-charged ligand/metal templates.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Single ligand type, single metal type (typical M12L24):
  python chgass.py bone.pdb \\
      --ligand-template L_corr.mol2 --metal-template P_corr.mol2

  # Mixed-ligand cage (e.g., LA + LB), Pd cage:
  python chgass.py bone.pdb \\
      --ligand-template LA_corr.mol2 --ligand-template LB_corr.mol2 \\
      --metal-template P_corr.mol2
        """,
    )
    ap.add_argument("pdb", help="Cage PDB (residues already split by pdb4munro.py)")
    ap.add_argument("--ligand-template", action="append", required=True,
                    help="Pre-charged ligand template mol2 (repeatable)")
    ap.add_argument("--metal-template", action="append", required=True,
                    help="Pre-charged metal template mol2 (repeatable)")
    ap.add_argument("--debug", action="store_true")
    args = ap.parse_args()

    print("=" * 70)
    print("CHGASS - per-residue MOL2 from pre-charged templates")
    print("=" * 70)

    if not os.path.isfile(args.pdb):
        print(f"ERROR: PDB not found: {args.pdb}")
        sys.exit(1)

    # 1. Load ligand templates
    print(f"\n1. Loading {len(args.ligand_template)} ligand template(s)")
    ligand_templates = []
    for path in args.ligand_template:
        if not os.path.isfile(path):
            print(f"   ERROR: template not found: {path}")
            sys.exit(1)
        tpl = load_ligand_template(path)
        ligand_templates.append(tpl)
        net = sum(a.charge for a in tpl['atoms'])
        print(f"   {os.path.basename(path)}: {len(tpl['atoms'])} atoms, "
              f"{len(tpl['bonds'])} bonds, net charge {net:+.4f}")

    # 2. Load metal templates -> dict element -> (charge, sybyl)
    print(f"\n2. Loading {len(args.metal_template)} metal template(s)")
    metal_charges = {}
    for path in args.metal_template:
        if not os.path.isfile(path):
            print(f"   ERROR: template not found: {path}")
            sys.exit(1)
        el, q, sybyl = load_metal_template(path)
        metal_charges[el] = (q, sybyl)
        print(f"   {os.path.basename(path)}: {el}, charge {q:+.4f}")

    # 3. Read cage PDB
    print(f"\n3. Reading cage PDB: {args.pdb}")
    residues = mgh.read_pdb(args.pdb)
    print(f"   {len(residues)} residue(s) to process")

    # 4. Process each residue
    print("\n4. Assigning residues:")
    all_atoms = []
    files_to_write = []
    skipped = []
    for resid in sorted(residues.keys()):
        atoms = residues[resid]
        resname = atoms[0].resname

        # Single-atom metal residue
        if len(atoms) == 1 and atoms[0].element in METALS:
            a = atoms[0]
            a.name = a.element
            entry = metal_charges.get(a.element)
            if entry is None:
                print(f"   {resname:>5} (resid {resid}): metal {a.element} "
                      f"with NO template - skipping")
                skipped.append((resname, resid, "no metal template"))
                continue
            q, sybyl = entry
            a.charge = q
            a.sybyl_type = sybyl
            all_atoms.append(a)
            files_to_write.append((f"{resname}.mol2", resname, [a], []))
            print(f"   {resname:>5} (resid {resid}): metal {a.element}, "
                  f"q = {q:+.4f}")
            continue

        # Ligand: match by isomorphism
        tpl, mapping = assign_ligand(resid, atoms, ligand_templates,
                                     debug=args.debug)
        if tpl is None:
            print(f"   {resname:>5} (resid {resid}): no ligand template matches "
                  f"({len(atoms)} atoms) - skipping")
            skipped.append((resname, resid, "no ligand template match"))
            continue
        new_atoms = build_residue_atoms(atoms, tpl, mapping, resname, resid)
        all_atoms.extend(new_atoms)
        files_to_write.append((f"{resname}.mol2", resname, new_atoms, tpl['bonds']))
        print(f"   {resname:>5} (resid {resid}): matched "
              f"{os.path.basename(tpl['path'])}, {len(new_atoms)} atoms")

    # 5. Apply MUNRO types (M1/M2/... on metals; Y1/Y2/... on coord N)
    print("\n5. Applying MUNRO atom types")
    mgh.apply_munro_types(all_atoms)

    # 6. Write per-residue MOL2s
    print("\n6. Writing per-residue MOL2 files")
    for fname, rname, ats, bnds in files_to_write:
        mgh.write_mol2(fname, rname, ats, bnds)
        print(f"   -> {fname}")

    # 7. Summary
    print("\n" + "=" * 70)
    print("Summary")
    print("=" * 70)
    print(f"  Wrote: {len(files_to_write)} MOL2 file(s)")
    if skipped:
        print(f"  Skipped: {len(skipped)} residue(s)")
        for r in skipped:
            print(f"    {r[0]} (resid {r[1]}): {r[2]}")
    print("=" * 70)

    return 0 if not skipped else 2


if __name__ == "__main__":
    sys.exit(main())
