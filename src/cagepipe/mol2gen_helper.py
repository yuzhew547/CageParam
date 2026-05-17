#!/usr/bin/env python3
"""
mol2gen_helper.py - Multi-template MOL2 generator with CHG support.

Importable module (no CLI). Driven by pdb4munro.py once it has produced
template MOL2 files via antechamber. Public entrypoint:

    generate_mol2_files(pdb_path, chg_path, template_mol2s, debug=False)

Reads:
  - pdb_path             : bone.pdb (residues with template-ordered atom names)
  - chg_path             : CHG file (ele x y z q per atom) with RESP charges
  - template_mol2s       : list of antechamber-produced template MOL2 files
                           carrying Sybyl atom types and explicit bond orders

Writes one MOL2 per residue in the cwd (e.g. LA1.mol2, P1.mol2, ...).
"""
import math
import os
from collections import defaultdict

import numpy as np

# ================= CONFIGURATION =================
BOND_CUTOFF = 1.90
METAL_BOND_CUTOFF = 2.5
# =================================================

COVALENT_RADII = {
    'H': 0.31, 'C': 0.76, 'N': 0.71, 'O': 0.66, 'F': 0.57,
    'P': 1.07, 'S': 1.05, 'CL': 1.02, 'BR': 1.20, 'I': 1.39,
    'B': 0.84, 'SI': 1.11, 'PD': 1.39, 'PT': 1.36, 'AU': 1.36, 'AG': 1.45,
}

METALS = {'PD', 'PT', 'AU', 'AG', 'FE', 'CO', 'NI', 'CU', 'ZN', 'RU', 'RH', 'IR'}


def get_bond_cutoff(el1, el2, tolerance=0.40):
    """Sum of Cordero covalent radii + tolerance.
    Tolerance was 0.45 A but caught a real non-bonded C-C contact at 1.941 A
    in tightly-packed amide cages, breaking template isomorphism. 0.40 A keeps
    all standard covalent bonds (longest typical organic single bond is C-I at
    ~2.14 A, well under any C-I cutoff) while rejecting close non-bonded contacts.
    """
    r1 = COVALENT_RADII.get(el1, 0.77)
    r2 = COVALENT_RADII.get(el2, 0.77)
    return r1 + r2 + tolerance


# --- 1. UTILITIES ---
def clean(text):
    if not text:
        return ""
    return text.strip().upper()


def get_unique_2char_code(index, prefixes):
    chars = [str(i) for i in range(0, 10)] + \
            [chr(i) for i in range(ord('A'), ord('Z') + 1)] + \
            [chr(i) for i in range(ord('a'), ord('z') + 1)]
    prefix_idx = index // len(chars)
    char_idx = index % len(chars)
    if prefix_idx >= len(prefixes):
        return "XX"
    return f"{prefixes[prefix_idx]}{chars[char_idx]}"


def dist_sq(a1, a2):
    return (a1.x - a2.x) ** 2 + (a1.y - a2.y) ** 2 + (a1.z - a2.z) ** 2


def dist(a1, a2):
    return math.sqrt(dist_sq(a1, a2))


# --- 2. DATA CLASSES ---
class Atom:
    def __init__(self, name, x, y, z, resname, resid, element=None):
        self.name = clean(name)
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)
        self.resname = clean(resname)
        self.resid = resid
        self.charge = 0.0
        self.sybyl_type = "Du"

        if element:
            self.element = element.upper()
        else:
            letters = "".join([c for c in self.name if c.isalpha()])
            if letters:
                if letters[:2] in ['BR', 'CL', 'PD', 'PT', 'AU', 'AG']:
                    self.element = letters[:2]
                else:
                    self.element = letters[0]
            else:
                self.element = "X"


# --- 3. GRAPH ENGINE (ISOMORPHISM) ---
def build_adjacency(atoms, use_covalent_radii=True):
    adj = defaultdict(list)
    n = len(atoms)
    for i in range(n):
        for j in range(i + 1, n):
            cutoff = get_bond_cutoff(atoms[i].element, atoms[j].element) \
                if use_covalent_radii else BOND_CUTOFF
            d = math.sqrt(dist_sq(atoms[i], atoms[j]))
            if d < cutoff:
                adj[i].append(j)
                adj[j].append(i)
    return adj


def get_extended_signature(idx, atoms, adj, depth=2):
    el = atoms[idx].element
    neighbors = adj[idx]
    n_els = tuple(sorted([atoms[n].element for n in neighbors]))
    if depth >= 2:
        second_level = []
        for n in neighbors:
            second_n_els = tuple(sorted(
                [atoms[nn].element for nn in adj[n] if nn != idx]
            ))
            second_level.append(second_n_els)
        second_level = tuple(sorted(second_level))
        return (el, len(neighbors), n_els, second_level)
    return (el, len(neighbors), n_els)


def solve_isomorphism(templ_atoms, templ_adj, target_atoms, target_adj, debug=False):
    if len(templ_atoms) != len(target_atoms):
        return None

    t_sigs = [get_extended_signature(i, templ_atoms, templ_adj) for i in range(len(templ_atoms))]
    tgt_sigs = [get_extended_signature(i, target_atoms, target_adj) for i in range(len(target_atoms))]

    mapping = {}
    used_targets = set()

    sig_counts = defaultdict(int)
    for sig in t_sigs:
        sig_counts[sig] += 1
    order = sorted(range(len(templ_atoms)), key=lambda i: sig_counts[t_sigs[i]])

    def backtrack(pos):
        if pos == len(templ_atoms):
            return True
        t_idx = order[pos]
        my_sig = t_sigs[t_idx]
        for cand_idx in range(len(target_atoms)):
            if cand_idx in used_targets:
                continue
            if tgt_sigs[cand_idx] != my_sig:
                continue
            valid = True
            for t_neighbor in templ_adj[t_idx]:
                if t_neighbor in mapping:
                    target_neighbor = mapping[t_neighbor]
                    if target_neighbor not in target_adj[cand_idx]:
                        valid = False
                        break
            if valid:
                mapping[t_idx] = cand_idx
                used_targets.add(cand_idx)
                if backtrack(pos + 1):
                    return True
                del mapping[t_idx]
                used_targets.remove(cand_idx)
        return False

    return mapping if backtrack(0) else None


# --- 4. READERS ---
def read_template(filepath):
    print(f"  Reading template: {filepath}")
    atoms = []
    bonds = []
    atom_types = {}

    with open(filepath, 'r') as f:
        lines = f.readlines()

    in_atom = False
    in_bond = False
    id_map = {}

    for line in lines:
        if line.startswith("@<TRIPOS>ATOM"):
            in_atom = True
            in_bond = False
            continue
        if line.startswith("@<TRIPOS>BOND"):
            in_atom = False
            in_bond = True
            continue
        if line.startswith("@<TRIPOS>"):
            in_atom = False
            in_bond = False
            continue

        parts = line.split()
        if not parts:
            continue

        if in_atom and len(parts) >= 6:
            idx = len(atoms)
            file_id = parts[0]
            name = parts[1]
            x, y, z = parts[2], parts[3], parts[4]
            sybyl_type = parts[5]
            id_map[file_id] = idx
            atom = Atom(name, x, y, z, "TMP", 1)
            atoms.append(atom)
            atom_types[clean(name)] = sybyl_type

        if in_bond and len(parts) >= 4:
            try:
                a1 = id_map[parts[1]]
                a2 = id_map[parts[2]]
                btype = parts[3]
                bonds.append((a1, a2, btype))
            except Exception:
                pass

    # Rebuild adjacency from coordinates with the same heuristic used for
    # CHG fragments. Ensures isomorphism on consistent bond perception, even
    # when the template MOL2 was produced via the OpenBabel fallback (which
    # has been observed to drop legitimate C-O ester bonds at ~1.46 A).
    adj_coords = build_adjacency(atoms, use_covalent_radii=True)
    coord_bonds = set()
    for i, nbrs in adj_coords.items():
        for j in nbrs:
            coord_bonds.add(tuple(sorted((i, j))))
    file_bonds = {tuple(sorted((a1, a2))) for a1, a2, _ in bonds}
    btype_lookup = {tuple(sorted((a1, a2))): bt for a1, a2, bt in bonds}

    merged_bonds = []
    for pair in sorted(coord_bonds | file_bonds):
        bt = btype_lookup.get(pair, "1")
        merged_bonds.append((pair[0], pair[1], bt))

    adj = defaultdict(list)
    for a1, a2, _ in merged_bonds:
        adj[a1].append(a2)
        adj[a2].append(a1)

    print(f"    -> {len(atoms)} atoms, {len(merged_bonds)} bonds "
          f"(file: {len(bonds)}, coord-derived: {len(coord_bonds)})")
    return atoms, adj, merged_bonds, atom_types


def read_pdb(filepath):
    print(f"\n  Reading PDB: {filepath}")
    residues = defaultdict(list)
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                name = line[12:16]
                resname = line[17:20]
                try:
                    resid = int(line[22:26])
                except ValueError:
                    resid = 1
                x, y, z = line[30:38], line[38:46], line[46:54]
                el = line[76:78].strip() if len(line) > 77 else None
                if clean(resname) not in ["WAT", "HOH", "TIP3", "TIP"]:
                    residues[resid].append(Atom(name, x, y, z, resname, resid, el))
    print(f"    -> {len(residues)} residues")
    return residues


def cluster_chg_into_molecules(filepath, exclude_metals=True):
    print(f"\n  Reading CHG file: {filepath}")
    chg_atoms = []
    with open(filepath, 'r') as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 5:
                try:
                    chg_atoms.append({
                        'element': parts[0].upper(),
                        'x': float(parts[1]),
                        'y': float(parts[2]),
                        'z': float(parts[3]),
                        'charge': float(parts[4]),
                    })
                except (ValueError, IndexError):
                    continue
    print(f"    -> Loaded {len(chg_atoms)} atoms from CHG")

    metal_atoms = []
    organic_atoms = []
    for atom in chg_atoms:
        if exclude_metals and atom['element'] in METALS:
            metal_atoms.append(atom)
        else:
            organic_atoms.append(atom)
    if exclude_metals:
        print(f"    -> Separated {len(metal_atoms)} metal atoms, "
              f"{len(organic_atoms)} organic atoms")

    molecules = []
    n = len(organic_atoms)
    visited = [False] * n

    def dfs(idx, current_mol):
        visited[idx] = True
        current_mol.append(organic_atoms[idx])
        for j in range(n):
            if not visited[j]:
                el1 = organic_atoms[idx]['element']
                el2 = organic_atoms[j]['element']
                cutoff = get_bond_cutoff(el1, el2)
                dx = organic_atoms[idx]['x'] - organic_atoms[j]['x']
                dy = organic_atoms[idx]['y'] - organic_atoms[j]['y']
                dz = organic_atoms[idx]['z'] - organic_atoms[j]['z']
                d = math.sqrt(dx * dx + dy * dy + dz * dz)
                if d < cutoff:
                    dfs(j, current_mol)

    for i in range(n):
        if not visited[i]:
            current_mol = []
            dfs(i, current_mol)
            molecules.append(current_mol)

    for metal in metal_atoms:
        molecules.append([metal])

    print(f"    -> Clustered into {len(molecules)} molecules")
    return molecules


# --- 5. TEMPLATE LIBRARY ---
class TemplateLibrary:
    def __init__(self):
        self.templates = []

    def add_template(self, mol2_file, resname_prefix=None):
        atoms, adj, bonds, types = read_template(mol2_file)
        if resname_prefix is None:
            basename = os.path.splitext(os.path.basename(mol2_file))[0]
            if basename.startswith('L'):
                if len(basename) >= 2 and basename[1].isalpha():
                    resname_prefix = basename[:2]
                else:
                    resname_prefix = basename[:1]
            else:
                resname_prefix = "L"
        self.templates.append({
            'prefix': resname_prefix,
            'atoms': atoms,
            'adj': adj,
            'bonds': bonds,
            'types': types,
            'mol2_file': mol2_file,
        })
        print(f"    Added template '{resname_prefix}' ({len(atoms)} atoms)")

    def match_pdb_to_template(self, pdb_atoms, debug=False):
        pdb_adj = build_adjacency(pdb_atoms, use_covalent_radii=True)
        for template in self.templates:
            if len(template['atoms']) != len(pdb_atoms):
                continue
            mapping = solve_isomorphism(
                template['atoms'], template['adj'],
                pdb_atoms, pdb_adj, debug=debug,
            )
            if mapping is not None:
                return template, mapping
        return None, None

    def match_chg_to_template(self, chg_mol, debug=False):
        chg_atoms = [Atom(a['element'], a['x'], a['y'], a['z'], "CHG", 1, a['element'])
                     for a in chg_mol]
        chg_adj = build_adjacency(chg_atoms, use_covalent_radii=True)
        for template in self.templates:
            if len(template['atoms']) != len(chg_atoms):
                continue
            mapping = solve_isomorphism(
                template['atoms'], template['adj'],
                chg_atoms, chg_adj, debug=debug,
            )
            if mapping is not None:
                return template, mapping
        return None, None


# --- 6. CHARGE ASSIGNMENT ---
def build_charge_maps(chg_molecules, template_library, debug=False):
    print("\n3. Validating CHG molecules against templates:")

    charges_by_template = defaultdict(lambda: defaultdict(list))
    metal_charges = defaultdict(list)

    for i, chg_mol in enumerate(chg_molecules):
        if len(chg_mol) == 1:
            element = chg_mol[0]['element']
            if element in METALS:
                metal_charges[element].append(chg_mol[0]['charge'])
                print(f"  Molecule {i+1}: Metal ({element})")
            continue

        template, mapping = template_library.match_chg_to_template(chg_mol, debug=debug)
        if template is not None:
            print(f"  Molecule {i+1}: Matched template '{template['prefix']}' "
                  f"({len(chg_mol)} atoms)")
            for t_idx in range(len(template['atoms'])):
                c_idx = mapping[t_idx]
                tpl_name = template['atoms'][t_idx].name
                charge = chg_mol[c_idx]['charge']
                charges_by_template[template['prefix']][tpl_name].append(charge)
        else:
            print(f"  Molecule {i+1}: No template match ({len(chg_mol)} atoms)")

    charge_maps = {}
    for template_prefix, atom_charges in charges_by_template.items():
        avg_charges = {}
        print(f"\n  Average charges for template '{template_prefix}':")
        print(f"    {'Atom':<6} {'Avg Charge':>10} {'Std Dev':>10} {'Count':>6}")
        print(f"    {'-'*6} {'-'*10} {'-'*10} {'-'*6}")
        for atom_name in sorted(atom_charges.keys()):
            charges = atom_charges[atom_name]
            avg = float(np.mean(charges))
            std = float(np.std(charges)) if len(charges) > 1 else 0.0
            avg_charges[atom_name] = avg
            print(f"    {atom_name:<6} {avg:>10.4f} {std:>10.4f} {len(charges):>6}")
        charge_maps[template_prefix] = avg_charges

    metal_charge_map = {}
    if metal_charges:
        print("\n  Metal charges:")
        for element, charges in metal_charges.items():
            avg = float(np.mean(charges))
            metal_charge_map[element] = avg
            print(f"    {element}: {avg:.6f} (from {len(charges)} atoms)")

    return charge_maps, metal_charge_map


# --- 7. PROCESSING ---
def apply_munro_types(all_atoms):
    print("\n  Applying MUNRO types (M, Y codes)...")
    metals = [a for a in all_atoms if a.element in METALS]
    metals.sort(key=lambda a: (a.resid, a.name))
    for i, m in enumerate(metals):
        m.sybyl_type = get_unique_2char_code(i, "MNO")
    print(f"    -> {len(metals)} metal atoms")

    nitrogens = [a for a in all_atoms if a.element == "N"]
    nitrogens.sort(key=lambda a: (a.resid, a.name))
    cutoff_sq = METAL_BOND_CUTOFF ** 2
    y_cnt = 0
    for n in nitrogens:
        connected = any(dist_sq(n, m) < cutoff_sq for m in metals)
        if connected:
            n.sybyl_type = get_unique_2char_code(y_cnt, "YZWVU")
            y_cnt += 1
    print(f"    -> {y_cnt} coordinating nitrogens")


def write_mol2(filename, resname, atoms, bonds):
    with open(filename, 'w') as f:
        f.write(f"@<TRIPOS>MOLECULE\n{resname}\n")
        f.write(f"{len(atoms)} {len(bonds)} 1 0 0\n")
        f.write("SMALL\nUSER_CHARGES\n\n")
        f.write("@<TRIPOS>ATOM\n")
        for i, a in enumerate(atoms, 1):
            f.write(f"{i:>7} {a.name:<8} {a.x:>10.4f} {a.y:>10.4f} {a.z:>10.4f} "
                    f"{a.sybyl_type:<5} {a.resid:>6} {resname:<5} {a.charge:>10.6f}\n")
        f.write("@<TRIPOS>BOND\n")
        for i, (a1_idx, a2_idx, btype) in enumerate(bonds, 1):
            f.write(f"{i:>6} {a1_idx+1:>5} {a2_idx+1:>5} {btype:>4}\n")
        f.write("@<TRIPOS>SUBSTRUCTURE\n")
        f.write(f"   1 {resname:<9} 1 TEMP              0 **** **** 0 ROOT\n")


# --- 8. PUBLIC ENTRY POINT ---
def generate_mol2_files(pdb_path, chg_path, template_mol2s, debug=False):
    """
    Generate per-residue MOL2 files from a PDB, CHG, and a list of template
    MOL2 files (typically antechamber output for each unique ligand template).
    Writes <resname>.mol2 in the current working directory.

    Returns a list of paths written.
    """
    print("=" * 70)
    print("MOL2GEN_HELPER - per-residue MOL2 generation")
    print("=" * 70)

    # 1. Templates
    print("\n1. Loading Templates:")
    library = TemplateLibrary()
    for mol2_file in template_mol2s:
        if mol2_file:
            library.add_template(mol2_file)
    if not library.templates:
        raise RuntimeError("No templates loaded.")
    print(f"\n  Total templates: {len(library.templates)}")

    # 2. CHG -> charge maps
    print("\n2. Reading CHG File:")
    chg_molecules = cluster_chg_into_molecules(chg_path, exclude_metals=True)
    charge_maps, metal_charge_map = build_charge_maps(chg_molecules, library, debug=debug)

    # 3. PDB
    print("\n4. Reading PDB:")
    pdb_residues = read_pdb(pdb_path)

    # 4. Match each residue to a template
    print("\n5. Processing Residues:")
    all_processed_atoms = []
    files_to_write = []
    for resid in sorted(pdb_residues.keys()):
        atoms = pdb_residues[resid]
        resname = atoms[0].resname

        if len(atoms) == 1 and atoms[0].element in METALS:
            atom = atoms[0]
            atom.name = atom.element
            atom.sybyl_type = "M0"
            if atom.element in metal_charge_map:
                atom.charge = metal_charge_map[atom.element]
            all_processed_atoms.append(atom)
            files_to_write.append((f"{resname}.mol2", resname, [atom], []))
            print(f"  {resname} (resid {resid}): Metal, charge = {atom.charge:.6f}")
            continue

        template, mapping = library.match_pdb_to_template(atoms, debug=debug)
        if template is None:
            print(f"  {resname} (resid {resid}): No template match - SKIPPED")
            continue

        new_atoms = []
        for t_idx in range(len(template['atoms'])):
            p_idx = mapping[t_idx]
            pdb_atom = atoms[p_idx]
            tpl_name = template['atoms'][t_idx].name
            new_atom = Atom(tpl_name, pdb_atom.x, pdb_atom.y, pdb_atom.z,
                            resname, resid, pdb_atom.element)
            new_atom.sybyl_type = template['types'].get(tpl_name, "Du")
            charge_map = charge_maps.get(template['prefix'], {})
            new_atom.charge = charge_map.get(tpl_name, 0.0)
            new_atoms.append(new_atom)
        all_processed_atoms.extend(new_atoms)
        files_to_write.append((f"{resname}.mol2", resname, new_atoms, template['bonds']))
        print(f"  {resname} (resid {resid}): Matched template '{template['prefix']}', "
              f"{len(new_atoms)} atoms")

    # 5. MUNRO types
    print("\n6. Applying MUNRO Types:")
    apply_munro_types(all_processed_atoms)

    # 6. Write
    print("\n7. Writing MOL2 Files:")
    written = []
    for fname, rname, ats, bnds in files_to_write:
        write_mol2(fname, rname, ats, bnds)
        written.append(fname)
        print(f"  -> {fname}")

    print("\n" + "=" * 70)
    print(f"Done. {len(written)} MOL2 files written.")
    print("=" * 70)
    return written
