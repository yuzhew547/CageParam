#!/usr/bin/env python3
"""
mol2gen_multi_chg.py - Multi-Template MOL2 Generator with CHG Support
Generates MOL2 files for systems with multiple ligand types using CHG charges.

Features:
- Multiple reference MOL2 files (for topology and atom types)
- CHG file for charge assignment
- Automatic topology validation
- Multi-template support
"""
import sys
import argparse
import math
import os
from collections import defaultdict
import numpy as np

# ================= CONFIGURATION =================
BOND_CUTOFF = 1.90       # Max distance to infer bonds in PDB
METAL_BOND_CUTOFF = 2.5  # Max distance for Metal-Ligand bonds
# =================================================

# Covalent radii for more accurate bond detection
COVALENT_RADII = {
    'H': 0.31, 'C': 0.76, 'N': 0.71, 'O': 0.66, 'F': 0.57,
    'P': 1.07, 'S': 1.05, 'CL': 1.02, 'BR': 1.20, 'I': 1.39,
    'B': 0.84, 'SI': 1.11, 'PD': 1.39, 'PT': 1.36, 'AU': 1.36, 'AG': 1.45
}

# Metals to exclude when clustering organic ligands
METALS = {'PD', 'PT', 'AU', 'AG', 'FE', 'CO', 'NI', 'CU', 'ZN', 'RU', 'RH', 'IR'}

def get_bond_cutoff(el1, el2, tolerance=0.45):
    """Get appropriate bond cutoff based on covalent radii"""
    r1 = COVALENT_RADII.get(el1, 0.77)
    r2 = COVALENT_RADII.get(el2, 0.77)
    return r1 + r2 + tolerance

# --- 1. UTILITIES ---
def clean(text):
    if not text: return ""
    return text.strip().upper()

def get_unique_2char_code(index, prefixes):
    chars = [str(i) for i in range(0, 10)] + \
            [chr(i) for i in range(ord('A'), ord('Z')+1)] + \
            [chr(i) for i in range(ord('a'), ord('z')+1)]
    real_idx = index 
    prefix_idx = real_idx // len(chars)
    char_idx = real_idx % len(chars)
    if prefix_idx >= len(prefixes): return "XX"
    return f"{prefixes[prefix_idx]}{chars[char_idx]}"

def dist_sq(a1, a2):
    return (a1.x - a2.x)**2 + (a1.y - a2.y)**2 + (a1.z - a2.z)**2

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
        
        # Element inference
        if element:
            self.element = element.upper()
        else:
            # Guess from name: extract letters, take first 1-2
            letters = "".join([c for c in self.name if c.isalpha()])
            if letters:
                # Common 2-letter elements
                if letters[:2] in ['BR', 'CL', 'PD', 'PT', 'AU', 'AG']:
                    self.element = letters[:2]
                else:
                    self.element = letters[0]
            else:
                self.element = "X"

# --- 3. GRAPH ENGINE (ISOMORPHISM) ---
def build_adjacency(atoms, use_covalent_radii=True):
    """Builds graph from coordinates using covalent radii for bond detection."""
    adj = defaultdict(list)
    n = len(atoms)
    
    for i in range(n):
        for j in range(i + 1, n):
            if use_covalent_radii:
                cutoff = get_bond_cutoff(atoms[i].element, atoms[j].element)
            else:
                cutoff = BOND_CUTOFF
            
            d = math.sqrt(dist_sq(atoms[i], atoms[j]))
            if d < cutoff:
                adj[i].append(j)
                adj[j].append(i)
    return adj

def get_extended_signature(idx, atoms, adj, depth=2):
    """Creates an extended signature including neighbors of neighbors."""
    el = atoms[idx].element
    neighbors = adj[idx]
    
    # First level neighbors
    n_els = tuple(sorted([atoms[n].element for n in neighbors]))
    
    if depth >= 2:
        # Second level - neighbors of neighbors
        second_level = []
        for n in neighbors:
            second_n_els = tuple(sorted([atoms[nn].element for nn in adj[n] if nn != idx]))
            second_level.append(second_n_els)
        second_level = tuple(sorted(second_level))
        return (el, len(neighbors), n_els, second_level)
    
    return (el, len(neighbors), n_els)

def solve_isomorphism(templ_atoms, templ_adj, target_atoms, target_adj, debug=False):
    """
    Maps Template Indices -> Target Indices using Graph Isomorphism (Backtracking).
    Returns: dict mapping templ_idx -> target_idx, or None if no match
    """
    if len(templ_atoms) != len(target_atoms):
        return None

    # Pre-calc signatures
    t_sigs = [get_extended_signature(i, templ_atoms, templ_adj) for i in range(len(templ_atoms))]
    tgt_sigs = [get_extended_signature(i, target_atoms, target_adj) for i in range(len(target_atoms))]
    
    mapping = {}
    used_targets = set()
    
    # Sort template atoms by signature uniqueness
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
            
            # Check Connectivity to already mapped neighbors
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

    if backtrack(0):
        return mapping
    return None

# --- 4. READERS ---
def read_template(filepath):
    """Read template mol2 file - topology, types, and bonds."""
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
            except: 
                pass
            
    # Build Adjacency from explicit bonds
    adj = defaultdict(list)
    for a1, a2, _ in bonds:
        adj[a1].append(a2)
        adj[a2].append(a1)
    
    print(f"    -> {len(atoms)} atoms, {len(bonds)} bonds")
    return atoms, adj, bonds, atom_types

def read_pdb(filepath):
    """Read PDB and group atoms by residue."""
    print(f"\n  Reading PDB: {filepath}")
    residues = defaultdict(list)
    
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                name = line[12:16]
                resname = line[17:20]
                try: 
                    resid = int(line[22:26])
                except: 
                    resid = 1
                x, y, z = line[30:38], line[38:46], line[46:54]
                
                el = line[76:78].strip() if len(line) > 77 else None
                
                if clean(resname) not in ["WAT", "HOH", "TIP3", "TIP"]: 
                    residues[resid].append(Atom(name, x, y, z, resname, resid, el))
    
    print(f"    -> {len(residues)} residues")
    return residues

def cluster_chg_into_molecules(filepath, exclude_metals=True):
    """
    Read CHG file and cluster atoms into separate molecules.
    """
    print(f"\n  Reading CHG file: {filepath}")
    
    chg_atoms = []
    with open(filepath, 'r') as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 5:
                try:
                    element = parts[0].upper()
                    x = float(parts[1])
                    y = float(parts[2])
                    z = float(parts[3])
                    charge = float(parts[4])
                    chg_atoms.append({
                        'element': element,
                        'x': x, 'y': y, 'z': z,
                        'charge': charge
                    })
                except (ValueError, IndexError):
                    continue
    
    print(f"    -> Loaded {len(chg_atoms)} atoms from CHG")
    
    # Separate metals from organic atoms
    metal_atoms = []
    organic_atoms = []
    
    for i, atom in enumerate(chg_atoms):
        if exclude_metals and atom['element'] in METALS:
            metal_atoms.append(atom)
        else:
            organic_atoms.append(atom)
    
    if exclude_metals:
        print(f"    -> Separated {len(metal_atoms)} metal atoms, {len(organic_atoms)} organic atoms")
    
    # Cluster organic atoms into molecules
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
                d = math.sqrt(dx*dx + dy*dy + dz*dz)
                if d < cutoff:
                    dfs(j, current_mol)
    
    for i in range(n):
        if not visited[i]:
            current_mol = []
            dfs(i, current_mol)
            molecules.append(current_mol)
    
    # Add metals as single-atom molecules
    for metal in metal_atoms:
        molecules.append([metal])
    
    print(f"    -> Clustered into {len(molecules)} molecules")
    return molecules

# --- 5. TEMPLATE LIBRARY ---
class TemplateLibrary:
    """Manages multiple ligand templates."""
    
    def __init__(self):
        self.templates = []  # List of template dicts
        
    def add_template(self, mol2_file, resname_prefix=None):
        """Add a template from MOL2 file."""
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
        
        template = {
            'prefix': resname_prefix,
            'atoms': atoms,
            'adj': adj,
            'bonds': bonds,
            'types': types,
            'mol2_file': mol2_file
        }
        
        self.templates.append(template)
        print(f"    Added template '{resname_prefix}' ({len(atoms)} atoms)")
        
    def match_pdb_to_template(self, pdb_atoms, debug=False):
        """Match PDB atoms to a template. Returns (template, mapping) or (None, None)."""
        pdb_adj = build_adjacency(pdb_atoms, use_covalent_radii=True)
        
        for template in self.templates:
            if len(template['atoms']) != len(pdb_atoms):
                continue
            
            mapping = solve_isomorphism(
                template['atoms'], template['adj'],
                pdb_atoms, pdb_adj,
                debug=debug
            )
            
            if mapping is not None:
                return template, mapping
        
        return None, None
    
    def match_chg_to_template(self, chg_mol, debug=False):
        """Match CHG molecule to a template. Returns (template, mapping) or (None, None)."""
        # Create Atom objects from CHG data
        chg_atoms = [Atom(a['element'], a['x'], a['y'], a['z'], "CHG", 1, a['element']) 
                     for a in chg_mol]
        chg_adj = build_adjacency(chg_atoms, use_covalent_radii=True)
        
        for template in self.templates:
            if len(template['atoms']) != len(chg_atoms):
                continue
            
            mapping = solve_isomorphism(
                template['atoms'], template['adj'],
                chg_atoms, chg_adj,
                debug=debug
            )
            
            if mapping is not None:
                return template, mapping
        
        return None, None

# --- 6. CHARGE ASSIGNMENT ---
def build_charge_maps(chg_molecules, template_library, debug=False):
    """
    Build charge maps for each template by matching CHG molecules.
    Returns dict: template_prefix -> {atom_name -> avg_charge}
    """
    print(f"\n3. Validating CHG molecules against templates:")
    
    charges_by_template = defaultdict(lambda: defaultdict(list))
    metal_charges = defaultdict(list)
    
    for i, chg_mol in enumerate(chg_molecules):
        # Handle metals
        if len(chg_mol) == 1:
            element = chg_mol[0]['element']
            if element in METALS:
                metal_charges[element].append(chg_mol[0]['charge'])
                print(f"  Molecule {i+1}: Metal ({element})")
            continue
        
        # Try to match to templates
        template, mapping = template_library.match_chg_to_template(chg_mol, debug=debug)
        
        if template is not None:
            print(f"  Molecule {i+1}: Matched template '{template['prefix']}' ({len(chg_mol)} atoms)")
            
            # Store charges with template atom names
            for t_idx in range(len(template['atoms'])):
                c_idx = mapping[t_idx]
                tpl_name = template['atoms'][t_idx].name
                charge = chg_mol[c_idx]['charge']
                charges_by_template[template['prefix']][tpl_name].append(charge)
        else:
            print(f"  Molecule {i+1}: No template match ({len(chg_mol)} atoms)")
    
    # Compute averages
    charge_maps = {}
    
    for template_prefix, atom_charges in charges_by_template.items():
        avg_charges = {}
        print(f"\n  Average charges for template '{template_prefix}':")
        print(f"    {'Atom':<6} {'Avg Charge':>10} {'Std Dev':>10} {'Count':>6}")
        print(f"    {'-'*6} {'-'*10} {'-'*10} {'-'*6}")
        
        for atom_name in sorted(atom_charges.keys()):
            charges = atom_charges[atom_name]
            avg = np.mean(charges)
            std = np.std(charges) if len(charges) > 1 else 0.0
            avg_charges[atom_name] = avg
            print(f"    {atom_name:<6} {avg:>10.4f} {std:>10.4f} {len(charges):>6}")
        
        charge_maps[template_prefix] = avg_charges
    
    # Metal charges
    metal_charge_map = {}
    if metal_charges:
        print(f"\n  Metal charges:")
        for element, charges in metal_charges.items():
            avg = np.mean(charges)
            metal_charge_map[element] = avg
            print(f"    {element}: {avg:.6f} (from {len(charges)} atoms)")
    
    return charge_maps, metal_charge_map

# --- 7. PROCESSING ---
def apply_munro_types(all_atoms):
    """Apply special types to metals and coordinating nitrogens."""
    print("\n  Applying MUNRO types (M, Y codes)...")
    
    # Find all metals
    metals = [a for a in all_atoms if a.element in METALS]
    metals.sort(key=lambda a: (a.resid, a.name))
    
    for i, m in enumerate(metals):
        m.sybyl_type = get_unique_2char_code(i, "MNO")
    
    print(f"    -> {len(metals)} metal atoms")

    # Find coordinating nitrogens
    nitrogens = [a for a in all_atoms if a.element == "N"]
    nitrogens.sort(key=lambda a: (a.resid, a.name))
    
    cutoff_sq = METAL_BOND_CUTOFF**2
    y_cnt = 0
    
    for n in nitrogens:
        connected = False
        for m in metals:
            if dist_sq(n, m) < cutoff_sq:
                connected = True
                break
        if connected:
            n.sybyl_type = get_unique_2char_code(y_cnt, "YZWVU")
            y_cnt += 1
    
    print(f"    -> {y_cnt} coordinating nitrogens")

def write_mol2(filename, resname, atoms, bonds):
    """Write MOL2 file."""
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

# --- MAIN ---
def main():
    parser = argparse.ArgumentParser(
        description="MOL2GEN Multi-Template with CHG Support",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Single template (backward compatible):
  python mol2gen_multi_chg.py bone.pdb cage.chg template.mol2
  
  # Multiple templates:
  python mol2gen_multi_chg.py bone.pdb cage.chg \\
      template1.mol2 template2.mol2 template3.mol2 template4.mol2
  
  # With explicit prefixes:
  python mol2gen_multi_chg.py bone.pdb cage.chg \\
      --templates LA:temp1.mol2 LB:temp2.mol2 LC:temp3.mol2
        """
    )
    
    parser.add_argument("pdb", help="PDB file")
    parser.add_argument("chg", help="CHG file with charges")
    
    # Template specification
    template_group = parser.add_mutually_exclusive_group(required=True)
    template_group.add_argument("templates_pos", nargs='*', metavar="TEMPLATE", 
                               help="Template MOL2 files (positional)")
    template_group.add_argument("--templates", nargs='+', 
                               help="Templates as PREFIX:FILE (e.g., LA:temp1.mol2)")
    
    parser.add_argument("--debug", action="store_true", help="Enable debug output")
    
    args = parser.parse_args()
    
    print("="*70)
    print("MOL2GEN Multi-Template with CHG Support")
    print("="*70)
    
    # 1. Load Templates
    print("\n1. Loading Templates:")
    library = TemplateLibrary()
    
    # Handle positional templates (backward compatible)
    if args.templates_pos:
        for mol2_file in args.templates_pos:
            if mol2_file:  # Skip empty strings
                library.add_template(mol2_file)
    
    # Handle explicit prefix templates
    if args.templates:
        for spec in args.templates:
            if ':' in spec:
                prefix, filepath = spec.split(':', 1)
                library.add_template(filepath, prefix)
            else:
                library.add_template(spec)
    
    if not library.templates:
        print("ERROR: No templates loaded!")
        sys.exit(1)
    
    print(f"\n  Total templates: {len(library.templates)}")
    
    # 2. Read CHG and build charge maps
    print("\n2. Reading CHG File:")
    chg_molecules = cluster_chg_into_molecules(args.chg, exclude_metals=True)
    charge_maps, metal_charge_map = build_charge_maps(chg_molecules, library, debug=args.debug)
    
    # 3. Read PDB
    print("\n4. Reading PDB:")
    pdb_residues = read_pdb(args.pdb)
    
    # 4. Process residues
    print("\n5. Processing Residues:")
    all_processed_atoms = []
    files_to_write = []
    
    for resid in sorted(pdb_residues.keys()):
        atoms = pdb_residues[resid]
        resname = atoms[0].resname
        
        # Is this a metal?
        if len(atoms) == 1 and atoms[0].element in METALS:
            atom = atoms[0]
            atom.name = atom.element
            atom.sybyl_type = "M0"  # Temp, will be fixed later
            
            # Assign charge from metal charge map
            if atom.element in metal_charge_map:
                atom.charge = metal_charge_map[atom.element]
            
            all_processed_atoms.append(atom)
            files_to_write.append((f"{resname}.mol2", resname, [atom], []))
            print(f"  {resname} (resid {resid}): Metal, charge = {atom.charge:.6f}")
            continue
        
        # Ligand residue - match to template
        template, mapping = library.match_pdb_to_template(atoms, debug=args.debug)
        
        if template is None:
            print(f"  {resname} (resid {resid}): No template match - SKIPPED")
            continue
        
        # Create reordered atoms matching template
        new_atoms = []
        for t_idx in range(len(template['atoms'])):
            p_idx = mapping[t_idx]
            pdb_atom = atoms[p_idx]
            
            tpl_name = template['atoms'][t_idx].name
            
            new_atom = Atom(tpl_name, pdb_atom.x, pdb_atom.y, pdb_atom.z, 
                           resname, resid, pdb_atom.element)
            
            new_atom.sybyl_type = template['types'].get(tpl_name, "Du")
            
            # Assign charge from charge map
            charge_map = charge_maps.get(template['prefix'], {})
            new_atom.charge = charge_map.get(tpl_name, 0.0)
            
            new_atoms.append(new_atom)
        
        all_processed_atoms.extend(new_atoms)
        files_to_write.append((f"{resname}.mol2", resname, new_atoms, template['bonds']))
        print(f"  {resname} (resid {resid}): Matched template '{template['prefix']}', {len(new_atoms)} atoms")
    
    # 5. Apply MUNRO types
    print("\n6. Applying MUNRO Types:")
    apply_munro_types(all_processed_atoms)
    
    # 6. Write files
    print("\n7. Writing MOL2 Files:")
    for fname, rname, ats, bnds in files_to_write:
        write_mol2(fname, rname, ats, bnds)
        print(f"  -> {fname}")
    
    print("\n" + "="*70)
    print("✓ Done!")
    print("="*70)

if __name__ == "__main__":
    main()
