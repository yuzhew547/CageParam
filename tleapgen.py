#!/usr/bin/env python3
"""
tleap_generator.py - Generate tleap.in files for Pd/Pt cage simulations

This module generates AMBER tleap input files by:
1. Parsing PDB files to identify residues (ligands, metals, anions)
2. Calculating metal-nitrogen connectivity to assign unique atom types
3. Writing a complete tleap.in file with all necessary sections

Usage:
    python tleap_generator.py -p input.pdb -o tleap.in [options]
"""

import argparse
import math
from collections import defaultdict, OrderedDict
import os


# --- GEOMETRY HELPERS ---
def dist(c1, c2):
    """Calculate Euclidean distance between two 3D coordinates."""
    return math.sqrt(sum((c1[i] - c2[i])**2 for i in range(3)))


def get_unique_2char_code(index, prefixes):
    """
    Generate unique 2-character codes for atom types.
    
    Args:
        index: 1-based index
        prefixes: String of prefix characters (e.g., "MNO" for metals, "YZWVU" for nitrogens)
    
    Returns:
        2-character code string
    """
    chars = [str(i) for i in range(0, 10)] + \
            [chr(i) for i in range(ord('A'), ord('Z')+1)] + \
            [chr(i) for i in range(ord('a'), ord('z')+1)]
    if index < 1:
        return "XX"
    real_idx = index - 1
    prefix_idx = real_idx // len(chars)
    char_idx = real_idx % len(chars)
    if prefix_idx >= len(prefixes):
        return "XX"
    return f"{prefixes[prefix_idx]}{chars[char_idx]}"


def parse_pdb(pdb_file):
    """
    Parse PDB file to extract atom information and residue data.
    
    Returns:
        Dictionary containing all parsed data
    """
    pd_residues, ml_residues, anion_residues = [], [], []
    pd_atoms, n_atoms, all_ligand_atoms = [], [], []
    metal_types = {}
    residue_names = {}
    residue_order = []  # Track order of all residues as they appear
    seen_res_keys = set()
    
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                chain = line[21].strip() or "-"
                try:
                    res_num = int(line[22:26].strip())
                except ValueError:
                    continue
                res_key = (chain, res_num)
                name = line[12:16].strip()
                res_name = line[17:20].strip().upper()
                
                # Handle 4-character residue names (columns 17-20 may overflow)
                if len(line) > 20:
                    # Try to get extended residue name
                    extended_res = line[17:21].strip().upper()
                    if extended_res and not extended_res[0].isdigit():
                        res_name = extended_res
                
                try:
                    coords = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
                except ValueError:
                    continue
                    
                atom_data = {'res_key': res_key, 'name': name, 'coords': coords, 'res_name': res_name}
                
                # Track residue order
                if res_key not in seen_res_keys:
                    seen_res_keys.add(res_key)
                    residue_order.append(res_key)
                    residue_names[res_key] = res_name
                
                # Metal residues (P1, P2, etc. but not PRO)
                if res_name.startswith("P") and not res_name.startswith("PRO"):
                    if res_key not in pd_residues:
                        pd_residues.append(res_key)
                    pd_atoms.append(atom_data)
                    metal_types[res_key] = name.upper()
                
                # Ligand residues (L1, L2, etc. but not LEU, LYS)
                elif res_name.startswith("L") and not res_name.startswith(("LEU", "LYS")):
                    if res_key not in ml_residues:
                        ml_residues.append(res_key)
                    all_ligand_atoms.append(atom_data)
                    if "N" in name.upper():
                        n_atoms.append(atom_data)
                
                # Anion residues (BF*, etc.)
                elif res_name.startswith("BF"):
                    if res_key not in anion_residues:
                        anion_residues.append(res_key)
    
    return {
        'pd_residues': pd_residues,
        'ml_residues': ml_residues,
        'anion_residues': anion_residues,
        'pd_atoms': pd_atoms,
        'n_atoms': n_atoms,
        'all_ligand_atoms': all_ligand_atoms,
        'metal_types': metal_types,
        'residue_names': residue_names,
        'residue_order': residue_order
    }


def calculate_connectivity(pdb_data, metal_n_cutoff=3.0):
    """
    Calculate metal-nitrogen connectivity to assign unique atom type codes.
    
    Uses the same algorithm as munro_frcmod.py to ensure consistency.
    """
    pd_residues = pdb_data['pd_residues']
    ml_residues = pdb_data['ml_residues']
    pd_atoms = pdb_data['pd_atoms']
    n_atoms = pdb_data['n_atoms']
    
    pd_map = {k: i + 1 for i, k in enumerate(pd_residues)}
    ml_map = {k: i + 1 for i, k in enumerate(ml_residues)}
    
    # Group coordination candidates by ligand index
    coord_candidates = defaultdict(list)
    for pd in pd_atoms:
        for n in n_atoms:
            if dist(pd['coords'], n['coords']) < metal_n_cutoff:
                l_idx = ml_map[n['res_key']]
                if n not in coord_candidates[l_idx]:
                    coord_candidates[l_idx].append(n)
    
    # Assign Y codes to coordinating nitrogens (same order as munro script)
    y_atoms_map = {}
    y_counter = 1
    for l_idx in sorted(coord_candidates.keys()):
        atoms = sorted(coord_candidates[l_idx], key=lambda x: x['name'])
        for atom in atoms:
            uid = (atom['res_key'], atom['name'])
            if uid not in y_atoms_map:
                y_code = get_unique_2char_code(y_counter, "YZWVU")
                y_atoms_map[uid] = {'code': y_code, 'data': atom}
                y_counter += 1
    
    # Build metal connectivity with bond information
    pd_connectivity = defaultdict(list)
    metal_codes = {}
    bond_info = []  # List of (ligand_mol_idx, n_atom_name, metal_mol_idx)
    
    for pd in pd_atoms:
        m_idx = pd_map[pd['res_key']]
        m_code = get_unique_2char_code(m_idx, "MNO")
        metal_codes[m_idx] = m_code
        
        for n in n_atoms:
            uid = (n['res_key'], n['name'])
            if uid in y_atoms_map and dist(pd['coords'], n['coords']) < metal_n_cutoff:
                info = y_atoms_map[uid]
                pd_connectivity[m_idx].append({
                    'y_code': info['code'],
                    'coords': n['coords'],
                    'data': n
                })
                
                # Store bond info for tleap
                l_idx = ml_map[n['res_key']]
                bond_info.append({
                    'ligand_idx': l_idx,
                    'n_name': n['name'],
                    'metal_idx': m_idx,
                    'metal_res_key': pd['res_key']
                })
    
    return {
        'pd_map': pd_map,
        'ml_map': ml_map,
        'y_atoms_map': y_atoms_map,
        'pd_connectivity': pd_connectivity,
        'metal_codes': metal_codes,
        'bond_info': bond_info
    }


def generate_atom_types(pdb_data, connectivity):
    """
    Generate addAtomTypes section content.
    """
    atom_types = []
    metal_types = pdb_data['metal_types']
    pd_residues = pdb_data['pd_residues']
    y_atoms_map = connectivity['y_atoms_map']
    metal_codes = connectivity['metal_codes']
    pd_map = connectivity['pd_map']
    
    # Add metal types (M1, M2, etc.)
    for res_key in pd_residues:
        m_idx = pd_map[res_key]
        m_code = metal_codes[m_idx]
        el_type = metal_types.get(res_key, "PD")
        element = "Pt" if "PT" in el_type.upper() else "Pd"
        atom_types.append(f'        {{ "{m_code}"  "{element}" "sp3" }}')
    
    # Add nitrogen types (Y1, Y2, etc.) - sorted by code
    for uid, info in sorted(y_atoms_map.items(), key=lambda x: x[1]['code']):
        y_code = info['code']
        atom_types.append(f'        {{ "{y_code}"  "N" "sp3" }}')
    
    return atom_types


def generate_bond_commands(pdb_data, connectivity):
    """
    Generate bond commands for metal-nitrogen connections.
    
    Format: bond mol.{ligand_residue_num}.{N_atom_name} mol.{metal_residue_num}.PD
    
    The residue numbers are based on the order in the PDB file (1-indexed).
    Ligands come first, then metals, then anions.
    """
    bonds = []
    ml_residues = pdb_data['ml_residues']
    pd_residues = pdb_data['pd_residues']
    bond_info = connectivity['bond_info']
    
    # Calculate the mol index for metals (ligands + metal index)
    num_ligands = len(ml_residues)
    
    # Group bonds by metal for organized output
    bonds_by_metal = defaultdict(list)
    for b in bond_info:
        bonds_by_metal[b['metal_idx']].append(b)
    
    # Generate bond commands sorted by metal index
    for m_idx in sorted(bonds_by_metal.keys()):
        metal_bonds = bonds_by_metal[m_idx]
        # Sort by ligand index then by N atom name
        metal_bonds.sort(key=lambda x: (x['ligand_idx'], x['n_name']))
        
        for b in metal_bonds:
            ligand_mol_idx = b['ligand_idx']
            metal_mol_idx = num_ligands + b['metal_idx']
            n_name = b['n_name']
            bonds.append(f"bond mol.{ligand_mol_idx}.{n_name} mol.{metal_mol_idx}.PD")
    
    return bonds


def get_unique_residue_names(pdb_data):
    """
    Get unique residue names organized by type.
    """
    residue_names = pdb_data['residue_names']
    
    ligands = set()
    metals = set()
    anions = set()
    
    for res_key, res_name in residue_names.items():
        if res_name.startswith("L") and not res_name.startswith(("LEU", "LYS")):
            ligands.add(res_name)
        elif res_name.startswith("P") and not res_name.startswith("PRO"):
            metals.add(res_name)
        elif res_name.startswith("BF"):
            anions.add(res_name)
    
    # Sort ligands and metals numerically
    def sort_key(name):
        # Extract number from name like L1, L2, L10, P1, P2, etc.
        prefix = ''.join(c for c in name if c.isalpha())
        num_str = ''.join(c for c in name if c.isdigit())
        num = int(num_str) if num_str else 0
        return (prefix, num)
    
    return {
        'ligands': sorted(ligands, key=sort_key),
        'metals': sorted(metals, key=sort_key),
        'anions': sorted(anions)
    }


def write_tleap_file(output_file, pdb_file, pdb_data, connectivity,
                     frcmod_file="munro.frcmod",
                     solvent_lib=None,
                     solvent_box_size=15.0,
                     output_prefix="ori"):
    """
    Write the complete tleap.in file.
    """
    residues = get_unique_residue_names(pdb_data)
    atom_types = generate_atom_types(pdb_data, connectivity)
    bonds = generate_bond_commands(pdb_data, connectivity)
    
    # Get PDB filename without path
    pdb_basename = os.path.basename(pdb_file)
    
    with open(output_file, 'w') as f:
        # Source force field files
        f.write("source leaprc.protein.ff19SB\n")
        f.write("source leaprc.gaff2\n")
        f.write("source leaprc.water.opc\n")
        f.write("\n")
        
        # Add atom types
        f.write("addAtomTypes {\n")
        for line in atom_types:
            f.write(line + "\n")
        f.write("}\n")
        f.write("\n")
        
        # Load ligand mol2 files
        for lig in residues['ligands']:
            f.write(f"{lig} = loadmol2 {lig}.mol2\n")
        
        # Load metal mol2 files
        for metal in residues['metals']:
            f.write(f"{metal} = loadmol2 {metal}.mol2\n")
        f.write("\n")
        
        # Load anion mol2 files
        for anion in residues['anions']:
            f.write(f"{anion} = loadmol2 {anion}.mol2\n")
        f.write("\n")
        
        # Load frcmod files for anions
        for anion in residues['anions']:
            f.write(f"loadamberparams {anion}.frcmod\n")
        
        # Load ion parameters and main frcmod
        f.write("loadamberparams frcmod.ionslm_126_opc\n")
        f.write(f"loadamberparams {frcmod_file}\n")
        f.write("\n")
        
        # Load PDB
        f.write(f"mol = loadpdb {pdb_basename}\n")
        f.write("\n")
        
        # Write bond commands
        for bond in bonds:
            f.write(bond + "\n")
        f.write("\n")
        
        # Save dry structure
        f.write(f"savepdb mol {output_prefix}_dry.pdb\n")
        f.write(f"saveamberparm mol {output_prefix}_dry.prmtop {output_prefix}_dry.inpcrd\n")
        
        # Solvation (if solvent library specified)
        if solvent_lib:
            f.write(f"loadoff {solvent_lib}\n")
            # Extract solvent name from library filename
            solvent_name = os.path.splitext(solvent_lib)[0]
            f.write(f"solvateBox mol {solvent_name} {solvent_box_size}\n")
            f.write(f"savepdb mol {output_prefix}_solv.pdb\n")
            f.write(f"saveamberparm mol {output_prefix}_solv.prmtop {output_prefix}_solv.inpcrd\n")
        
        f.write("quit\n")


def main():
    parser = argparse.ArgumentParser(
        description="Generate tleap.in files for Pd/Pt cage simulations",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Basic usage
    python tleap_generator.py -p cage.pdb -o tleap.in
    
    # With custom solvent
    python tleap_generator.py -p cage.pdb -o tleap.in --solvent dmso_acrylate.lib
    
    # With custom output prefix and box size
    python tleap_generator.py -p cage.pdb -o tleap.in --prefix my_system --box-size 20.0
"""
    )
    
    parser.add_argument("-p", "--pdb", required=True,
                        help="Input PDB file")
    parser.add_argument("-o", "--output", default="tleap.in",
                        help="Output tleap.in file (default: tleap.in)")
    parser.add_argument("-f", "--frcmod", default="munro.frcmod",
                        help="Main frcmod file name (default: munro.frcmod)")
    parser.add_argument("--solvent", default=None,
                        help="Solvent library file (e.g., dmso_acrylate.lib)")
    parser.add_argument("--box-size", type=float, default=15.0,
                        help="Solvent box size in Angstroms (default: 15.0)")
    parser.add_argument("--prefix", default="ori",
                        help="Output file prefix (default: ori)")
    parser.add_argument("--cutoff", type=float, default=3.0,
                        help="Metal-nitrogen distance cutoff in Angstroms (default: 3.0)")
    
    args = parser.parse_args()
    
    # Validate input file
    if not os.path.exists(args.pdb):
        print(f"Error: PDB file '{args.pdb}' not found.")
        return 1
    
    print(f"Parsing PDB file: {args.pdb}")
    pdb_data = parse_pdb(args.pdb)
    
    print(f"Found {len(pdb_data['pd_residues'])} metal residues")
    print(f"Found {len(pdb_data['ml_residues'])} ligand residues")
    print(f"Found {len(pdb_data['anion_residues'])} anion residues")
    print(f"Found {len(pdb_data['n_atoms'])} nitrogen atoms")
    
    print("Calculating metal-nitrogen connectivity...")
    connectivity = calculate_connectivity(pdb_data, metal_n_cutoff=args.cutoff)
    
    print(f"Identified {len(connectivity['y_atoms_map'])} coordinating nitrogens")
    print(f"Generated {len(connectivity['bond_info'])} metal-nitrogen bonds")
    
    print(f"Writing tleap.in file: {args.output}")
    write_tleap_file(
        output_file=args.output,
        pdb_file=args.pdb,
        pdb_data=pdb_data,
        connectivity=connectivity,
        frcmod_file=args.frcmod,
        solvent_lib=args.solvent,
        solvent_box_size=args.box_size,
        output_prefix=args.prefix
    )
    
    print(f"(U^(I)^U) tleap.in file saved to {args.output}")
    return 0


if __name__ == "__main__":
    exit(main())
