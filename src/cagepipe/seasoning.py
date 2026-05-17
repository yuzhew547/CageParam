#!/usr/bin/env python3
import math
import sys
import shutil
import string
import os
import random
import re

# ========= DEFAULT SETTINGS =========
DEFAULT_CLASH_DIST = 4.0   
DEFAULT_BUFFER     = 6.0   
MAX_ATTEMPTS       = 10000 
# ====================================

def get_user_input(prompt, default=None, is_file=False, is_int=False):
    while True:
        val_str = f" [{default}]" if default is not None else ""
        user_val = input(f"{prompt}{val_str}: ").strip()
        if not user_val:
            if default is not None: return default
            else:
                print("Error: This value is required.")
                continue
        if is_file and not os.path.exists(user_val):
            print(f"Error: File '{user_val}' not found.")
            continue
        if is_int:
            try: return int(user_val)
            except ValueError:
                print("Error: Please enter a valid integer.")
                continue
        return user_val

def parse_pdb_line(line):
    # Handle files with missing first character
    if line.startswith("TOM"): line = "A" + line
    return {
        'record': line[:6].strip(),
        'serial': int(line[6:11]),
        'atom_name': line[12:16], # Keep spaces for alignment logic later
        'alt_loc': line[16],      # <--- This was the culprit (Col 17)
        'res_name': line[17:20].strip(),
        'chain': line[21],
        'res_seq': int(line[22:26]),
        'x': float(line[30:38]),
        'y': float(line[38:46]),
        'z': float(line[46:54]),
        'line_rest': line[54:].rstrip() 
    }

def read_pdb_atoms(filename):
    atoms, headers = [], []
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM", "TOM")):
                try: atoms.append(parse_pdb_line(line))
                except: continue
            elif line.startswith(("TER", "END", "CONECT", "MASTER")): continue
            else: headers.append(line)
    return atoms, headers

def get_mol2_resname(mol2_file):
    """Extract residue name from MOL2 file (line 2 after @<TRIPOS>MOLECULE)"""
    with open(mol2_file, 'r') as f:
        for line in f:
            if line.startswith("@<TRIPOS>MOLECULE"):
                # Next line is the molecule/residue name
                resname = f.readline().strip()
                return resname
    return None

def get_geometry_stats(atoms):
    if not atoms: return (0,0,0), 0.0
    cx = sum(a['x'] for a in atoms) / len(atoms)
    cy = sum(a['y'] for a in atoms) / len(atoms)
    cz = sum(a['z'] for a in atoms) / len(atoms)
    max_r = max(math.sqrt((a['x']-cx)**2 + (a['y']-cy)**2 + (a['z']-cz)**2) for a in atoms)
    return (cx, cy, cz), max_r

def check_clash(new_coords, fixed_atoms, cutoff):
    cutoff_sq = cutoff * cutoff
    for (nx, ny, nz) in new_coords:
        for atom in fixed_atoms:
            dx, dy, dz = nx - atom['x'], ny - atom['y'], nz - atom['z']
            if (dx*dx + dy*dy + dz*dz) < cutoff_sq: return True
    return False

def write_atom_line(atom, serial, res_seq, res_name, coords):
    """
    Formats an atom dictionary into a STRICT PDB line (80 chars).
    Handles Atom Name alignment and Chain ID spacing.
    """
    x, y, z = coords
    
    # 1. Fix Atom Name Alignment (Cols 13-16)
    raw_name = atom['atom_name'].strip()
    if len(raw_name) < 4:
        atom_name_fmt = f" {raw_name:<3}" # " C14"
    else:
        atom_name_fmt = f"{raw_name:<4}"  # "1HH3"

    # 2. Fix Chain ID (Col 22)
    chain_id = atom['chain'].strip()
    if not chain_id: chain_id = " "
    else: chain_id = chain_id[0]

    return (
        f"{atom['record']:<6}"       # Cols 1-6
        f"{serial:>5d}"              # Cols 7-11
        " "                          # Col  12
        f"{atom_name_fmt[:4]}"       # Cols 13-16
        f"{atom['alt_loc']:1}"       # Col  17
        f"{res_name:>3}"             # Cols 18-20
        " "                          # Col  21
        f"{chain_id}"                # Col  22
        f"{res_seq:>4d}"             # Cols 23-26
        "    "                       # Cols 27-30
        f"{x:8.3f}{y:8.3f}{z:8.3f}"  # Cols 31-54
        f"{atom['line_rest']}\n"     # Cols 55+
    )

def write_mol2(template_file, output_file, old_resname, new_resname):
    """
    Copy MOL2 file and replace residue name.
    Uses word-boundary matching to avoid partial replacements.
    """
    if template_file == output_file: return
    with open(template_file, 'r') as f:
        content = f.read()
    # Use word boundary replacement to avoid partial matches
    content = re.sub(r'\b' + re.escape(old_resname) + r'\b', new_resname, content)
    with open(output_file, 'w') as f:
        f.write(content)

def name_generator(used_names):
    """
    Generate unique 3-letter residue names for anions.
    Sequence: BFA, BFB, BFC, ... BFZ, BGA, BGB, ... BZZ
    """
    alphabet = string.ascii_uppercase
    # 3-letter names starting with BF (BFA, BFB, ... BFZ)
    for char in alphabet:
        name = f"BF{char}"
        if name not in used_names: yield name
    # 3-letter names B + [G-Z] + [A-Z] (BGA, BGB, ... BZZ)
    for char2 in alphabet:
        if char2 <= 'F': continue  # Skip BA*, BB*, BC*, BD*, BE*, BF* (already done BF*)
        for char3 in alphabet:
            name = f"B{char2}{char3}"
            if name not in used_names: yield name

def main():
    print("--- Robust Anion Placer (Aligned) ---")
    
    sys_pdb_file = get_user_input("Input System PDB (Cage)", "bone.pdb", is_file=True)
    out_pdb_file = get_user_input("Output PDB Name", "tastybone.pdb")
    n_anions = get_user_input("How many anions to add?", 24, is_int=True)
    
    templ_pdb_file = get_user_input("Template Anion PDB", "BFA.pdb", is_file=True)
    templ_frc_file = get_user_input("Template Anion FRCMOD", "BFA.frcmod", is_file=True)
    templ_mol2_file = get_user_input("Template Anion MOL2", "BFA.mol2", is_file=True)

    sys_atoms, sys_headers = read_pdb_atoms(sys_pdb_file)
    templ_atoms, _ = read_pdb_atoms(templ_pdb_file)
    
    if not templ_atoms or not sys_atoms:
        print("Error: Missing atoms."); sys.exit(1)

    # Get residue name from MOL2 file (not PDB!) for proper replacement
    template_resname_mol2 = get_mol2_resname(templ_mol2_file)
    if not template_resname_mol2:
        print(f"Error: Could not read residue name from {templ_mol2_file}")
        sys.exit(1)
    print(f"   Template MOL2 residue name: {template_resname_mol2}")
    
    cage_center, cage_radius = get_geometry_stats(sys_atoms)
    templ_center, _ = get_geometry_stats(templ_atoms)
    
    existing_resnames = set(a['res_name'] for a in sys_atoms)
    name_gen = name_generator(existing_resnames)
    
    min_dist = cage_radius + DEFAULT_BUFFER
    max_dist = min_dist + 15.0
    
    current_atoms = list(sys_atoms)
    added_lines = []
    last_serial = max((a['serial'] for a in sys_atoms), default=0)
    last_res_seq = max((a['res_seq'] for a in sys_atoms), default=0)
    
    count = 0
    for _ in range(n_anions):
        try: new_resname = next(name_gen)
        except StopIteration: break
            
        last_res_seq += 1
        placed = False
        
        for attempt in range(MAX_ATTEMPTS):
            theta = random.uniform(0, 2 * math.pi)
            phi = math.acos(random.uniform(-1, 1))
            r = random.uniform(min_dist, max_dist)
            
            rx = cage_center[0] + r * math.sin(phi) * math.cos(theta)
            ry = cage_center[1] + r * math.sin(phi) * math.sin(theta)
            rz = cage_center[2] + r * math.cos(phi)
            
            dx, dy, dz = (rx - templ_center[0], ry - templ_center[1], rz - templ_center[2])
            cand_coords = [(ta['x'] + dx, ta['y'] + dy, ta['z'] + dz) for ta in templ_atoms]
            
            if not check_clash(cand_coords, current_atoms, DEFAULT_CLASH_DIST):
                placed = True
                ion_lines = []
                for i, ta in enumerate(templ_atoms):
                    last_serial += 1
                    cx, cy, cz = cand_coords[i]
                    
                    new_atom = ta.copy()
                    new_atom['alt_loc'] = ' '  # <--- FIXED: Clear AltLoc from template
                    new_atom.update({'x': cx, 'y': cy, 'z': cz})
                    
                    current_atoms.append(new_atom)
                    ion_lines.append(write_atom_line(new_atom, last_serial, last_res_seq, new_resname, (cx, cy, cz)))
                
                added_lines.extend(ion_lines)
                added_lines.append("TER\n")
                
                # File Generation
                new_pdb = f"{new_resname}.pdb"
                new_mol2 = f"{new_resname}.mol2"
                new_frcmod = f"{new_resname}.frcmod"

                if new_pdb != templ_pdb_file:
                    with open(new_pdb, "w") as f:
                        f.writelines(ion_lines); f.write("TER\nEND\n")
                if new_mol2 != templ_mol2_file:
                    # Use residue name from MOL2 file for replacement
                    write_mol2(templ_mol2_file, new_mol2, template_resname_mol2, new_resname)
                if new_frcmod != templ_frc_file:
                    shutil.copy(templ_frc_file, new_frcmod)

                count += 1
                if count % 5 == 0 or count == n_anions:
                    print(f"   ... Placed {count}/{n_anions} ({new_resname})")
                break
        
        if not placed: break

    with open(out_pdb_file, 'w') as f:
        for h in sys_headers: f.write(h)
        for a in sys_atoms:
            f.write(write_atom_line(a, a['serial'], a['res_seq'], a['res_name'], (a['x'], a['y'], a['z'])))
        f.write("TER\n")
        f.writelines(added_lines)
        f.write("END\n")
    print(f"Done! Saved to: {out_pdb_file} o==[::::]==o !!!(o(I)oU) ")

if __name__ == "__main__":
    main()
