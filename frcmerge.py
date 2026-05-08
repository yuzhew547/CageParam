import os

def parse_atom_key(line, section):
    """
    Extracts the 'definition key' from a parameter line.
    Examples:
      BOND "c3-hc 300.0 ..." -> ('c3', 'hc')
      ANGL "c3-c3-hc ..."    -> ('c3', 'c3', 'hc')
    """
    parts = line.strip().split()
    if not parts: return None
    
    # MASS section keys are just the atom type (first column)
    if section == 'MASS':
        return parts[0]
    
    # Parameter sections use definitions like "t1-t2-t3" or "t1-t2"
    # FRCMOD standard is typically dash-separated in the first column.
    label = parts[0]
    if '-' not in label: 
        return None # Skip header lines or malformed lines
        
    atoms = label.split('-')
    
    if section == 'BOND':
        # Key: Sorted tuple of 2 atoms (order doesn't matter)
        if len(atoms) < 2: return None
        return tuple(sorted(atoms[:2]))
        
    elif section == 'ANGL':
        # Key: (min(a,c), b, max(a,c)) -> Center b is fixed, outer atoms sorted
        if len(atoms) < 3: return None
        a, b, c = atoms[:3]
        return (min(a, c), b, max(a, c))
        
    elif section in ['DIHE', 'IMPR']:
        # Key: min( (a,b,c,d), (d,c,b,a) ) -> Order matters, but reversible
        if len(atoms) < 4: return None
        fwd = tuple(atoms[:4])
        rev = tuple(atoms[:4][::-1])
        return min(fwd, rev)
        
    elif section == 'NONB':
        return parts[0]
        
    return None

def merge_ligand_params(final_params, ligand_frcmod_path):
    """
    Reads a ligand frcmod file and adds parameters to final_params
    ONLY IF they are not already defined.
    """
    if not ligand_frcmod_path:
        return
    
    if not os.path.exists(ligand_frcmod_path):
        print(f"(Uo(I)oU)!!! Warning: Ligand FRCMOD '{ligand_frcmod_path}' not found. Skipping merge.")
        return

    print(f"Merging missing parameters from {ligand_frcmod_path}...")

    # 1. Index Existing Parameters in Munro
    # We parse the sets currently in final_params to see what Munro has already generated.
    existing_keys = {s: set() for s in ['MASS', 'BOND', 'ANGL', 'DIHE', 'IMPR', 'NONB']}
    
    for section, lines in final_params.items():
        for line in lines:
            key = parse_atom_key(line, section)
            if key:
                existing_keys[section].add(key)

    # 2. Parse Ligand File & Merge
    current_section = None
    
    with open(ligand_frcmod_path, 'r') as f:
        for line in f:
            line_strip = line.strip()
            if not line_strip: continue
            
            # Detect Section Headers
            if line_strip.startswith('MASS'): current_section = 'MASS'; continue
            if line_strip.startswith('BOND'): current_section = 'BOND'; continue
            if line_strip.startswith(('ANGL', 'ANGLE')): current_section = 'ANGL'; continue
            if line_strip.startswith('DIHE'): current_section = 'DIHE'; continue
            if line_strip.startswith('IMPR'): current_section = 'IMPR'; continue
            if line_strip.startswith('NONB'): current_section = 'NONB'; continue
            
            if not current_section: continue
            
            # Parse the key for this line
            key = parse_atom_key(line, current_section)
            if not key: continue
            
            # CHECK: Do we already have this parameter?
            if key not in existing_keys[current_section]:
                # If not, add it to our final set
                # We append a tag so you know it came from the ligand file
                final_params[current_section].add(f"{line.rstrip():<60}  SOURCE_LIGAND")
                
                # IMPORTANT: For Bonds/Angles, duplicate definitions are errors.
                # So we mark this key as 'seen'.
                # For Dihedrals, AMBER allows multiple terms (lines) for the same atoms.
                # If Munro has *no* definition, we want ALL terms from the ligand file.
                # So we generally do NOT block subsequent lines for Dihedrals.
                if current_section not in ['DIHE', 'IMPR']:
                    existing_keys[current_section].add(key)
