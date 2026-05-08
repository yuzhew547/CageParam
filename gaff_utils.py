class GaffDatabase:
    def __init__(self, gaff_dat_path):
        self.angles = {}  # Key: (t1, center, t2), Value: (k, eq)
        self.impropers = {} # Key: tuple(sorted(atoms)), Value: (k, phase, per)
        if gaff_dat_path:
            self.load_database(gaff_dat_path)

    def load_database(self, path):
        print(f"Loading GAFF parameters from {path}...")
        with open(path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line: continue
                parts = line.split()
                if not parts: continue

                # --- ANGLE (3 atoms) ---
                if '-' in parts[0] and len(parts[0].split('-')) == 3:
                    try:
                        sub = parts[0].split('-')
                        self.add_angle(sub[0], sub[1], sub[2], float(parts[1]), float(parts[2]))
                    except: pass
                
                # --- IMPROPER (4 atoms) ---
                elif '-' in parts[0] and len(parts[0].split('-')) == 4:
                    try:
                        sub = parts[0].split('-')
                        self.add_improper(sub, float(parts[1]), float(parts[2]), float(parts[3]))
                    except: pass

                # --- STANDARD SPACE SEPARATED ---
                elif len(parts) >= 5:
                    if len(parts) == 5: 
                         try:
                             self.add_angle(parts[0], parts[1], parts[2], float(parts[3]), float(parts[4]))
                         except: pass
                    elif len(parts) >= 7:
                        try:
                            atoms = [parts[0], parts[1], parts[2], parts[3]]
                            self.add_improper(atoms, float(parts[4]), float(parts[5]), float(parts[6]))
                        except: pass

    def add_angle(self, t1, c, t2, k, eq):
        self.angles[(t1, c, t2)] = (k, eq)
        self.angles[(t2, c, t1)] = (k, eq)

    def add_improper(self, atoms, k, phase, per):
        self.impropers[tuple(atoms)] = (k, phase, per)

    def search_angle(self, t1, center, t2):
        return self.angles.get((t1, center, t2))

    def search_improper(self, a, b, c, d):
        # Try exact match
        key = (a, b, c, d)
        if key in self.impropers: return self.impropers[key]
        
        # Try Wildcards (X) - standard AMBER improper often has X-X-c-d
        candidates = [
            ("X", "X", c, d),
            ("X", b, c, d),
            (a, "X", c, d),
            (a, b, c, "X")
        ]
        for cand in candidates:
            if cand in self.impropers: return self.impropers[cand]
        return None

class Mol2TypeMap:
    def __init__(self, mol2_path):
        self.name_to_type = {} 
        if mol2_path: self.load_mol2(mol2_path)

    def load_mol2(self, path):
        print(f"Loading Atom Types from {path}...")
        parsing_atoms = False  # Standardized variable name
        try:
            with open(path, 'r') as f:
                for line in f:
                    if line.startswith("@<TRIPOS>ATOM"):
                        parsing_atoms = True
                        continue
                    elif line.startswith("@<TRIPOS>BOND"):
                        parsing_atoms = False
                        continue
                    
                    if parsing_atoms: # Now correctly references parsing_atoms
                        parts = line.split()
                        if len(parts) >= 6:
                            self.name_to_type[parts[1]] = parts[5]
        except FileNotFoundError:
            print(f"Error: Mol2 file {path} not found.")

    def get_type(self, atom_name):
        return self.name_to_type.get(atom_name, "??")
