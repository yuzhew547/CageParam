#!/usr/bin/env python3
import argparse
import math
from collections import defaultdict
import sys
import os
import glob

try:
    from .gaff_utils import GaffDatabase, Mol2TypeMap            # package mode
except ImportError:
    try:
        from gaff_utils import GaffDatabase, Mol2TypeMap         # script-mode fallback
    except ImportError:
        print("Error: gaff_utils.py not found.")
        sys.exit(1)

try:
    from . import gaff_typing                                    # package mode
except ImportError:
    try:
        import gaff_typing                                       # script-mode fallback
    except ImportError:
        gaff_typing = None  # only required if --auto-from-pdb is used

def _stabilize_multiterm_dihe(sorted_lines):
    """AMBER multi-term DIHE/IMPR convention: continuation lines have PN<0,
    the terminator has PN>0, and they MUST appear continuation-first. The
    output is alphabetically sorted (which groups same atom-type signatures
    together via the leading 11-char field), but within a group the line
    with PN>0 sorts ahead of PN<0 — wrong. Re-sort within each group so
    PN<0 (largest |PN| first) precedes PN>0; otherwise tleap mis-parses and
    crashes with `Realloc: Cannot allocate memory!` in varArray.c.
    """
    out = []
    i = 0
    n = len(sorted_lines)
    while i < n:
        sig = sorted_lines[i][:11]
        j = i
        while j < n and sorted_lines[j][:11] == sig:
            j += 1
        group = sorted_lines[i:j]
        if len(group) > 1:
            def _key(line):
                parts = line[11:].split()
                try:
                    pn = float(parts[3])
                except (IndexError, ValueError):
                    return (1, 0.0)
                # PN<0 (continuation) first, larger |PN| first within the
                # continuation block; PN>0 (terminator) last.
                return (0, -abs(pn)) if pn < 0 else (1, abs(pn))
            group.sort(key=_key)
        out.extend(group)
        i = j
    return out


def merge_ligand_params(params, filepath):
    if not filepath or not os.path.exists(filepath):
        return
    print(f"  Merging from {filepath}...")
    section_map = {
        'MASS': 'MASS', 'BOND': 'BOND', 'BONDS': 'BOND',
        'ANGL': 'ANGL', 'ANGLE': 'ANGL', 'ANGLES': 'ANGL',
        'DIHE': 'DIHE', 'DIHEDRAL': 'DIHE', 'DIHEDRALS': 'DIHE',
        'IMPR': 'IMPR', 'IMPROPER': 'IMPR', 'IMPROPERS': 'IMPR',
        'NONB': 'NONB', 'NONBON': 'NONB', 'NONBONDED': 'NONB'
    }
    current_section = None
    with open(filepath, 'r') as f:
        for line in f:
            raw_line = line.strip()
            if not raw_line: continue
            first_word = raw_line.split()[0].upper()
            if first_word in section_map:
                current_section = section_map[first_word]
                continue
            if current_section:
                params[current_section].add(raw_line)

def load_parameters(filepath):
    general_params = {
        # Bond defaults: Y-{metal} (legacy/pyridine values) and Y-{metal}-{ligand class}
        'Y-M':              {'k': 95.1,   'eq': 2.06},   # fallback for unknown metals
        'Y-PD':             {'k': 92.0,   'eq': 2.06},   # pyridine N - Pd (legacy alias)
        'Y-PT':             {'k': 100.0,  'eq': 2.06},   # pyridine N - Pt (legacy alias)
        'Y-PD-pyridine':    {'k': 92.0,   'eq': 2.06},
        'Y-PT-pyridine':    {'k': 100.0,  'eq': 2.06},
        'Y-PD-imidazole':   {'k': 100.0,  'eq': 2.07},
        'Y-PT-imidazole':   {'k': 112.0,  'eq': 2.06},
        # Y-M-Y angle defaults
        'Y-M-Y':              {'k': 125.0, 'eq': '90/180'},
        'Y-PD-Y-pyridine':    {'k': 125.0, 'eq': '90/180'},
        'Y-PT-Y-pyridine':    {'k': 125.0, 'eq': '90/180'},
        'Y-PD-Y-imidazole':   {'k': 130.0, 'eq': '90/180'},
        'Y-PT-Y-imidazole':   {'k': 135.0, 'eq': '90/180'},
        # c-Y-M angle (same default for both ligand classes)
        'ca-Y-M': {'k': 150.0, 'eq': 120},
    }
    specific_params = {}
    if not filepath: return general_params, specific_params
    try:
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"): continue
                parts = line.split()
                key = parts[0]
                try:
                    k_val = float(parts[1])
                    eq_val = parts[2]
                    if eq_val != "90/180": eq_val = float(eq_val)
                    param_data = {'k': k_val, 'eq': eq_val}
                    if key in general_params: general_params[key] = param_data
                    else: specific_params[key] = param_data
                except (IndexError, ValueError): continue
    except FileNotFoundError: pass
    return general_params, specific_params

def get_active_param(keys_to_check, general_key, general_params, specific_params):
    for key in keys_to_check:
        if key in specific_params: return specific_params[key], True
    return general_params[general_key], False

def find_smallest_ring(start_atom, neighbors_dict):
    """Smallest ring containing start_atom, by BFS between each pair of its neighbors.
    Returns a list of atom names (including start_atom) or None if Y is not in a ring."""
    from collections import deque
    nbrs = neighbors_dict.get(start_atom, [])
    best = None
    for i in range(len(nbrs)):
        for j in range(i + 1, len(nbrs)):
            n1, n2 = nbrs[i], nbrs[j]
            parent = {n1: None}
            queue = deque([n1])
            while queue:
                cur = queue.popleft()
                if cur == n2: break
                for nb in neighbors_dict.get(cur, []):
                    if nb == start_atom or nb in parent: continue
                    parent[nb] = cur
                    queue.append(nb)
            if n2 in parent:
                path = []
                x = n2
                while x is not None:
                    path.append(x); x = parent[x]
                ring = [start_atom] + path
                if best is None or len(ring) < len(best):
                    best = ring
    return best

def classify_ring(start_atom, neighbors_dict, type_map):
    """Classify the ring containing start_atom (Y) into 'pyridine', 'imidazole',
    or 'unknown' based on (ring size, N count). Returns (class_name, ring_atoms)."""
    ring = find_smallest_ring(start_atom, neighbors_dict)
    if ring is None: return 'unknown', None
    n_count = sum(1 for a in ring
                  if type_map.get_type(a) and type_map.get_type(a).startswith('n'))
    size = len(ring)
    if size == 6:
        return 'pyridine', ring        # incl. pyrazine/pyrimidine — all pure-aromatic, Y=nb
    if size == 5 and n_count >= 2:
        return 'imidazole', ring       # incl. pyrazole/triazole — non-pure-aromatic, Y=nc/nd
    return 'unknown', ring

def pick_effective_base_type(ring_class, neighbor_types, base_type, fallback_types):
    """Pick Y's effective GAFF base type. Imidazole-class rings force nc/nd
    partnering with cd/cc carbons; everything else stays at base_type."""
    if ring_class == 'imidazole':
        if 'cd' in neighbor_types and 'nc' in fallback_types: return 'nc'
        if 'cc' in neighbor_types and 'nd' in fallback_types: return 'nd'
    return base_type

def find_bond_with_fallback(gaff_db, neighbor_type, base_type, fallback_types):
    """Search a bond in GAFF as (neighbor_type, base_type), iterating fallbacks.
    Reorders the fallback list so cd-neighbor prefers nc, cc-neighbor prefers nd
    (matching GAFF's cc/cd <-> nd/nc partnering). Also tries cp/cq -> ca on
    the neighbor side. Returns ((k, eq), used_type, source_tag) or (None, None, None).
    """
    equiv = {'cp': 'ca', 'cq': 'ca'}
    preferred = []
    if neighbor_type == 'cd' and 'nc' in fallback_types: preferred.append('nc')
    if neighbor_type == 'cc' and 'nd' in fallback_types: preferred.append('nd')
    ordered_fallback = preferred + [b for b in fallback_types if b not in preferred]
    candidates = [base_type] + [b for b in ordered_fallback if b != base_type]
    for bt in candidates:
        hit = gaff_db.search_bond(neighbor_type, bt)
        if hit:
            tag = "GAFF" if bt == base_type else f"GAFF_fallback({bt})"
            return hit, bt, tag
        alt = equiv.get(neighbor_type, neighbor_type)
        if alt != neighbor_type:
            hit = gaff_db.search_bond(alt, bt)
            if hit:
                tag = (f"GAFF_equiv({alt}-{bt})" if bt == base_type
                       else f"GAFF_fallback_equiv({alt}-{bt})")
                return hit, bt, tag
    return None, None, None

def find_angle_with_fallback(gaff_db, t1, t2, t3, sub_index, fallback_types):
    """Search an angle in GAFF where types[sub_index] is the placeholder for Y.

    Iterates: the user's base_type (already in types[sub_index]) -> each in
    fallback_types. The fallback list is reordered by ring context so that
    cd-neighbored Ys prefer nc and cc-neighbored Ys prefer nd (the standard
    imidazole cc/cd <-> nd/nc partnering in GAFF). On each attempt also tries
    cp/cq -> ca equivalences on the non-placeholder slots. Returns
    ((k, eq), source_tag) or (None, None).
    """
    types = [t1, t2, t3]
    primary = types[sub_index]
    equiv = {'cp': 'ca', 'cq': 'ca'}
    other_types = [t for k, t in enumerate(types) if k != sub_index]
    preferred = []
    if 'cd' in other_types and 'nc' in fallback_types: preferred.append('nc')
    if 'cc' in other_types and 'nd' in fallback_types: preferred.append('nd')
    ordered_fallback = preferred + [b for b in fallback_types if b not in preferred]
    candidates = [primary] + [b for b in ordered_fallback if b != primary]
    for bt in candidates:
        trial = list(types); trial[sub_index] = bt
        hit = gaff_db.search_angle(*trial)
        if hit:
            tag = "GAFF" if bt == primary else f"GAFF_fallback({bt})"
            return hit, tag
        alt = list(trial); changed = False
        for k in range(3):
            if k == sub_index: continue
            new = equiv.get(alt[k], alt[k])
            if new != alt[k]:
                alt[k] = new; changed = True
        if changed:
            hit = gaff_db.search_angle(*alt)
            if hit:
                fmt = "-".join(alt)
                etag = (f"GAFF_equiv({fmt})" if bt == primary
                        else f"GAFF_fallback_equiv({fmt})")
                return hit, etag
    return None, None

def dist(c1, c2):
    return math.sqrt(sum((c1[i]-c2[i])**2 for i in range(3)))

def calculate_angle(c_a, c_b, c_c):
    v1 = [c_a[i]-c_b[i] for i in range(3)]
    v2 = [c_c[i]-c_b[i] for i in range(3)]
    dot = sum(v1[i]*v2[i] for i in range(3))
    mag = math.sqrt(sum(v**2 for v in v1)) * math.sqrt(sum(v**2 for v in v2))
    return 0.0 if mag == 0 else math.degrees(math.acos(max(min(dot/mag, 1.0), -1.0)))

def get_unique_2char_code(index, prefixes):
    chars = [str(i) for i in range(0, 10)] + \
            [chr(i) for i in range(ord('A'), ord('Z')+1)] + \
            [chr(i) for i in range(ord('a'), ord('z')+1)]
    if index < 1: return "XX"
    real_idx = index - 1
    prefix_idx = real_idx // len(chars)
    char_idx = real_idx % len(chars)
    if prefix_idx >= len(prefixes): return "XX"
    return f"{prefixes[prefix_idx]}{chars[char_idx]}"

def build_organic_connectivity(all_atoms, cutoff=1.9):
    adj = defaultdict(lambda: defaultdict(list))
    atoms_by_res = defaultdict(list)
    for atom in all_atoms:
        atoms_by_res[atom['res_key']].append(atom)
    for res_key, atoms in atoms_by_res.items():
        n = len(atoms)
        for i in range(n):
            for j in range(i + 1, n):
                a1, a2 = atoms[i], atoms[j]
                if dist(a1['coords'], a2['coords']) < cutoff:
                    adj[res_key][a1['name']].append(a2['name'])
                    adj[res_key][a2['name']].append(a1['name'])
    return adj

class LigandTypeManager:
    def __init__(self, gaff_db, base_type="nb"):
        self.gaff_db = gaff_db
        self.base_type = base_type
        self.type_maps = {}
        
    def add_ligand_type(self, mol2_file, resname_pattern=None):
        if resname_pattern is None:
            basename = os.path.splitext(os.path.basename(mol2_file))[0]
            if basename.startswith('L'):
                if len(basename) >= 2 and basename[1].isalpha():
                    resname_pattern = basename[:2]
                else:
                    resname_pattern = basename[:1]
            else:
                resname_pattern = "L"
        print(f"  Loading ligand type '{resname_pattern}' from {mol2_file}")
        self.type_maps[resname_pattern] = Mol2TypeMap(mol2_file)
        
    def get_type_map(self, resname):
        if resname in self.type_maps: return self.type_maps[resname]
        if len(resname) >= 2:
            if resname[1].isalpha(): prefix = resname[:2]
            else: prefix = resname[:1]
        else: prefix = resname
        if prefix in self.type_maps: return self.type_maps[prefix]
        for pattern, type_map in self.type_maps.items():
            if resname.startswith(pattern): return type_map
        print(f"  WARNING: No type map found for residue '{resname}'")
        return None

def main():
    parser = argparse.ArgumentParser(description="MUNRO Multi-Template")
    parser.add_argument("-p", "--pdb", required=True)
    parser.add_argument("-g", "--gaff", default=None,
                        help="Path to gaff2.dat or gaff.dat. "
                             "Defaults to the gaff2.dat bundled with cagepipe.")
    parser.add_argument("-o", "--output", default="munro.frcmod")
    ligand_group = parser.add_mutually_exclusive_group(required=True)
    ligand_group.add_argument("-m", "--mol2")
    ligand_group.add_argument("--ligands", nargs='+')
    ligand_group.add_argument("--auto", action='store_true')
    ligand_group.add_argument("--auto-from-pdb", nargs='*', metavar="TEMPLATE_PDB",
                              help="Run antechamber + parmchk2 on template PDB(s) "
                                   "(or auto-glob *template*.pdb if no files given), "
                                   "then load the resulting .mol2/.frcmod automatically.")
    parser.add_argument("-c", "--custom")
    parser.add_argument("--base-type", default="nb", help="GAFF type that Y represents (default: nb for pyridine N)")
    parser.add_argument("--fallback-types", nargs='*', default=['nc', 'nd'],
                        help="GAFF base-type fallbacks tried when --base-type yields no "
                             "angle match (default: nc nd, suitable for imidazole-type N).")
    parser.add_argument("-l", "--ligand-frcmod", dest='single_ligand_frcmod')
    parser.add_argument("--merge-ligands", nargs='+')
    parser.add_argument("--net-charge", type=int, default=0,
                        help="Default net charge passed to antechamber (default: 0)")
    parser.add_argument("--charge", nargs='+', default=[], metavar="PREFIX:CHARGE",
                        help="Per-prefix net charges, e.g. LA:0 LB:-1")
    parser.add_argument("--atom-type", default=None,
                        help="Atom type style for antechamber/parmchk2 "
                             "(default: auto-derive from -g basename — "
                             "gaff2.dat -> gaff2, gaff.dat -> gaff).")
    parser.add_argument("--force-typing", action="store_true",
                        help="Re-run antechamber/parmchk2 even if cached outputs are fresh")
    args = parser.parse_args()

    # Default --gaff to the bundled gaff2.dat (works both in script and package modes).
    if args.gaff is None:
        try:
            from importlib.resources import files as _ir_files
            args.gaff = str(_ir_files("cagepipe").joinpath("data", "gaff2.dat"))
        except Exception:
            args.gaff = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                     "data", "gaff2.dat")
        print(f"  -g not given; using bundled GAFF: {args.gaff}")

    # Auto-derive antechamber/parmchk2 atom-type style from the gaff data file
    # name unless the user passed --atom-type explicitly.
    if args.atom_type is None:
        gaff_base = os.path.basename(args.gaff).lower()
        args.atom_type = "gaff2" if "gaff2" in gaff_base else "gaff"
        print(f"  --atom-type auto-set to '{args.atom_type}' "
              f"(from -g {os.path.basename(args.gaff)})")

    print("="*70)
    print("MUNRO Multi-Template - Force Field Parameter Generator")
    print("="*70)

    print("\n1. Loading GAFF Database:")
    gaff_db = GaffDatabase(args.gaff)
    gen_params, spec_params = load_parameters(args.custom)

    print("\n2. Loading Ligand Templates:")
    ligand_mgr = LigandTypeManager(gaff_db, args.base_type)
    auto_generated_frcmods = []  # files produced by gaff_typing, queued for merge
    if args.mol2:
        ligand_mgr.add_ligand_type(args.mol2, "L")
    elif args.ligands:
        for mol2_file in args.ligands:
            ligand_mgr.add_ligand_type(mol2_file)
    elif args.auto:
        mol2_files = sorted(glob.glob("L*.mol2"))
        if not mol2_files:
            print("  ERROR: No L*.mol2 files found!"); sys.exit(1)
        for mol2_file in mol2_files:
            ligand_mgr.add_ligand_type(mol2_file)
    elif args.auto_from_pdb is not None:
        if gaff_typing is None:
            print("  ERROR: gaff_typing.py is required for --auto-from-pdb"); sys.exit(1)
        templates = list(args.auto_from_pdb)
        if not templates:
            templates = gaff_typing.discover_templates()
        if not templates:
            print("  ERROR: No template PDB files supplied or found "
                  "(looked for *template*.pdb)"); sys.exit(1)
        charges_by_prefix = {}
        for spec in args.charge:
            if ":" in spec:
                k, v = spec.split(":", 1)
                charges_by_prefix[k] = int(v)
        print(f"  Running antechamber + parmchk2 on {len(templates)} template(s):")
        for t in templates:
            print(f"    - {t}")
        results = gaff_typing.generate_for_templates(
            templates,
            output_dir=".",
            net_charge=args.net_charge,
            charges_by_prefix=charges_by_prefix,
            atom_type=args.atom_type,
            force=args.force_typing,
        )
        for mol2_file, frcmod_file, prefix in results:
            ligand_mgr.add_ligand_type(mol2_file, prefix)
            auto_generated_frcmods.append(frcmod_file)

    # PARSE PDB
    print(f"\n3. Reading PDB: {args.pdb}")
    pd_residues, ml_residues = [], []
    pd_atoms, n_atoms, all_ligand_atoms = [], [], []
    metal_types = {}
    ligand_resnames = set()

    with open(args.pdb, 'r') as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                res_key = (line[21].strip() or "-", int(line[22:26].strip()))
                name = line[12:16].strip()
                coords = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
                res_name = line[17:20].strip().upper()
                atom_data = {'res_key': res_key, 'name': name, 'coords': coords, 'res_name': res_name}
                if res_name.startswith("P") and not res_name.startswith("PRO"):
                    if res_key not in pd_residues: pd_residues.append(res_key)
                    pd_atoms.append(atom_data)
                    metal_types[res_key] = name.upper()
                elif res_name.startswith("L") and not res_name.startswith(("LEU","LYS")):
                    if res_key not in ml_residues: ml_residues.append(res_key)
                    all_ligand_atoms.append(atom_data)
                    ligand_resnames.add(res_name)
                    if "N" in name: n_atoms.append(atom_data)

    pd_map = {k: i+1 for i, k in enumerate(pd_residues)}
    ml_map = {k: i+1 for i, k in enumerate(ml_residues)}
    print(f"  Metals: {len(pd_atoms)} in {len(pd_residues)} residues")
    print(f"  Ligands: {len(all_ligand_atoms)} in {len(ml_residues)} residues")
    print(f"  Unique ligand types: {sorted(ligand_resnames)}")

    # BUILD CONNECTIVITY
    print("\n4. Building Connectivity:")
    adj_list = build_organic_connectivity(all_ligand_atoms, cutoff=1.9)
    
    y_atoms_map = {}
    pd_connectivity = defaultdict(list)
    coord_candidates = defaultdict(list)
    for pd in pd_atoms:
        for n in n_atoms:
            if dist(pd['coords'], n['coords']) < 3.0:
                l_idx = ml_map[n['res_key']]
                coord_candidates[l_idx].append(n)

    y_counter = 1
    for l_idx in sorted(coord_candidates.keys()):
        atoms = sorted(coord_candidates[l_idx], key=lambda x: x['name'])
        for atom in atoms:
            uid = (atom['res_key'], atom['name'])
            if uid not in y_atoms_map:
                y_code = get_unique_2char_code(y_counter, "YZWVU")
                y_atoms_map[uid] = {'code': y_code, 'data': atom}
                y_counter += 1

    for pd in pd_atoms:
        m_idx = pd_map[pd['res_key']]
        for n in n_atoms:
            uid = (n['res_key'], n['name'])
            if uid in y_atoms_map and dist(pd['coords'], n['coords']) < 3.0:
                info = y_atoms_map[uid]
                pd_connectivity[m_idx].append({'y_code': info['code'], 'coords': n['coords'], 'data': n})

    print(f"  Y atoms: {len(y_atoms_map)}")
    print(f"  Metal centers: {len(pd_connectivity)}")

    # Classify each Y by its ring and lock in its effective GAFF base type.
    y_ring_class = {}      # y_code -> 'pyridine' | 'imidazole' | 'unknown'
    y_effective_base = {}  # y_code -> GAFF base used for bond/angle GAFF lookups
    for uid, info in y_atoms_map.items():
        y_code = info['code']
        res_key, y_name = uid
        res_name = info['data']['res_name']
        type_map = ligand_mgr.get_type_map(res_name)
        if type_map is None:
            y_ring_class[y_code] = 'unknown'
            y_effective_base[y_code] = args.base_type
            continue
        cls, _ring = classify_ring(y_name, adj_list[res_key], type_map)
        neighbor_types = [type_map.get_type(n) for n in adj_list[res_key].get(y_name, [])]
        y_ring_class[y_code] = cls
        y_effective_base[y_code] = pick_effective_base_type(
            cls, neighbor_types, args.base_type, args.fallback_types)
    cls_counts = defaultdict(int)
    for cls in y_ring_class.values(): cls_counts[cls] += 1
    print(f"  Ring classes: {dict(cls_counts)}")

    # GENERATE PARAMETERS
    print("\n5. Generating Parameters:")
    params = {'MASS': set(), 'BOND': set(), 'ANGL': set(), 'DIHE': set(), 'IMPR': set(), 'NONB': set()}

    # MASS
    for res_key, idx in pd_map.items():
        m_code = get_unique_2char_code(idx, 'MNO')
        el_type = metal_types.get(res_key, "PD")
        mass_val = "195.08" if "PT" in el_type else "106.42"
        mass_comm = "Pt ion" if "PT" in el_type else "Pd ion"
        params['MASS'].add(f"{m_code:<2} {mass_val:<20} {mass_comm}")

    for uid, info in y_atoms_map.items():
        params['MASS'].add(f"{info['code']:<2} 14.01         0.530               Sp2 N")

    # BOND
    inv_pd_map = {v: k for k, v in pd_map.items()}
    for m_idx, y_list in pd_connectivity.items():
        m_code = get_unique_2char_code(m_idx, "MNO")
        metal_el = metal_types.get(inv_pd_map.get(m_idx), "PD").upper()
        for y_item in y_list:
            yc = y_item['y_code']
            cls = y_ring_class.get(yc, 'unknown')
            candidates = [f"Y-{metal_el}-{cls}", f"Y-{metal_el}", "Y-M"]
            general_key = next((c for c in candidates if c in gen_params), "Y-M")
            p, is_spec = get_active_param([f"{yc}-{m_code}", f"{m_code}-{yc}"],
                                          general_key, gen_params, spec_params)
            label = ("User-defined" if is_spec
                     else (f"Default({cls})" if cls in ('pyridine','imidazole') else "Default"))
            params['BOND'].add(f"{yc:<2}-{m_code:<2}   {p['k']:<5}   {p['eq']:<7}     {label}")

    for uid, info in y_atoms_map.items():
        y_code = info['code']
        res_key = uid[0]
        y_name = uid[1]
        res_name = info['data']['res_name']
        type_map = ligand_mgr.get_type_map(res_name)
        if type_map is None: continue
        neighbors = adj_list[res_key].get(y_name, [])
        neighbor_types = [type_map.get_type(n) for n in neighbors]
        eff_base = y_effective_base.get(y_code, args.base_type)
        ctx_note = "" if eff_base == args.base_type else f"_ctx({eff_base})"
        for n_name, t_n in zip(neighbors, neighbor_types):
            param, _, source = find_bond_with_fallback(
                gaff_db, t_n, eff_base, args.fallback_types)
            if param:
                k, eq = param
                params['BOND'].add(f"{t_n:<2}-{y_code:<2}   {k:<5.1f}   {eq:<6.4f}    {source}{ctx_note}")
            else:
                params['BOND'].add(f"{t_n:<2}-{y_code:<2}   488.0   1.339     GAFF_Ref({t_n}-{eff_base})")
                print(f"  WARNING: No GAFF bond for {t_n}-{eff_base} "
                      f"(fallbacks tried: {args.fallback_types}) "
                      f"({n_name}-{y_name} in {res_name}); using default 488.0/1.339")

    # ANGLE
    for m_idx, y_list in pd_connectivity.items():
        m_code = get_unique_2char_code(m_idx, "MNO")
        metal_el = metal_types.get(inv_pd_map.get(m_idx), "PD").upper()
        m_coords = next(p['coords'] for p in pd_atoms if pd_map[p['res_key']] == m_idx)
        for i in range(len(y_list)):
            for j in range(i + 1, len(y_list)):
                y1, y2 = y_list[i], y_list[j]
                cls1 = y_ring_class.get(y1['y_code'], 'unknown')
                cls2 = y_ring_class.get(y2['y_code'], 'unknown')
                common_cls = cls1 if cls1 == cls2 else 'mixed'
                candidates = []
                if common_cls in ('pyridine', 'imidazole'):
                    candidates.append(f"Y-{metal_el}-Y-{common_cls}")
                candidates.append("Y-M-Y")
                general_key = next((c for c in candidates if c in gen_params), "Y-M-Y")
                p, is_spec = get_active_param(
                    [f"{y1['y_code']}-{m_code}-{y2['y_code']}"],
                    general_key, gen_params, spec_params)
                eq = 90.0
                if p['eq'] == "90/180":
                    ang = calculate_angle(y1['coords'], m_coords, y2['coords'])
                    eq = 180.0 if ang > 135 else 90.0
                else: eq = float(p['eq'])
                label = ("User-defined" if is_spec
                         else (f"Default({common_cls})" if common_cls in ('pyridine','imidazole') else "Default"))
                params['ANGL'].add(f"{y1['y_code']:<2}-{m_code:<2}-{y2['y_code']:<2}   {p['k']:.2f}     {eq:.2f}    {label}")

    # Y-M-neighbor angles
    for uid, info in y_atoms_map.items():
        y_code = info['code']
        res_key = uid[0]
        y_name = uid[1]
        m_code = next(get_unique_2char_code(m, "MNO") for m, ys in pd_connectivity.items() if any(y['y_code']==y_code for y in ys))
        res_name = info['data']['res_name']
        type_map = ligand_mgr.get_type_map(res_name)
        if type_map is None: continue
        neighbors = adj_list[res_key].get(y_name, [])
        for n_name in neighbors:
            t_n = type_map.get_type(n_name)
            p, is_spec = get_active_param([f"{t_n}-{y_code}-{m_code}"], "ca-Y-M", gen_params, spec_params)
            params['ANGL'].add(f"{t_n:<2}-{y_code:<2}-{m_code:<2}   {p['k']:.2f}     {p['eq']:.2f}    {'User-defined' if is_spec else 'Default'}")

    # Internal angles - THIS IS WHERE THE BUG IS
    for uid, info in y_atoms_map.items():
        y_code = info['code']
        res_key = uid[0]
        y_name = uid[1]
        res_name = info['data']['res_name']
        type_map = ligand_mgr.get_type_map(res_name)
        if type_map is None: continue
        neighbors = adj_list[res_key].get(y_name, [])
        
        eff_base = y_effective_base.get(y_code, args.base_type)
        ctx_note = "" if eff_base == args.base_type else f"_ctx({eff_base})"

        # Angles around Y (Y is the center)
        for i in range(len(neighbors)):
            for j in range(i+1, len(neighbors)):
                n1, n2 = neighbors[i], neighbors[j]
                t1, t2 = type_map.get_type(n1), type_map.get_type(n2)
                param, source = find_angle_with_fallback(
                    gaff_db, t1, eff_base, t2, 1, args.fallback_types)
                if param:
                    params['ANGL'].add(f"{t1:<2}-{y_code:<2}-{t2:<2}   {param[0]:.2f}     {param[1]:.2f}    {source}{ctx_note}")
                else:
                    # Default: aromatic-edge angle ~70 kcal/mol/rad^2, 120 deg.
                    # tleap rejects the unit if any angle is undefined, so we
                    # emit a sane default instead of dropping the term.
                    params['ANGL'].add(f"{t1:<2}-{y_code:<2}-{t2:<2}   70.00     120.00    Default_no_GAFF({t1}-{eff_base}-{t2})")
                    print(f"  WARNING: No GAFF angle for {t1}-{eff_base}-{t2} "
                          f"(fallbacks tried: {args.fallback_types}) "
                          f"({n1}-{y_name}-{n2} in {res_name}); using default 70.0/120.0")

        # Angles where Y is an endpoint (tnn-t1-Y)
        for n1 in neighbors:
            t1 = type_map.get_type(n1)
            nn_list = adj_list[res_key].get(n1, [])
            for nn in nn_list:
                if nn == y_name: continue
                tnn = type_map.get_type(nn)
                param, source = find_angle_with_fallback(
                    gaff_db, tnn, t1, eff_base, 2, args.fallback_types)
                if param:
                    params['ANGL'].add(f"{tnn:<2}-{t1:<2}-{y_code:<2}   {param[0]:.2f}     {param[1]:.2f}    {source}{ctx_note}")
                else:
                    params['ANGL'].add(f"{tnn:<2}-{t1:<2}-{y_code:<2}   70.00     120.00    Default_no_GAFF({tnn}-{t1}-{eff_base})")
                    print(f"  WARNING: No GAFF angle for {tnn}-{t1}-{eff_base} "
                          f"(fallbacks tried: {args.fallback_types}) "
                          f"({nn}-{n1}-{y_name} in {res_name}); using default 70.0/120.0")

    # DIHE
    z_param = "1    0.00          0.00   2.0      Generic_M_Interface"
    for uid, info in y_atoms_map.items():
        yc = info['code']
        res_key = uid[0]
        y_name = uid[1]
        m_code = next(get_unique_2char_code(m, "MNO") for m, ys in pd_connectivity.items() if any(y['y_code']==yc for y in ys))
        params['DIHE'].add(f"X -{yc:<2}-{m_code:<2}-X    {z_param}")
        params['DIHE'].add(f"X -{m_code:<2}-{yc:<2}-X    {z_param}")
        res_name = info['data']['res_name']
        type_map = ligand_mgr.get_type_map(res_name)
        if type_map is None: continue
        neighbors = adj_list[res_key].get(y_name, [])
        for n in neighbors:
            t_n = type_map.get_type(n)
            params['DIHE'].add(f"X -{t_n:<2}-{yc:<2}-X    {z_param}")
            params['DIHE'].add(f"X -{yc:<2}-{t_n:<2}-X    {z_param}")

    # IMPR
    for uid, info in y_atoms_map.items():
        y_code = info['code']
        res_key = uid[0]
        y_name = uid[1]
        res_name = info['data']['res_name']
        type_map = ligand_mgr.get_type_map(res_name)
        if type_map is None: continue
        c_alphas = adj_list[res_key].get(y_name, [])
        for ca_name in c_alphas:
            ca_neighbors = adj_list[res_key].get(ca_name, [])
            candidates = [n for n in ca_neighbors if n != y_name]
            if len(candidates) < 2: continue
            t1 = type_map.get_type(candidates[0])
            t2 = type_map.get_type(candidates[1])
            subst_types = ('h', 'f', 'cl', 'br', 'c3', 'oh', 'os')
            if t1.startswith(subst_types) and not t2.startswith(subst_types):
                x_type, ring_type = t1, t2
            elif t2.startswith(subst_types) and not t1.startswith(subst_types):
                x_type, ring_type = t2, t1
            else:
                x_type, ring_type = t1, t2
            t_ca = type_map.get_type(ca_name)
            param = gaff_db.search_improper(args.base_type, ring_type, t_ca, x_type)
            if param:
                k, ph, per = param
                params['IMPR'].add(f"{y_code:<2}-{ring_type:<2}-{t_ca:<2}-{x_type:<2}   {k:<5}        {ph:<5}        {per:<5}  GAFF")
            else:
                params['IMPR'].add(f"{y_code:<2}-{ring_type:<2}-{t_ca:<2}-{x_type:<2}   1.1          180.0        2.0    Default")

    # VDW
    for res_key, idx in pd_map.items():
        m_code = get_unique_2char_code(idx, 'MNO')
        el_type = metal_types.get(res_key, "PD")
        if "PT" in el_type: r_vdw, e_vdw, comm = 1.2190, 0.0015090300, "Pt2+"
        else: r_vdw, e_vdw, comm = 1.2690, 0.0032106800, "Pd2+"
        params['NONB'].add(f"  {m_code:<2}           {r_vdw:.4f}  {e_vdw:.10f}       {comm}")

    vdw_r, vdw_e = 1.8240, 0.1700
    for uid, info in y_atoms_map.items():
        params['NONB'].add(f"  {info['code']:<2}           {vdw_r:.4f}  {vdw_e:.4f}             GAFF({args.base_type})")

    # WRITE OUTPUT
    # 6. Merge ligand frcmod files
    frcmod_files = []
    if args.single_ligand_frcmod:
        frcmod_files.append(args.single_ligand_frcmod)
    if args.merge_ligands:
        frcmod_files.extend(args.merge_ligands)
    if auto_generated_frcmods:
        frcmod_files.extend(auto_generated_frcmods)
    if frcmod_files:
        print("\n6. Merging Ligand FRCMOD files:")
        for frcmod_path in frcmod_files:
            merge_ligand_params(params, frcmod_path)

    print(f"\n7. Writing Output: {args.output}")
    with open(args.output, 'w') as f:
        f.write(f"REMARK GENERATED BY MUNRO MULTI-TEMPLATE v2.0\n")
        for section in ['MASS', 'BOND', 'ANGL', 'DIHE', 'IMPR', 'NONB']:
            if params[section]:
                f.write(f"{section}\n")
                sorted_lines = sorted(list(params[section]))
                if section in ('DIHE', 'IMPR'):
                    sorted_lines = _stabilize_multiterm_dihe(sorted_lines)
                for line in sorted_lines:
                    f.write(line + "\n")
                f.write("\n")

    print("\nSummary:")
    for s in ['MASS','BOND','ANGL','DIHE','IMPR','NONB']:
        print(f"  {s}: {len(params[s])}")

if __name__ == "__main__":
    main()
