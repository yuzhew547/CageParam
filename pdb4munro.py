#!/usr/bin/env python3
"""
pdb4munro.py - Auto-detect ligands and produce bone.pdb

Reads a system file (.pdb or .xyz; XYZ is converted via OpenBabel),
identifies discrete ligand fragments by graph connectivity (the BOND_CUTOFF
of 1.90 A naturally excludes metal-ligand bonds such as Pd-N at ~2.0 A),
clusters the fragments into unique types via graph isomorphism, and writes:
  - <prefix>tempK_template.pdb for each unique ligand type
  - bone.pdb with metals separated and ligands renamed
"""
import os
import sys
import argparse
import shutil
import subprocess
from collections import defaultdict

# ==================== CONFIGURATION ====================
BOND_CUTOFF = 1.90  # Angstrom; organic bonds, Pd-N (~2.0 A) excluded
METAL_ELEMENTS = {"PD", "PT", "AU", "AG", "NI", "CU", "ZN", "FE", "CO", "RU", "RH", "IR"}
OBABEL_CANDIDATES = [
    "/home/gridsan/ywang6/sft/build/bin/obabel",
    "obabel",
    "babel",
]
# =======================================================


# --------------- File loading (PDB / XYZ) ---------------

def find_obabel():
    for p in OBABEL_CANDIDATES:
        if os.path.isfile(p) and os.access(p, os.X_OK):
            return p
        which = shutil.which(p)
        if which:
            return which
    return None


def obabel_env():
    """Build env vars so the local OpenBabel build can find its plugins."""
    env = os.environ.copy()
    libdir = "/home/gridsan/ywang6/sft/build/lib"
    datadir = "/home/gridsan/ywang6/sft/openbabel-openbabel-2-4-0/data"
    if os.path.isdir(libdir):
        env.setdefault("BABEL_LIBDIR", libdir)
    if os.path.isdir(datadir):
        env.setdefault("BABEL_DATADIR", datadir)
    return env


def xyz_to_pdb(xyz_path, pdb_path):
    obabel = find_obabel()
    if obabel is None:
        raise RuntimeError(
            "OpenBabel (obabel) not found. Install it or update OBABEL_CANDIDATES."
        )
    print(f"   Converting {xyz_path} -> {pdb_path}\n   using {obabel}")
    # OpenBabel CLI needs `-O<file>` joined, not separated.
    res = subprocess.run(
        [obabel, xyz_path, f"-O{pdb_path}"],
        capture_output=True, text=True, env=obabel_env(),
    )
    if res.returncode != 0 or not os.path.isfile(pdb_path):
        raise RuntimeError(
            f"obabel failed (rc={res.returncode}):\n"
            f"STDOUT: {res.stdout}\nSTDERR: {res.stderr}"
        )


def parse_pdb(filename):
    """Return list of atom dicts parsed from a PDB file."""
    atoms = []
    with open(filename, 'r') as f:
        for line in f:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue
            element = line[76:78].strip().upper()
            atom_name = line[12:16].strip()
            if not element:
                element = ''.join(c for c in atom_name if c.isalpha()).upper()
            atoms.append({
                'element': element,
                'name': atom_name,
                'x': x, 'y': y, 'z': z,
                'orig_idx': len(atoms),
            })
    return atoms


def load_input(filename):
    """Load atoms from .pdb or .xyz (xyz routed through OpenBabel)."""
    ext = os.path.splitext(filename)[1].lower()
    if ext == ".pdb":
        return parse_pdb(filename)
    if ext == ".xyz":
        base = os.path.splitext(os.path.basename(filename))[0]
        pdb_path = f"{base}_obabel.pdb"
        xyz_to_pdb(filename, pdb_path)
        return parse_pdb(pdb_path)
    raise ValueError(f"Unsupported file extension '{ext}'. Use .pdb or .xyz")


# --------------- Connectivity / isomorphism ---------------

def get_dist_sq(a1, a2):
    return (a1['x'] - a2['x'])**2 + (a1['y'] - a2['y'])**2 + (a1['z'] - a2['z'])**2


def build_adjacency(atoms):
    """Brute-force adjacency by distance cutoff."""
    adj = defaultdict(list)
    cutoff_sq = BOND_CUTOFF**2
    n = len(atoms)
    for i in range(n):
        for j in range(i + 1, n):
            if get_dist_sq(atoms[i], atoms[j]) < cutoff_sq:
                adj[i].append(j)
                adj[j].append(i)
    return adj


def get_signature(atom_idx, atoms, adj):
    el = atoms[atom_idx]['element']
    n_els = sorted(atoms[n]['element'] for n in adj[atom_idx])
    return (el, len(n_els), tuple(n_els))


def solve_isomorphism(templ_atoms, templ_adj, target_atoms, target_adj):
    """Backtracking search for a 1:1 atom mapping templ -> target."""
    if len(templ_atoms) != len(target_atoms):
        return None

    templ_sigs = [get_signature(i, templ_atoms, templ_adj) for i in range(len(templ_atoms))]
    target_sigs = [get_signature(i, target_atoms, target_adj) for i in range(len(target_atoms))]

    mapping = {}
    used_targets = set()

    def backtrack(t_idx):
        if t_idx == len(templ_atoms):
            return True
        t_sig = templ_sigs[t_idx]
        for cand in range(len(target_atoms)):
            if cand in used_targets or target_sigs[cand] != t_sig:
                continue
            ok = True
            for nbr in templ_adj[t_idx]:
                if nbr in mapping and mapping[nbr] not in target_adj[cand]:
                    ok = False
                    break
            if ok:
                mapping[t_idx] = cand
                used_targets.add(cand)
                if backtrack(t_idx + 1):
                    return True
                del mapping[t_idx]
                used_targets.remove(cand)
        return False

    return mapping if backtrack(0) else None


# --------------- Ligand fragment discovery ---------------

def find_connected_components(atoms, adj):
    visited = set()
    components = []
    for i in range(len(atoms)):
        if i in visited:
            continue
        q = [i]
        visited.add(i)
        comp = []
        while q:
            curr = q.pop(0)
            comp.append(curr)
            for n in adj[curr]:
                if n not in visited:
                    visited.add(n)
                    q.append(n)
        components.append(comp)
    return components


def auto_rename_by_element(atoms):
    """Rename atoms in place to C1, C2, N1, ..., H1, ..."""
    counts = defaultdict(int)
    for atom in atoms:
        el = atom['element'].upper()
        prefix = el.title() if len(el) > 1 else el
        counts[el] += 1
        atom['name'] = f"{prefix}{counts[el]}"


def discover_templates(ligand_fragments):
    """
    ligand_fragments: list of (atoms_list, local_adj)

    Cluster fragments by graph isomorphism. Each unique fragment becomes a
    template; subsequent matching fragments are mapped onto it.

    Returns:
        templates   - list of {index, name, atoms, adj, num_atoms, src_ligand}
        assignments - list aligned with ligand_fragments;
                      each entry = (template_index, mapping[t_idx -> local_idx])
    """
    templates = []
    assignments = []
    for lig_idx, (atoms, adj) in enumerate(ligand_fragments):
        matched = False
        for tpl in templates:
            mapping = solve_isomorphism(tpl['atoms'], tpl['adj'], atoms, adj)
            if mapping is not None:
                assignments.append((tpl['index'], mapping))
                matched = True
                break
        if matched:
            continue
        # Seed a new template from this fragment.
        tpl_atoms = [dict(a) for a in atoms]
        auto_rename_by_element(tpl_atoms)
        tpl_adj = build_adjacency(tpl_atoms)
        idx = len(templates)
        templates.append({
            'index': idx,
            'name': f"temp{idx + 1}",
            'atoms': tpl_atoms,
            'adj': tpl_adj,
            'num_atoms': len(tpl_atoms),
            'src_ligand': lig_idx,
        })
        assignments.append((idx, {i: i for i in range(len(atoms))}))
    return templates, assignments


# --------------- PDB writing helpers ---------------

def fmt_atom_name(name):
    clean = name.strip()
    return f" {clean:<3}" if len(clean) < 4 else f"{clean:<4}"


def write_template_pdb(filename, atoms, resname="TMP", resseq=1, chain="A"):
    with open(filename, 'w') as f:
        for i, atom in enumerate(atoms):
            name_fmt = fmt_atom_name(atom['name'])[:4]
            f.write(
                f"HETATM"
                f"{i+1:>5d} "
                f"{name_fmt}"
                f" "
                f"{resname:<3}"
                f" "
                f"{chain}"
                f"{resseq:>4d}    "
                f"{atom['x']:8.3f}{atom['y']:8.3f}{atom['z']:8.3f}"
                f"  1.00  0.00          "
                f"{atom['element']:>2}\n"
            )
        f.write("END\n")


def get_residue_prefix(template_index, total_templates):
    """L (1 type) | LA/LB/LC/LD (<=4 types) | L1/L2/... (>=5 types)."""
    if total_templates == 1:
        return "L"
    if total_templates <= 4:
        return f"L{'ABCD'[template_index]}"
    return f"L{template_index + 1}"


# --------------- Main pipeline ---------------

def main():
    parser = argparse.ArgumentParser(
        description="PDB4MUNRO - automatic ligand detection from a metal cage.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Inputs:    .pdb or .xyz (xyz converted with OpenBabel)
Method:    fragments are split by 1.90 A bond cutoff (excludes Pd-N etc.)
           and clustered by graph isomorphism into unique ligand types.

Examples:
  python pdb4munro.py cage.pdb
  python pdb4munro.py cage.xyz --output bone.pdb
        """,
    )
    parser.add_argument('input', help='Input system file (.pdb or .xyz)')
    parser.add_argument('-o', '--output', default='bone.pdb',
                        help='Output PDB filename (default: bone.pdb)')
    parser.add_argument('--naming', choices=['auto', 'sequential'],
                        default='auto',
                        help='auto: LA/LB/LC/LD per template type; '
                             'sequential: L1, L2, L3, ...')
    parser.add_argument('--debug', action='store_true', help='Verbose output')
    args = parser.parse_args()

    print("=" * 70)
    print("PDB4MUNRO - Auto Ligand Detection")
    print("=" * 70)

    # 1. Load atoms
    print(f"\n1. Loading: {args.input}")
    sys_atoms = load_input(args.input)
    if not sys_atoms:
        print("ERROR: no atoms parsed from input")
        sys.exit(1)
    print(f"   Read {len(sys_atoms)} atoms.")

    # 2. Separate metals and ligand atoms
    metals, ligand_pool = [], []
    for a in sys_atoms:
        a['element'] = a['element'].upper()
        (metals if a['element'] in METAL_ELEMENTS else ligand_pool).append(a)
    metal_kinds = sorted({m['element'] for m in metals})
    print(f"   {len(metals)} metal atoms ({metal_kinds}), "
          f"{len(ligand_pool)} ligand atoms.")

    # 3. Find connected ligand fragments
    print("\n2. Identifying ligand fragments by connectivity")
    lig_adj = build_adjacency(ligand_pool)
    components = find_connected_components(ligand_pool, lig_adj)
    ligand_fragments = []
    for comp in components:
        atoms_list = [ligand_pool[i] for i in comp]
        local_adj = build_adjacency(atoms_list)
        ligand_fragments.append((atoms_list, local_adj))
    print(f"   Found {len(ligand_fragments)} discrete ligand fragments.")

    size_hist = defaultdict(int)
    for atoms, _ in ligand_fragments:
        size_hist[len(atoms)] += 1
    for sz, count in sorted(size_hist.items()):
        print(f"     {count} fragment(s) with {sz} atoms")

    # 4. Cluster fragments by isomorphism
    print("\n3. Clustering fragments by graph isomorphism")
    templates, assignments = discover_templates(ligand_fragments)
    print(f"   Found {len(templates)} unique ligand type(s).")
    for tpl in templates:
        prefix = get_residue_prefix(tpl['index'], len(templates))
        print(f"     {tpl['name']} -> prefix {prefix}: "
              f"{tpl['num_atoms']} atoms (seed = ligand #{tpl['src_ligand']+1})")

    # 5. Save each template
    for tpl in templates:
        prefix = get_residue_prefix(tpl['index'], len(templates))
        out = f"{prefix}{tpl['name']}_template.pdb"
        write_template_pdb(out, tpl['atoms'], resname="TMP", resseq=1)
        print(f"     -> wrote {out}")

    # 6. Build bone.pdb
    print(f"\n4. Writing {args.output}")
    output_lines = []
    template_residue_counter = defaultdict(int)

    for i, ((atoms, _), (tpl_idx, mapping)) in enumerate(zip(ligand_fragments, assignments)):
        tpl = templates[tpl_idx]
        template_residue_counter[tpl_idx] += 1

        if args.naming == 'sequential':
            res_name = f"L{i + 1}"
        else:
            prefix = get_residue_prefix(tpl_idx, len(templates))
            res_name = f"{prefix}{template_residue_counter[tpl_idx]}"
        res_seq = i + 1

        for t_idx in range(len(tpl['atoms'])):
            local_idx = mapping[t_idx]
            messy = atoms[local_idx]
            atom_name = tpl['atoms'][t_idx]['name']
            el = tpl['atoms'][t_idx]['element']
            name_fmt = fmt_atom_name(atom_name)[:4]
            output_lines.append(
                f"HETATM"
                f"{t_idx + 1:>5d} "
                f"{name_fmt}"
                f" "
                f"{res_name:<3}"
                f" "
                f"A"
                f"{res_seq:>4d}    "
                f"{messy['x']:8.3f}{messy['y']:8.3f}{messy['z']:8.3f}"
                f"  1.00  0.00          "
                f"{el:>2}\n"
            )
        output_lines.append("TER\n")

    metal_start_res = len(ligand_fragments) + 1
    for i, m in enumerate(metals):
        el = m['element']
        res_seq = metal_start_res + i
        res_name = f"P{i + 1}"
        name_fmt = fmt_atom_name(el)[:4]
        output_lines.append(
            f"HETATM"
            f"{i + 1:>5d} "
            f"{name_fmt}"
            f" "
            f"{res_name:<3}"
            f" "
            f"A"
            f"{res_seq:>4d}    "
            f"{m['x']:8.3f}{m['y']:8.3f}{m['z']:8.3f}"
            f"  1.00  0.00          "
            f"{el:>2}\n"
        )
        output_lines.append("TER\n")

    output_lines.append("END\n")
    with open(args.output, 'w') as f:
        f.writelines(output_lines)

    # 7. Summary
    print("\n" + "=" * 70)
    print("Summary")
    print("=" * 70)
    print(f"  Total ligand fragments:    {len(ligand_fragments)}")
    print(f"  Distinct ligand types:     {len(templates)}")
    print(f"  Metal atoms:               {len(metals)}")
    for tpl in templates:
        prefix = get_residue_prefix(tpl['index'], len(templates))
        count = template_residue_counter[tpl['index']]
        print(f"    {tpl['name']} ({tpl['num_atoms']} atoms): "
              f"{count} occurrence(s) -> {prefix}1..{prefix}{count}")
    print(f"  Output: {args.output}")
    print("=" * 70)


if __name__ == "__main__":
    main()
