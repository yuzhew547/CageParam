#!/usr/bin/env python3
"""
respfit.py - RESP charge fitting for metal-organic cages.

Two modes, auto-selected from the arguments:

  WHOLE-CAGE (default, backward-compatible):
      python respfit.py <molden> [-o output.chg]
    Runs Multiwfn 2-stage RESP on a single molden and emits a CHG file
    consumable by `pdb4munro.py --chg`. Use for M2L4 (and M6L12 if the
    QM converged on the whole cage).

  DIFFERENTIAL (M12L24 fallback):
      python respfit.py --free <src> --cluster <src>
                        --ligand-output L_corr.mol2
                        --metal-output  P_corr.mol2
    Computes q_cage[i] = q_bound[i] + (q_bound[partner(i)] - q_free[i])
    where partner(i) is the ligand's internal mirror image, q_bound is
    averaged over the ligand instances in the Pd(II)L4 cluster, and q_free
    comes from the free-L source. --free / --cluster auto-dispatch on
    extension: .mol2 reads column-9 charges, .chg reads Multiwfn or legacy
    single-column, .molden runs whole-cage RESP first.

Pre-flight classification:
  Both modes inspect every molden input by counting metals and organic
  fragments and report a classification (free_ligand / cluster / cage).
  Use --classify-only to print the report without running RESP.

Auto-classification rules:
  * 0 metals                                -> free_ligand
  * 1 metal,  >=2 organic fragments         -> cluster   (suggest --cluster)
  * >=2 metals AND >=2 organic fragments    -> cage      (whole-cage RESP path)
"""
import argparse
import math
import os
import shutil
import subprocess
import sys
from collections import defaultdict

import numpy as np

try:
    from . import mol2gen_helper as mgh          # package mode
except ImportError:
    import mol2gen_helper as mgh                  # script-mode fallback


DEFAULT_MULTIWFN = os.environ.get(
    "CAGEPIPE_MULTIWFN",
    "/home/gridsan/ywang6/sft/Multiwfn_3.8_bin_Linux_noGUI/Multiwfn_noGUI",
)

# Atomic numbers for elements we may encounter.
TRUE_Z = {
    'H': 1, 'HE': 2, 'LI': 3, 'BE': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8,
    'F': 9, 'NE': 10, 'NA': 11, 'MG': 12, 'AL': 13, 'SI': 14, 'P': 15,
    'S': 16, 'CL': 17, 'AR': 18, 'K': 19, 'CA': 20, 'SC': 21, 'TI': 22,
    'V': 23, 'CR': 24, 'MN': 25, 'FE': 26, 'CO': 27, 'NI': 28, 'CU': 29,
    'ZN': 30, 'GA': 31, 'GE': 32, 'AS': 33, 'SE': 34, 'BR': 35, 'KR': 36,
    'RB': 37, 'SR': 38, 'Y': 39, 'ZR': 40, 'NB': 41, 'MO': 42, 'TC': 43,
    'RU': 44, 'RH': 45, 'PD': 46, 'AG': 47, 'CD': 48, 'IN': 49, 'SN': 50,
    'SB': 51, 'TE': 52, 'I': 53, 'XE': 54, 'CS': 55, 'BA': 56, 'LA': 57,
    'HF': 72, 'TA': 73, 'W': 74, 'RE': 75, 'OS': 76, 'IR': 77, 'PT': 78,
    'AU': 79, 'HG': 80,
}

# Core electrons replaced by the ECP. LANL2DZ/LANL2TZ-style defaults for the
# common cage metals; expand as new ligands appear.
ECP_CORES = {
    'PD': 28, 'PT': 60, 'AU': 60, 'AG': 28,
    'RU': 28, 'RH': 28, 'IR': 60, 'OS': 60,
}

# MK ESP-fitting radius (Angstrom) supplied to Multiwfn's RESP prompt for
# elements not in its built-in table. Values are from AMBER's
# dat/antechamber/ESPPARM.DAT, "MK" column (~ Bondi vdW * 0.95).
# Note: this is NOT the LJ Rmin/2 used in frcmod nonbonded sections (those
# come from Li/Merz ion sets and are ~1.21-1.31 for Pd/Pt — see
# frcmod.ions234lm_126_*).
MK_RADIUS = {
    'PD': 1.55,   # AMBER ESPPARM.DAT
    'PT': 1.66,   # AMBER ESPPARM.DAT
    'AG': 1.63,   # AMBER ESPPARM.DAT
    'AU': 1.58,   # AMBER ESPPARM.DAT
    # Ru/Rh/Os/Ir are 0.00 in ESPPARM.DAT (no AMBER value); estimates below.
    'RU': 1.30, 'RH': 1.30, 'IR': 1.30, 'OS': 1.30,
}

METALS = mgh.METALS


# ============================================================
# Multiwfn driving (whole-cage RESP from a molden)
# ============================================================

def parse_molden(path):
    """Return (header_lines, atom_records).

    atom_records is a list of dicts with keys: line_idx, symbol, index, z_in,
    x, y, z, raw. Coordinate units are not parsed (we never modify them).
    """
    with open(path) as f:
        lines = f.readlines()

    atoms = []
    in_atoms = False
    for i, line in enumerate(lines):
        s = line.strip()
        if s.lower().startswith('[atoms]'):
            in_atoms = True
            continue
        if in_atoms:
            if s.startswith('['):
                in_atoms = False
                continue
            if not s:
                continue
            parts = line.split()
            if len(parts) < 6:
                continue
            try:
                atoms.append({
                    'line_idx': i,
                    'symbol': parts[0],
                    'index': int(parts[1]),
                    'z_in': int(parts[2]),
                    'x': float(parts[3]),
                    'y': float(parts[4]),
                    'z': float(parts[5]),
                    'raw': line,
                })
            except ValueError:
                continue
    return lines, atoms


def patch_molden(src, dst):
    """Write `dst` with ECP-corrected nuclear charges. Returns count patched."""
    lines, atoms = parse_molden(src)
    if not atoms:
        raise RuntimeError(f"No atoms parsed from [Atoms] section of {src}")

    patched = 0
    overrides = {}
    for a in atoms:
        el = a['symbol'].upper()
        true_z = TRUE_Z.get(el)
        core = ECP_CORES.get(el)
        if true_z is None or core is None:
            continue
        eff_z = true_z - core
        if a['z_in'] == true_z:
            overrides[a['line_idx']] = eff_z
            patched += 1
        elif a['z_in'] != eff_z and a['z_in'] != true_z:
            print(f"  ! atom {a['index']} ({el}): unexpected Z={a['z_in']} "
                  f"(true={true_z}, ecp_eff={eff_z}); leaving as-is")

    if patched == 0:
        if src != dst:
            shutil.copy(src, dst)
        return 0

    with open(dst, 'w') as f:
        for i, line in enumerate(lines):
            if i in overrides:
                parts = line.split()
                parts[2] = str(overrides[i])
                f.write(f"{parts[0]:<3} {int(parts[1]):>3} {parts[2]:>4} "
                        f"{float(parts[3]):>12.5f}{float(parts[4]):>12.5f}"
                        f"{float(parts[5]):>12.5f}\n")
            else:
                f.write(line)
    return patched


def build_multiwfn_input(elements_present):
    """Multiwfn stdin sequence: main menu -> RESP -> radii -> save -> quit.

    Multiwfn prompts for the MK fitting radius once per element it does not
    know; the order tracks the first appearance in the molecule, which (for
    a well-ordered TeraChem molden) matches `elements_present`."""
    seq = ['7', '18', '1']
    seen = set()
    for el in elements_present:
        e = el.upper()
        if e in MK_RADIUS and e not in seen:
            seq.append(f"{MK_RADIUS[e]}")
            seen.add(e)
    seq.append('y')
    seq.append('q')
    return '\n'.join(seq) + '\n'


def write_local_settings(nthreads):
    """Write a minimal settings.ini in cwd to override Multiwfn's nthreads."""
    with open('settings.ini', 'w') as f:
        f.write("// Auto-generated by respfit.py\n")
        f.write(f"  nthreads= {int(nthreads)}\n")


def detect_nthreads(requested):
    if requested is not None:
        return int(requested)
    slurm = os.environ.get('SLURM_CPUS_PER_TASK')
    if slurm and slurm.isdigit():
        return int(slurm)
    cpu = os.cpu_count()
    return cpu if cpu else 4


def run_multiwfn(multiwfn_bin, molden_path, input_path, log_path):
    """Stream `input_path` into Multiwfn; capture stdout/stderr to `log_path`."""
    with open(input_path) as fin, open(log_path, 'w') as fout:
        rc = subprocess.run(
            [multiwfn_bin, molden_path, '-silent'],
            stdin=fin, stdout=fout, stderr=subprocess.STDOUT,
        ).returncode
    return rc


def run_resp_on_molden(molden_path, output_chg, multiwfn=DEFAULT_MULTIWFN,
                       nthreads=None, keep_aux=False):
    """End-to-end RESP from a molden: ECP-patch -> Multiwfn -> CHG.

    Returns the path to the produced CHG file (== output_chg)."""
    if not os.path.isfile(molden_path):
        raise RuntimeError(f"molden not found: {molden_path}")
    if not os.path.isfile(multiwfn):
        raise RuntimeError(f"Multiwfn binary not found: {multiwfn}")

    print(f"\n[RESP] Reading molden: {molden_path}")
    _, atoms = parse_molden(molden_path)
    if not atoms:
        raise RuntimeError("no atoms parsed from [Atoms] section")
    print(f"       {len(atoms)} atoms")

    print("[RESP] Checking ECP nuclear charges")
    work_molden = "_respfit_work.molden"
    n_patched = patch_molden(molden_path, work_molden)
    if n_patched:
        ecp_seen = sorted({a['symbol'].upper() for a in atoms
                           if a['symbol'].upper() in ECP_CORES
                           and a['z_in'] == TRUE_Z[a['symbol'].upper()]})
        print(f"       Patched {n_patched} atom(s) for ECP: {ecp_seen}")
    else:
        print("       No ECP patch needed (column-3 Z already effective)")

    elements_in_order = []
    seen_el = set()
    _, work_atoms = parse_molden(work_molden)
    for a in work_atoms:
        e = a['symbol'].upper()
        if e not in seen_el:
            seen_el.add(e)
            elements_in_order.append(e)

    print("[RESP] Generating Multiwfn input")
    resp_in = build_multiwfn_input(elements_in_order)
    in_path = "_respfit.in"
    with open(in_path, 'w') as f:
        f.write(resp_in)

    nthr = detect_nthreads(nthreads)
    write_local_settings(nthr)
    print(f"[RESP] Running Multiwfn (nthreads={nthr}, this may take a while)")
    log_path = "_respfit.log"
    rc = run_multiwfn(multiwfn, work_molden, in_path, log_path)
    print(f"       rc={rc}, log -> {log_path}")
    if rc != 0:
        print("       Multiwfn returned non-zero; inspect the log.")

    base = os.path.splitext(work_molden)[0]
    candidate = f"{base}.chg"
    if not os.path.isfile(candidate):
        candidate = None
        for f in os.listdir('.'):
            if f.lower().endswith('.chg'):
                candidate = f
                break
    if not candidate or not os.path.isfile(candidate):
        raise RuntimeError(f"no .chg file produced. Inspect {log_path}.")

    if os.path.abspath(candidate) != os.path.abspath(output_chg):
        shutil.copy(candidate, output_chg)
    print(f"[RESP] -> {output_chg}")

    if not keep_aux:
        for f in (work_molden, in_path, log_path, candidate, 'settings.ini'):
            if f and f != output_chg and os.path.isfile(f):
                try:
                    os.remove(f)
                except OSError:
                    pass

    return output_chg


# ============================================================
# Charge / geometry parsers (shared by both modes)
# ============================================================

def parse_pdb_for_chg(pdb_path):
    """Return [(element, x, y, z, atom_name)] in PDB atom order."""
    out = []
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            try:
                x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
            except ValueError:
                continue
            name = line[12:16].strip()
            el = line[76:78].strip().upper()
            if not el:
                el = "".join(c for c in name if c.isalpha()).upper()
                if el[:2] in ("BR", "CL", "PD", "PT", "AU", "AG"):
                    el = el[:2]
                else:
                    el = el[:1]
            out.append((el, x, y, z, name))
    return out


def read_chg(chg_path, pdb_path=None):
    """Return list of dicts {element, x, y, z, charge, name}.

    Multiwfn format (5+ tokens, first non-numeric): ele x y z q.
    Legacy single-column: charges only; needs pdb_path for atom order."""
    with open(chg_path) as f:
        first = ""
        for line in f:
            s = line.strip()
            if s:
                first = s
                break
    parts_first = first.split()
    multiwfn = False
    if len(parts_first) >= 5:
        try:
            float(parts_first[0])
        except ValueError:
            multiwfn = True

    if multiwfn:
        out = []
        with open(chg_path) as f:
            for line in f:
                p = line.split()
                if len(p) < 5:
                    continue
                try:
                    out.append({
                        'element': p[0].upper(),
                        'x': float(p[1]), 'y': float(p[2]), 'z': float(p[3]),
                        'charge': float(p[4]),
                        'name': None,
                    })
                except ValueError:
                    continue
        return out

    if pdb_path is None:
        raise RuntimeError(
            f"{chg_path} looks like single-column charges; supply the matching "
            f"--*-pdb so atom order/elements/coords are known."
        )
    charges = []
    with open(chg_path) as f:
        for line in f:
            for tok in line.split():
                try:
                    charges.append(float(tok))
                except ValueError:
                    pass
    pdb_atoms = parse_pdb_for_chg(pdb_path)
    if len(charges) != len(pdb_atoms):
        raise RuntimeError(
            f"{chg_path} has {len(charges)} charges but {pdb_path} has "
            f"{len(pdb_atoms)} atoms"
        )
    out = []
    for (el, x, y, z, name), q in zip(pdb_atoms, charges):
        out.append({'element': el, 'x': x, 'y': y, 'z': z,
                    'charge': q, 'name': name})
    return out


def read_mol2_full(mol2_path):
    """Parse a full mol2: atoms (with coords/charges/types) and bonds.

    Returns (atoms, bonds, atom_types).
        atoms: list of dicts {name, element, x, y, z, charge, sybyl, file_id}
        bonds: list of (a_idx_0based, b_idx_0based, btype)
        atom_types: {atom_name_upper: sybyl_type}
    """
    atoms = []
    bonds = []
    types = {}
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
                el = "".join(c for c in name if c.isalpha()).upper()
                if el[:2] in ("BR", "CL", "PD", "PT", "AU", "AG"):
                    el = el[:2]
                else:
                    el = el[:1]
                idx = len(atoms)
                id_map[file_id] = idx
                atoms.append({
                    'name': name, 'element': el, 'x': x, 'y': y, 'z': z,
                    'charge': charge, 'sybyl': sybyl, 'file_id': file_id,
                })
                types[name.strip().upper()] = sybyl
            elif section == "@<TRIPOS>BOND" and len(p) >= 4:
                try:
                    a = id_map[p[1]]; b = id_map[p[2]]
                except KeyError:
                    continue
                bonds.append((a, b, p[3]))
    return atoms, bonds, types


def to_helper_atoms(records, resname="TMP"):
    """Convert {element,x,y,z,name?} records to mgh.Atom objects."""
    out = []
    for r in records:
        nm = r.get('name') or r['element']
        a = mgh.Atom(nm, r['x'], r['y'], r['z'], resname, 1, element=r['element'])
        a.charge = r.get('charge', 0.0)
        out.append(a)
    return out


def cluster_records(records):
    """Group connected organic atoms via covalent-radii adjacency.

    Returns (organic_groups, metal_records). Each organic_group is a list of
    record dicts (one connected fragment); metals are returned individually."""
    organics = [r for r in records if r['element'] not in METALS]
    metals = [r for r in records if r['element'] in METALS]

    helper_atoms = to_helper_atoms(organics)
    adj = mgh.build_adjacency(helper_atoms, use_covalent_radii=True)
    n = len(organics)
    visited = [False] * n
    groups = []
    for i in range(n):
        if visited[i]:
            continue
        stack = [i]
        comp = []
        while stack:
            j = stack.pop()
            if visited[j]:
                continue
            visited[j] = True
            comp.append(organics[j])
            for k in adj[j]:
                if not visited[k]:
                    stack.append(k)
        groups.append(comp)
    return groups, metals


# ============================================================
# Auto-classification
# ============================================================

def classify_records(records):
    """Return a dict describing what `records` looks like:

      {'n_atoms': int,
       'n_metals': int,
       'metal_elements': sorted list,
       'n_organic_fragments': int,
       'organic_fragment_sizes': sorted list,
       'classification': 'free_ligand' | 'cluster' | 'cage' | 'unknown'}
    """
    organic_groups, metals = cluster_records(records)
    sizes = sorted([len(g) for g in organic_groups])
    metal_els = sorted({m['element'] for m in metals})

    if len(metals) == 0:
        cls = "free_ligand"
    elif len(metals) >= 2 and len(organic_groups) >= 2:
        cls = "cage"
    elif len(metals) == 1 and len(organic_groups) >= 2:
        cls = "cluster"
    elif len(metals) >= 1 and len(organic_groups) == 1:
        # Could be: a chelated single ligand (metal embedded) or a tightly-packed
        # fragment we failed to split. Treat as cluster-like but flag.
        cls = "cluster"
    else:
        cls = "unknown"

    return {
        'n_atoms': len(records),
        'n_metals': len(metals),
        'metal_elements': metal_els,
        'n_organic_fragments': len(organic_groups),
        'organic_fragment_sizes': sizes,
        'classification': cls,
    }


def classify_source(spec, pdb=None):
    """Classify any source file (.mol2 / .chg / .molden / .pdb).

    For .molden, this peeks at [Atoms] without running RESP — cheap and
    always safe to call upfront."""
    ext = os.path.splitext(spec)[1].lower()
    records = []
    if ext == ".mol2":
        atoms, _, _ = read_mol2_full(spec)
        records = [{'element': a['element'], 'x': a['x'], 'y': a['y'], 'z': a['z'],
                    'charge': a['charge'], 'name': a['name']} for a in atoms]
    elif ext == ".chg":
        records = read_chg(spec, pdb_path=pdb)
    elif ext == ".molden":
        _, atoms = parse_molden(spec)
        for a in atoms:
            records.append({'element': a['symbol'].upper(),
                            'x': a['x'], 'y': a['y'], 'z': a['z'],
                            'charge': 0.0, 'name': None})
    elif ext == ".pdb":
        for el, x, y, z, name in parse_pdb_for_chg(spec):
            records.append({'element': el, 'x': x, 'y': y, 'z': z,
                            'charge': 0.0, 'name': name})
    else:
        raise RuntimeError(f"Unsupported extension {ext!r} for {spec}")
    return classify_records(records)


def print_classification(label, info):
    cls = info['classification']
    metals = info['metal_elements']
    sizes = info['organic_fragment_sizes']
    print(f"   {label}: {info['n_atoms']} atoms, "
          f"{info['n_metals']} metal(s){' ' + ','.join(metals) if metals else ''}, "
          f"{info['n_organic_fragments']} organic fragment(s)"
          f"{' ' + str(sizes) if sizes else ''}")
    print(f"   {label}: classified as '{cls}'")


# ============================================================
# Differential mode helpers
# ============================================================

def load_charge_source(spec, pdb=None, tag="src", multiwfn=DEFAULT_MULTIWFN,
                       nthreads=None):
    """Resolve --free / --cluster: dispatch on extension and return either
        ('mol2', records, bonds, types)        if a mol2
        ('records', records, None, None)        for chg/molden sources
    """
    ext = os.path.splitext(spec)[1].lower()
    if ext == ".mol2":
        atoms, bonds, types = read_mol2_full(spec)
        recs = [{'element': a['element'], 'x': a['x'], 'y': a['y'], 'z': a['z'],
                 'charge': a['charge'], 'name': a['name']} for a in atoms]
        return 'mol2', recs, bonds, types
    if ext == ".chg":
        return 'records', read_chg(spec, pdb_path=pdb), None, None
    if ext == ".molden":
        out_chg = f"_respfit_{tag}.chg"
        run_resp_on_molden(spec, out_chg, multiwfn=multiwfn, nthreads=nthreads,
                           keep_aux=False)
        return 'records', read_chg(out_chg, pdb_path=None), None, None
    raise RuntimeError(f"Unsupported extension {ext!r} for {spec}")


def find_self_automorphism(tpl_atoms, tpl_adj, max_autos=64):
    """Return a non-identity automorphism mapping tpl_idx -> tpl_idx as a
    dict, or None if the template has only the identity automorphism.

    Enumerates up to `max_autos` automorphisms; picks the one that swaps the
    most atoms (the dominant symmetry — usually the C2 mirror of a ditopic
    ligand)."""
    t_sigs = [mgh.get_extended_signature(i, tpl_atoms, tpl_adj)
              for i in range(len(tpl_atoms))]
    n = len(tpl_atoms)

    sig_counts = defaultdict(int)
    for s in t_sigs:
        sig_counts[s] += 1
    order = sorted(range(n), key=lambda i: sig_counts[t_sigs[i]])

    found = []
    mapping = {}
    used = set()

    def backtrack(pos):
        if len(found) >= max_autos:
            return
        if pos == n:
            if not all(mapping[k] == k for k in mapping):
                found.append(dict(mapping))
            return
        t_idx = order[pos]
        my_sig = t_sigs[t_idx]
        for cand in range(n):
            if cand in used or t_sigs[cand] != my_sig:
                continue
            ok = True
            for nbr in tpl_adj[t_idx]:
                if nbr in mapping and mapping[nbr] not in tpl_adj[cand]:
                    ok = False
                    break
            if ok:
                mapping[t_idx] = cand
                used.add(cand)
                backtrack(pos + 1)
                if len(found) >= max_autos:
                    return
                del mapping[t_idx]
                used.remove(cand)

    backtrack(0)
    if not found:
        return None
    found.sort(key=lambda m: -sum(1 for k, v in m.items() if k != v))
    return found[0]


def map_records_to_template(records, tpl_atoms, tpl_adj, label=""):
    """Match a set of atom records to the template by isomorphism."""
    if len(records) != len(tpl_atoms):
        print(f"   [{label}] size mismatch: {len(records)} records vs "
              f"{len(tpl_atoms)} template atoms")
        return None
    tgt_atoms = to_helper_atoms(records)
    tgt_adj = mgh.build_adjacency(tgt_atoms, use_covalent_radii=True)
    return mgh.solve_isomorphism(tpl_atoms, tpl_adj, tgt_atoms, tgt_adj)


def write_diff_mol2(path, resname, atoms, bonds, header_note="RESP-corrected"):
    """Write a mol2 in the same shape as mol2gen_helper.write_mol2."""
    with open(path, 'w') as f:
        f.write(f"@<TRIPOS>MOLECULE\n{resname}\n")
        f.write(f"{len(atoms)} {len(bonds)} 1 0 0\n")
        f.write(f"SMALL\nUSER_CHARGES\n{header_note}\n\n")
        f.write("@<TRIPOS>ATOM\n")
        for i, a in enumerate(atoms, 1):
            f.write(f"{i:>7} {a['name']:<8} {a['x']:>10.4f} {a['y']:>10.4f} "
                    f"{a['z']:>10.4f} {a['sybyl']:<5} {1:>6} {resname:<5} "
                    f"{a['charge']:>10.6f}\n")
        f.write("@<TRIPOS>BOND\n")
        for i, (a, b, bt) in enumerate(bonds, 1):
            f.write(f"{i:>6} {a+1:>5} {b+1:>5} {bt:>4}\n")
        f.write("@<TRIPOS>SUBSTRUCTURE\n")
        f.write(f"   1 {resname:<9} 1 TEMP              0 **** **** 0 ROOT\n")


# ============================================================
# Mode entry points
# ============================================================

def run_whole_cage(args):
    """Whole-cage mode: single molden -> single CHG."""
    print("=" * 70)
    print(f"RESPFIT (whole-cage) - {args.molden}")
    print("=" * 70)

    print("\n1. Pre-flight classification")
    info = classify_source(args.molden)
    print_classification(os.path.basename(args.molden), info)
    if info['classification'] == 'free_ligand':
        print("   NOTE: 0 metals detected. If you intended differential mode, "
              "re-run with `--free <this> --cluster <pdl4.molden>`.")
    elif info['classification'] == 'cluster':
        print("   NOTE: this looks like a small cluster (1 metal + multiple "
              "ligand fragments). For M12L24 workflows use:")
        print("     respfit.py --free <free.mol2_or_molden> "
              f"--cluster {args.molden} --ligand-output L_corr.mol2 "
              "--metal-output P_corr.mol2")
    elif info['classification'] == 'cage':
        print("   OK: looks like a whole cage; whole-cage RESP is the right "
              "path.")
    if args.classify_only:
        print("\n--classify-only set; not running Multiwfn.")
        return 0

    print("\n2. Running RESP on molden")
    run_resp_on_molden(args.molden, args.output, multiwfn=args.multiwfn,
                       nthreads=args.nthreads, keep_aux=args.keep_aux)

    print("\n" + "=" * 70)
    print(f"Done. CHG: {args.output}")
    print("=" * 70)
    return 0


def run_differential(args):
    """Differential mode: free + cluster -> corrected ligand & metal mol2."""
    print("=" * 70)
    print("RESPFIT (differential) - "
          "q_cage = q_bound + (q_bound[mirror] - q_free)")
    print("=" * 70)

    # Pre-flight classification of inputs
    print("\n1. Pre-flight classification")
    free_info = classify_source(args.free, pdb=args.free_pdb)
    print_classification("free   ", free_info)
    cluster_info = classify_source(args.cluster, pdb=args.cluster_pdb)
    print_classification("cluster", cluster_info)
    if free_info['classification'] != 'free_ligand':
        print(f"   WARNING: --free has {free_info['n_metals']} metal(s); "
              f"expected 0. Free-L charges should come from a metal-free source.")
    if cluster_info['classification'] not in ('cluster', 'cage'):
        print(f"   WARNING: --cluster has {cluster_info['n_metals']} metal(s); "
              f"expected >=1.")
    if cluster_info['classification'] == 'cage':
        print("   NOTE: --cluster looks like a whole cage. Differential mode "
              "still works, but if QM converged on the whole cage you can "
              "use whole-cage mode directly:")
        print(f"     respfit.py {args.cluster}")
    if args.classify_only:
        print("\n--classify-only set; not running RESP.")
        return 0

    # 1. Load free template (atom names, sybyl types, bonds)
    tpl_path = args.free_template or (args.free if args.free.endswith(".mol2")
                                      else None)
    if tpl_path is None:
        print("ERROR: --free-template is required when --free is not a .mol2")
        return 1
    if not os.path.isfile(tpl_path):
        print(f"ERROR: free template not found: {tpl_path}")
        return 1
    print(f"\n2. Loading free-L template (topology): {tpl_path}")
    tpl_atoms_raw, tpl_bonds, tpl_types = read_mol2_full(tpl_path)
    tpl_helper_atoms = to_helper_atoms(tpl_atoms_raw)
    tpl_adj = mgh.build_adjacency(tpl_helper_atoms, use_covalent_radii=True)
    n_tpl = len(tpl_helper_atoms)
    print(f"   Template: {n_tpl} atoms, {len(tpl_bonds)} bonds")

    # 2. Free-L charges
    print(f"\n3. Loading free-L charges: {args.free}")
    _, free_recs, _, _ = load_charge_source(args.free, pdb=args.free_pdb,
                                             tag="free", multiwfn=args.multiwfn,
                                             nthreads=args.nthreads)
    if len(free_recs) != n_tpl:
        groups, _ = cluster_records(free_recs)
        groups.sort(key=len, reverse=True)
        if not groups or len(groups[0]) != n_tpl:
            print(f"ERROR: free input has {len(free_recs)} atoms; largest "
                  f"organic fragment has {(len(groups[0]) if groups else 0)} "
                  f"(template wants {n_tpl}).")
            return 1
        free_recs = groups[0]
        print(f"   Picked largest organic fragment ({n_tpl} atoms) from free input")
    map_free = map_records_to_template(free_recs, tpl_helper_atoms, tpl_adj, "free")
    if map_free is None:
        print("ERROR: free-L did not match the template by graph isomorphism")
        return 1
    q_free = {t: free_recs[map_free[t]]['charge'] for t in range(n_tpl)}
    print(f"   q_free populated for {n_tpl} atoms; sum = {sum(q_free.values()):+.4f}")

    # 3. Cluster
    print(f"\n4. Loading cluster: {args.cluster}")
    _, cluster_recs, _, _ = load_charge_source(args.cluster, pdb=args.cluster_pdb,
                                                tag="cluster", multiwfn=args.multiwfn,
                                                nthreads=args.nthreads)
    groups, metals = cluster_records(cluster_recs)
    print(f"   {len(cluster_recs)} atoms -> {len(groups)} organic fragment(s) "
          f"+ {len(metals)} metal atom(s)")
    matched_groups = []
    for gi, g in enumerate(groups):
        if len(g) != n_tpl:
            continue
        m = map_records_to_template(g, tpl_helper_atoms, tpl_adj, f"frag{gi+1}")
        if m is None:
            continue
        matched_groups.append((g, m))
    if not matched_groups:
        print("ERROR: no cluster ligand fragment matched the template")
        return 1
    print(f"   Matched {len(matched_groups)} ligand instance(s) in cluster")

    q_bound = {}
    for t in range(n_tpl):
        vals = [g[m[t]]['charge'] for g, m in matched_groups]
        q_bound[t] = float(np.mean(vals))
    print(f"   q_bound populated; sum (per ligand) = {sum(q_bound.values()):+.4f}")

    # 4. Symmetry partner
    print("\n5. Finding ligand internal symmetry (graph automorphism)")
    auto = find_self_automorphism(tpl_helper_atoms, tpl_adj)
    if auto is None:
        print("   WARNING: no non-identity automorphism found; falling back to "
              "identity (q_cage = q_bound; no mirror correction)")
        partner = {i: i for i in range(n_tpl)}
    else:
        n_swap = sum(1 for k, v in auto.items() if k != v)
        print(f"   Found automorphism that swaps {n_swap}/{n_tpl} atoms.")
        partner = auto

    # 5. Apply formula
    print("\n6. Applying q_cage = q_bound[i] + (q_bound[partner(i)] - q_free[i])")
    q_cage = {i: q_bound[i] + (q_bound[partner[i]] - q_free[i])
              for i in range(n_tpl)}
    print(f"   q_cage sum = {sum(q_cage.values()):+.4f}")

    # 6. Write outputs
    print(f"\n7. Writing {args.ligand_output} (resname {args.ligand_resname})")
    out_atoms = []
    for i, a in enumerate(tpl_atoms_raw):
        out_atoms.append({'name': a['name'], 'x': a['x'], 'y': a['y'], 'z': a['z'],
                          'sybyl': a['sybyl'], 'charge': q_cage[i]})
    write_diff_mol2(args.ligand_output, args.ligand_resname, out_atoms, tpl_bonds,
                    header_note=f"RESP differential; net = "
                                f"{sum(q_cage.values()):+.4f}")

    if metals:
        m_q = float(np.mean([m['charge'] for m in metals]))
        m_el = metals[0]['element']
        print(f"\n8. Writing {args.metal_output} (resname {args.metal_resname}, "
              f"{m_el}, q = {m_q:+.4f})")
        m_atom = {'name': m_el, 'x': metals[0]['x'], 'y': metals[0]['y'],
                  'z': metals[0]['z'], 'sybyl': "M0", 'charge': m_q}
        write_diff_mol2(args.metal_output, args.metal_resname, [m_atom], [],
                        header_note="Metal RESP charge from cluster")
    else:
        print("\n8. No metal atom in cluster - skipping --metal-output")

    print("\n" + "=" * 70)
    print("Done.")
    print("=" * 70)
    return 0


# ============================================================
# CLI
# ============================================================

def main():
    ap = argparse.ArgumentParser(
        description="RESP charge fitting for metal-organic cages "
                    "(whole-cage and differential modes auto-selected).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Whole-cage mode (M2L4, small M6L12):
  python respfit.py cage.molden -o cage.chg

Differential mode (M12L24 or any cage too large for whole-cage RESP):
  python respfit.py --free MOL.mol2 --cluster pdl4.molden \\
      --free-template MOL.mol2 \\
      --ligand-output L_corr.mol2 --metal-output P_corr.mol2

Inspect inputs without running RESP:
  python respfit.py cage.molden --classify-only
  python respfit.py --free MOL.mol2 --cluster pdl4.molden --classify-only
        """,
    )

    # Whole-cage positional
    ap.add_argument("molden", nargs="?",
                    help="Whole-cage molden (whole-cage mode). Omit when "
                         "using --free + --cluster.")
    ap.add_argument("-o", "--output", default="input.chg",
                    help="(whole-cage) Output CHG path (default: input.chg)")

    # Differential
    ap.add_argument("--free", default=None,
                    help="(differential) Free ligand source: .mol2 / .chg / .molden")
    ap.add_argument("--cluster", default=None,
                    help="(differential) Pd(II)L4 cluster source: .mol2 / .chg / .molden")
    ap.add_argument("--free-pdb", default=None,
                    help="PDB matching --free (only needed for legacy single-column .chg)")
    ap.add_argument("--cluster-pdb", default=None,
                    help="PDB matching --cluster (only needed for legacy single-column .chg)")
    ap.add_argument("--free-template", default=None,
                    help="Antechamber-typed mol2 of the free ligand "
                         "(atom names/types/bonds). Defaults to --free if it's a .mol2.")
    ap.add_argument("--ligand-output", default="L_corr.mol2",
                    help="(differential) Output corrected ligand mol2")
    ap.add_argument("--metal-output", default="P_corr.mol2",
                    help="(differential) Output metal mol2")
    ap.add_argument("--ligand-resname", default="ML1",
                    help="Residue name in --ligand-output (default: ML1)")
    ap.add_argument("--metal-resname", default="P1",
                    help="Residue name in --metal-output (default: P1)")

    # Shared
    ap.add_argument("--multiwfn", default=DEFAULT_MULTIWFN,
                    help="Path to Multiwfn binary.")
    ap.add_argument("--nthreads", default=None,
                    help="OpenMP threads for Multiwfn (default: $SLURM_CPUS_PER_TASK).")
    ap.add_argument("--keep-aux", action="store_true",
                    help="(whole-cage) keep intermediate molden / log / Multiwfn input.")
    ap.add_argument("--classify-only", action="store_true",
                    help="Inspect inputs and report classification, then exit "
                         "without running RESP or fitting charges.")
    ap.add_argument("--debug", action="store_true")

    args = ap.parse_args()

    # Mode dispatch
    diff_flags = (args.free is not None) or (args.cluster is not None)
    whole_flag = args.molden is not None

    if diff_flags and whole_flag:
        ap.error("Provide either a positional <molden> (whole-cage) OR "
                 "--free + --cluster (differential), not both.")

    if diff_flags:
        if args.free is None or args.cluster is None:
            ap.error("Differential mode requires both --free and --cluster.")
        return run_differential(args)

    if whole_flag:
        return run_whole_cage(args)

    ap.error("No input provided. Pass a positional <molden> for whole-cage RESP "
             "or --free + --cluster for differential RESP.")


if __name__ == "__main__":
    sys.exit(main())
