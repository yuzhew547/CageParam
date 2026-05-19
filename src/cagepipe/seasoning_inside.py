#!/usr/bin/env python3
"""filling - parametrize a guest, place it inside the cage cavity, and
optionally surround the cage with N counterions in one pass.

Console-script name: ``filling`` (legacy module name: seasoning_inside).

Strategy (guest, inside):
  1. Align the guest's longest principal axis with the cage's longest axis
     (Pt-Pt vector when present, otherwise PCA on heavy atoms) so an elongated
     guest lies along the natural channel.
  2. Try random rotations + small jitter at the cavity center.
  3. Use a tight clash cutoff (default 2.0 A); relax progressively if no fit.

Strategy (counterions, outside): when --counterions N > 0, place N copies of
a counterion template in random shells around the cage (radius = cage radius +
buffer .. + buffer + 15 A). Residue names cycle BFA, BFB, ..., BFZ, BGA, ...
(matching the standalone ``seasoning`` script). Per-anion ``BFx.{pdb,mol2,
frcmod}`` are written next to the cage by re-using the counterion template;
the guest is treated as a fixed obstacle so counterions never overlap it.

Inputs:
  --cage bone.pdb              cage PDB (no guest)
  --templ GUEST.pdb|.xyz       guest molecule; .xyz auto-converted via OpenBabel
  --resname GS1                guest residue name (default GS1)
  [--autoparam]                also produce GUEST.mol2 + GUEST.frcmod via
                               antechamber (AM1-BCC, GAFF2) + parmchk2
  [--charge N]                 net charge passed to antechamber for the guest
  [--counterions N]            also place N counterions around the cage
  [--counterion-templ BFA.pdb] counterion PDB template (resname cycles BFA..)
  [--counterion-mol2 BFA.mol2] counterion MOL2 to be renamed per copy
  [--counterion-frcmod BFA.frcmod]   counterion FRCMOD (copied per copy)

Output:
  tastybone.pdb            cage + 1 guest (GS1) inside + N counterions outside
  [GUEST.mol2, GUEST.frcmod]   when --autoparam is set
  [BFA.pdb/mol2/frcmod, BFB...]    one set per placed counterion
"""
import argparse, math, os, random, re, shutil, string, subprocess, sys, tempfile
from itertools import combinations

try:
    import numpy as np
except ImportError:
    sys.stderr.write("numpy is required for filling/seasoning_inside.py\n")
    sys.exit(2)


# ---------- External tool discovery (mirrors pdb4munro) ----------
OBABEL_CANDIDATES = [
    "/home/gridsan/ywang6/sft/build/bin/obabel",
    "obabel",
    "babel",
]
ANTECHAMBER_CANDIDATES = [
    "/home/gridsan/ywang6/.conda/envs/AmberTools25/bin/antechamber",
    "antechamber",
]
PARMCHK2_CANDIDATES = [
    "/home/gridsan/ywang6/.conda/envs/AmberTools25/bin/parmchk2",
    "parmchk2",
]


def _find_exe(candidates):
    for p in candidates:
        if os.path.isfile(p) and os.access(p, os.X_OK):
            return p
        which = shutil.which(p)
        if which:
            return which
    return None


def _obabel_env():
    env = os.environ.copy()
    libdir = "/home/gridsan/ywang6/sft/build/lib"
    datadir = "/home/gridsan/ywang6/sft/openbabel-openbabel-2-4-0/data"
    if os.path.isdir(libdir):
        env.setdefault("BABEL_LIBDIR", libdir)
    if os.path.isdir(datadir):
        env.setdefault("BABEL_DATADIR", datadir)
    return env


def xyz_to_pdb(xyz_path, pdb_path):
    obabel = _find_exe(OBABEL_CANDIDATES)
    if obabel is None:
        raise RuntimeError(
            "OpenBabel (obabel) not found. Install it or extend OBABEL_CANDIDATES."
        )
    print(f"  obabel: {xyz_path} -> {pdb_path}")
    res = subprocess.run(
        [obabel, xyz_path, f"-O{pdb_path}"],
        capture_output=True, text=True, env=_obabel_env(),
    )
    if res.returncode != 0 or not os.path.isfile(pdb_path):
        raise RuntimeError(
            f"obabel failed (rc={res.returncode}):\nSTDOUT:\n{res.stdout}\nSTDERR:\n{res.stderr}"
        )


def run_autoparam(pdb_path, stem, charge, resname, atom_type="gaff2"):
    """antechamber + parmchk2 -> <stem>.mol2 and <stem>.frcmod next to pdb_path."""
    antechamber = _find_exe(ANTECHAMBER_CANDIDATES)
    if antechamber is None:
        raise RuntimeError(
            "antechamber not found. Activate AmberTools or extend ANTECHAMBER_CANDIDATES."
        )
    parmchk2 = _find_exe(PARMCHK2_CANDIDATES)
    if parmchk2 is None:
        raise RuntimeError(
            "parmchk2 not found. Activate AmberTools or extend PARMCHK2_CANDIDATES."
        )

    out_dir = os.path.dirname(os.path.abspath(pdb_path)) or "."
    final_mol2   = os.path.join(out_dir, f"{stem}.mol2")
    final_frcmod = os.path.join(out_dir, f"{stem}.frcmod")

    # Run antechamber + parmchk2 in a scratch tempdir so sqm.*, ATOMTYPE.INF,
    # ANTECHAMBER_*.AC, NEWPDB.PDB etc. never land in the cage directory.
    with tempfile.TemporaryDirectory(prefix=f"filling_{stem}_") as scratch:
        scratch_mol2   = os.path.join(scratch, f"{stem}.mol2")
        scratch_frcmod = os.path.join(scratch, f"{stem}.frcmod")
        scratch_pdb    = os.path.join(scratch, os.path.basename(pdb_path))
        shutil.copy2(pdb_path, scratch_pdb)

        print(f"  antechamber: AM1-BCC, GAFF2, nc={charge}, resname={resname}")
        ac = subprocess.run(
            [antechamber,
             "-i", scratch_pdb,    "-fi", "pdb",
             "-o", scratch_mol2,   "-fo", "mol2",
             "-c", "bcc", "-nc", str(int(charge)),
             "-at", atom_type, "-rn", resname,
             "-pf", "y"],
            capture_output=True, text=True, cwd=scratch,
        )
        if ac.returncode != 0 or not os.path.isfile(scratch_mol2):
            raise RuntimeError(
                f"antechamber failed (rc={ac.returncode}):\nSTDOUT:\n{ac.stdout}\nSTDERR:\n{ac.stderr}"
            )

        print(f"  parmchk2: GAFF2 frcmod")
        pc = subprocess.run(
            [parmchk2, "-i", scratch_mol2, "-f", "mol2",
             "-o", scratch_frcmod, "-s", atom_type],
            capture_output=True, text=True, cwd=scratch,
        )
        if pc.returncode != 0 or not os.path.isfile(scratch_frcmod):
            raise RuntimeError(
                f"parmchk2 failed (rc={pc.returncode}):\nSTDOUT:\n{pc.stdout}\nSTDERR:\n{pc.stderr}"
            )

        shutil.copy2(scratch_mol2,   final_mol2)
        shutil.copy2(scratch_frcmod, final_frcmod)

    return final_mol2, final_frcmod


def parse_pdb_line(line):
    if line.startswith("TOM"):
        line = "A" + line
    return {
        'record':   line[:6].strip(),
        'serial':   int(line[6:11]),
        'atom_name': line[12:16],
        'alt_loc':  line[16],
        'res_name': line[17:20].strip(),
        'chain':    line[21],
        'res_seq':  int(line[22:26]),
        'x': float(line[30:38]),
        'y': float(line[38:46]),
        'z': float(line[46:54]),
        'line_rest': line[54:].rstrip(),
    }


def read_pdb_atoms(path):
    atoms, headers = [], []
    with open(path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM", "TOM")):
                try:
                    atoms.append(parse_pdb_line(line))
                except Exception:
                    continue
            elif line.startswith(("TER", "END", "CONECT", "MASTER")):
                continue
            else:
                headers.append(line)
    return atoms, headers


def write_atom_line(atom, serial, res_seq, res_name, coords):
    x, y, z = coords
    raw_name = atom['atom_name'].strip()
    if len(raw_name) < 4:
        atom_name_fmt = f" {raw_name:<3}"
    else:
        atom_name_fmt = f"{raw_name:<4}"
    chain_id = atom['chain'].strip() or " "
    chain_id = chain_id[0]
    return (
        f"{atom['record']:<6}"
        f"{serial:>5d}"
        " "
        f"{atom_name_fmt[:4]}"
        f"{atom['alt_loc']:1}"
        f"{res_name:>3}"
        " "
        f"{chain_id}"
        f"{res_seq:>4d}"
        "    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}"
        f"{atom['line_rest']}\n"
    )


def coords_array(atoms):
    return np.array([[a['x'], a['y'], a['z']] for a in atoms], dtype=float)


def principal_axis(xyz):
    """Return unit vector of the largest-variance axis (PCA)."""
    c = xyz.mean(axis=0)
    M = xyz - c
    _, _, vt = np.linalg.svd(M, full_matrices=False)
    return c, vt[0] / np.linalg.norm(vt[0])


def rotation_align(v_from, v_to):
    """3x3 rotation matrix aligning unit vector v_from onto v_to."""
    a = v_from / np.linalg.norm(v_from)
    b = v_to   / np.linalg.norm(v_to)
    v = np.cross(a, b)
    s = np.linalg.norm(v)
    c = float(np.dot(a, b))
    if s < 1e-8:
        if c > 0:
            return np.eye(3)
        # 180-degree flip; pick any perpendicular axis
        axis = np.array([1.0, 0.0, 0.0])
        if abs(a[0]) > 0.9:
            axis = np.array([0.0, 1.0, 0.0])
        axis -= a * float(np.dot(axis, a))
        axis /= np.linalg.norm(axis)
        return rotation_axis_angle(axis, math.pi)
    K = np.array([[0, -v[2], v[1]],
                  [v[2], 0, -v[0]],
                  [-v[1], v[0], 0]])
    return np.eye(3) + K + K @ K * ((1 - c) / (s * s))


def rotation_axis_angle(axis, angle):
    axis = axis / np.linalg.norm(axis)
    c, s = math.cos(angle), math.sin(angle)
    x, y, z = axis
    C = 1 - c
    return np.array([
        [c + x*x*C,    x*y*C - z*s, x*z*C + y*s],
        [y*x*C + z*s,  c + y*y*C,   y*z*C - x*s],
        [z*x*C - y*s,  z*y*C + x*s, c + z*z*C  ],
    ])


def random_rotation(rng):
    """Uniformly distributed rotation matrix (Shoemake / quaternion)."""
    u1, u2, u3 = rng.random(), rng.random(), rng.random()
    q = np.array([
        math.sqrt(1 - u1) * math.sin(2 * math.pi * u2),
        math.sqrt(1 - u1) * math.cos(2 * math.pi * u2),
        math.sqrt(u1)     * math.sin(2 * math.pi * u3),
        math.sqrt(u1)     * math.cos(2 * math.pi * u3),
    ])
    x, y, z, w = q
    return np.array([
        [1 - 2*(y*y + z*z), 2*(x*y - z*w),     2*(x*z + y*w)],
        [2*(x*y + z*w),     1 - 2*(x*x + z*z), 2*(y*z - x*w)],
        [2*(x*z - y*w),     2*(y*z + x*w),     1 - 2*(x*x + y*y)],
    ])


def min_pair_dist(cand, cage_xyz):
    """Smallest distance from any candidate atom to any cage atom."""
    diff = cand[:, None, :] - cage_xyz[None, :, :]
    d2 = (diff * diff).sum(axis=-1)
    return math.sqrt(float(d2.min()))


def cage_center_and_axis(cage_atoms):
    """Use Pt-Pt midpoint and vector if >=2 Pt atoms; else PCA centroid+axis."""
    pt_xyz = np.array([[a['x'], a['y'], a['z']] for a in cage_atoms
                       if a['atom_name'].strip().upper().startswith('PT')], dtype=float)
    xyz = coords_array(cage_atoms)
    if len(pt_xyz) >= 2:
        center = pt_xyz.mean(axis=0)
        axis = pt_xyz[-1] - pt_xyz[0]
        axis = axis / np.linalg.norm(axis)
        return center, axis, ('Pt-Pt midpoint', float(np.linalg.norm(pt_xyz[-1] - pt_xyz[0])))
    center, axis = principal_axis(xyz)
    return center, axis, ('PCA centroid', None)


def place_inside(cage_atoms, templ_atoms, clash=2.0, jitter=1.0,
                 axial_slide=2.0, attempts=20000, seed=None, verbose=True):
    rng = random.Random(seed)
    cage_xyz = coords_array(cage_atoms)
    tmpl_xyz = coords_array(templ_atoms)

    cage_c, cage_axis, (center_kind, pt_pt) = cage_center_and_axis(cage_atoms)
    tmpl_c, tmpl_axis  = principal_axis(tmpl_xyz)

    if verbose:
        d_cage = np.linalg.norm(cage_xyz - cage_c, axis=1)
        msg = f"  cavity center: {center_kind}; nearest cage atom: {d_cage.min():.2f} A"
        if pt_pt is not None:
            msg += f"; Pt-Pt = {pt_pt:.2f} A"
        print(msg)

    # base orientation: align PFOA long axis with cage long axis
    R_align = rotation_align(tmpl_axis, cage_axis)
    base = (tmpl_xyz - tmpl_c) @ R_align.T

    cutoffs = [clash]
    while cutoffs[-1] > 1.5:
        cutoffs.append(round(cutoffs[-1] - 0.1, 2))

    best_overall = None
    for cut in cutoffs:
        best = None
        for k in range(attempts):
            # rotations: 50% pure around-axis, 25% small tilt, 25% full random
            r = rng.random()
            if r < 0.50:
                R_extra = rotation_axis_angle(cage_axis, rng.uniform(0, 2 * math.pi))
            elif r < 0.75:
                # spin around axis + small tilt off axis (<=20 deg)
                R_spin = rotation_axis_angle(cage_axis, rng.uniform(0, 2 * math.pi))
                # random perpendicular tilt
                perp = np.cross(cage_axis, np.array([1.0, 0.0, 0.0]))
                if np.linalg.norm(perp) < 1e-6:
                    perp = np.cross(cage_axis, np.array([0.0, 1.0, 0.0]))
                perp /= np.linalg.norm(perp)
                R_tilt = rotation_axis_angle(perp, rng.uniform(-math.radians(20), math.radians(20)))
                R_extra = R_tilt @ R_spin
            else:
                R_extra = random_rotation(rng)
            cand = base @ R_extra.T
            # jitter: small radial offset + axial slide
            jx = rng.uniform(-jitter, jitter)
            jy = rng.uniform(-jitter, jitter)
            jz = rng.uniform(-jitter, jitter)
            slide = rng.uniform(-axial_slide, axial_slide) * cage_axis
            cand = cand + (cage_c + np.array([jx, jy, jz]) + slide)
            d = min_pair_dist(cand, cage_xyz)
            if d >= cut:
                if verbose:
                    print(f"  placed: cutoff={cut:.2f} A, min cage-guest dist={d:.2f} A, attempt {k+1}")
                return cand.tolist(), d, cut
            if best is None or d > best[1]:
                best = (cand.tolist(), d)
        if verbose:
            print(f"  no fit at cutoff {cut:.2f} A (best min-dist {best[1]:.2f} A); relaxing.")
        if best_overall is None or best[1] > best_overall[1]:
            best_overall = best
    return best_overall[0], best_overall[1], cutoffs[-1]


# ---------- Counterion placement (merged from seasoning.py) ----------
COUNTERION_CLASH    = 4.0
COUNTERION_BUFFER   = 6.0
COUNTERION_SHELL_W  = 15.0
COUNTERION_ATTEMPTS = 10000


def _geometry_stats(atoms):
    cx = sum(a['x'] for a in atoms) / len(atoms)
    cy = sum(a['y'] for a in atoms) / len(atoms)
    cz = sum(a['z'] for a in atoms) / len(atoms)
    r = max(math.sqrt((a['x']-cx)**2 + (a['y']-cy)**2 + (a['z']-cz)**2) for a in atoms)
    return (cx, cy, cz), r


def _mol2_resname(mol2_file):
    with open(mol2_file) as f:
        for line in f:
            if line.startswith("@<TRIPOS>MOLECULE"):
                return f.readline().strip()
    return None


def _clash_with(cand_coords, fixed_atoms, cutoff):
    c2 = cutoff * cutoff
    for nx, ny, nz in cand_coords:
        for a in fixed_atoms:
            dx, dy, dz = nx - a['x'], ny - a['y'], nz - a['z']
            if dx*dx + dy*dy + dz*dz < c2:
                return True
    return False


def _counterion_name_gen(used):
    """BFA, BFB, ..., BFZ, BGA, ..., BZZ (matches seasoning.name_generator)."""
    alpha = string.ascii_uppercase
    for c in alpha:
        name = f"BF{c}"
        if name not in used:
            yield name
    for c2 in alpha:
        if c2 <= 'F':
            continue
        for c3 in alpha:
            name = f"B{c2}{c3}"
            if name not in used:
                yield name


def place_counterions_outside(cage_atoms, current_atoms,
                               templ_pdb, templ_mol2, templ_frcmod,
                               n_anions, start_serial, start_res_seq,
                               seed=None, verbose=True):
    """Place N counterion copies in a shell around the cage.

    `current_atoms` is the running list of fixed atoms (cage + guest +
    previously-placed counterions); new placements are appended in-place.
    Returns (out_pdb_lines, last_serial, last_res_seq, placed_count).
    """
    rng = random.Random(seed if seed is not None else 0)
    templ_atoms, _ = read_pdb_atoms(templ_pdb)
    if not templ_atoms:
        raise RuntimeError(f"Counterion template empty: {templ_pdb}")

    tmpl_resname = _mol2_resname(templ_mol2)
    if not tmpl_resname:
        raise RuntimeError(f"Cannot read residue name from {templ_mol2}")

    cage_center, cage_radius = _geometry_stats(cage_atoms)
    tmpl_center, _ = _geometry_stats(templ_atoms)
    used = set(a['res_name'] for a in current_atoms)

    min_d = cage_radius + COUNTERION_BUFFER
    max_d = min_d + COUNTERION_SHELL_W

    gen = _counterion_name_gen(used)
    last_serial, last_res = start_serial, start_res_seq
    out_lines = []
    placed_count = 0

    for _ in range(n_anions):
        try:
            new_resname = next(gen)
        except StopIteration:
            break
        last_res += 1
        placed = False
        for _attempt in range(COUNTERION_ATTEMPTS):
            theta = rng.uniform(0, 2 * math.pi)
            phi = math.acos(rng.uniform(-1, 1))
            r = rng.uniform(min_d, max_d)
            rx = cage_center[0] + r * math.sin(phi) * math.cos(theta)
            ry = cage_center[1] + r * math.sin(phi) * math.sin(theta)
            rz = cage_center[2] + r * math.cos(phi)
            dx, dy, dz = rx - tmpl_center[0], ry - tmpl_center[1], rz - tmpl_center[2]
            cand = [(ta['x'] + dx, ta['y'] + dy, ta['z'] + dz) for ta in templ_atoms]
            if _clash_with(cand, current_atoms, COUNTERION_CLASH):
                continue
            ion_lines = []
            for i, ta in enumerate(templ_atoms):
                last_serial += 1
                cx, cy, cz = cand[i]
                a = ta.copy()
                a['alt_loc'] = ' '
                a['x'], a['y'], a['z'] = cx, cy, cz
                current_atoms.append(a)
                ion_lines.append(write_atom_line(a, last_serial, last_res,
                                                 new_resname, (cx, cy, cz)))
            out_lines.extend(ion_lines)
            out_lines.append("TER\n")

            # Per-anion {pdb, mol2, frcmod}. Skip writes that would overwrite
            # the user-supplied templates with placement-specific data.
            target_pdb = f"{new_resname}.pdb"
            if os.path.abspath(target_pdb) != os.path.abspath(templ_pdb):
                with open(target_pdb, "w") as f:
                    f.writelines(ion_lines)
                    f.write("TER\nEND\n")
            target_mol2 = f"{new_resname}.mol2"
            if os.path.abspath(target_mol2) != os.path.abspath(templ_mol2):
                with open(templ_mol2) as f:
                    content = f.read()
                content = re.sub(r'\b' + re.escape(tmpl_resname) + r'\b',
                                 new_resname, content)
                with open(target_mol2, "w") as f:
                    f.write(content)
            target_frcmod = f"{new_resname}.frcmod"
            if os.path.abspath(target_frcmod) != os.path.abspath(templ_frcmod):
                shutil.copy(templ_frcmod, target_frcmod)

            placed_count += 1
            placed = True
            if verbose and (placed_count % 5 == 0 or placed_count == n_anions):
                print(f"   ... counterions placed {placed_count}/{n_anions} ({new_resname})")
            break
        if not placed:
            if verbose:
                print(f"   counterion {new_resname}: no placement after "
                      f"{COUNTERION_ATTEMPTS} attempts; stopping.")
            break

    return out_lines, last_serial, last_res, placed_count


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--cage',    default='bone.pdb',     help='cage PDB (no guest)')
    ap.add_argument('--out',     default='tastybone.pdb',help='output PDB')
    ap.add_argument('--templ',   default='BFA.pdb',
                    help='guest template; .pdb or .xyz (.xyz auto-converted via OpenBabel)')
    ap.add_argument('--resname', default='GS1',          help='residue name for placed guest (default GS1)')
    ap.add_argument('--autoparam', action='store_true',
                    help='also run antechamber + parmchk2 to generate <stem>.mol2 and <stem>.frcmod')
    ap.add_argument('--charge',  type=int,   default=0,
                    help='net integer charge for antechamber when --autoparam (default 0)')
    ap.add_argument('--clash',   type=float, default=2.0, help='initial clash cutoff for guest placement (A)')
    ap.add_argument('--jitter',  type=float, default=1.0, help='max +/- radial offset from centroid (A)')
    ap.add_argument('--axial',   type=float, default=2.0, help='max +/- slide along Pt-Pt axis (A)')
    ap.add_argument('--attempts',type=int,   default=20000)
    ap.add_argument('--seed',    type=int,   default=0)
    # Counterion options (merged from seasoning).
    ap.add_argument('--counterions',        type=int, default=0,
                    help='also place N counterion copies around the cage (default 0; off)')
    ap.add_argument('--counterion-templ',   default='BFA.pdb',
                    help='counterion PDB template (default BFA.pdb); resname cycles BFA, BFB, ...')
    ap.add_argument('--counterion-mol2',    default='BFA.mol2',
                    help='counterion MOL2 template (default BFA.mol2)')
    ap.add_argument('--counterion-frcmod',  default='BFA.frcmod',
                    help='counterion FRCMOD template (default BFA.frcmod)')
    args = ap.parse_args()

    # 1. Resolve the guest template: xyz -> pdb via OpenBabel, pdb passes through.
    templ_path = args.templ
    if not os.path.isfile(templ_path):
        sys.stderr.write(f"Guest template not found: {templ_path}\n")
        sys.exit(1)
    stem, ext = os.path.splitext(os.path.basename(templ_path))
    ext_lower = ext.lower()
    if ext_lower == '.xyz':
        pdb_from_xyz = os.path.join(os.path.dirname(templ_path) or '.', f"{stem}.pdb")
        xyz_to_pdb(templ_path, pdb_from_xyz)
        templ_path = pdb_from_xyz
    elif ext_lower != '.pdb':
        sys.stderr.write(f"Unsupported guest format '{ext}'; use .pdb or .xyz\n")
        sys.exit(1)

    # 2. Optionally parametrize the guest (antechamber + parmchk2).
    if args.autoparam:
        print("Generating GAFF2 parameters for guest:")
        mol2, frcmod = run_autoparam(templ_path, stem, args.charge, args.resname)
        print(f"  -> {mol2}\n  -> {frcmod}")

    cage_atoms, cage_headers = read_pdb_atoms(args.cage)
    tmpl_atoms, _ = read_pdb_atoms(templ_path)
    if not cage_atoms or not tmpl_atoms:
        sys.stderr.write("Missing atoms in cage or template.\n")
        sys.exit(1)

    print(f"Cage atoms: {len(cage_atoms)}   guest template atoms: {len(tmpl_atoms)}")
    new_coords, min_d, cut_used = place_inside(
        cage_atoms, tmpl_atoms,
        clash=args.clash, jitter=args.jitter, axial_slide=args.axial,
        attempts=args.attempts, seed=args.seed,
    )
    if min_d < cut_used:
        print(f"WARNING: best-effort placement; min cage-guest distance {min_d:.2f} A < target {cut_used:.2f} A")

    last_serial = max((a['serial'] for a in cage_atoms), default=0)
    last_res    = max((a['res_seq'] for a in cage_atoms), default=0)

    # Materialise the guest's placed atoms (so counterion clash check sees them).
    guest_lines = []
    guest_atoms = []
    last_res += 1
    guest_res_seq = last_res
    for i, ta in enumerate(tmpl_atoms):
        last_serial += 1
        a = ta.copy()
        a['alt_loc'] = ' '
        cx, cy, cz = new_coords[i]
        a['x'], a['y'], a['z'] = cx, cy, cz
        guest_atoms.append(a)
        guest_lines.append(write_atom_line(a, last_serial, guest_res_seq,
                                           args.resname, (cx, cy, cz)))

    # Counterions (optional): place N copies in shells around the cage.
    counterion_lines = []
    placed_counterions = 0
    if args.counterions > 0:
        for need in (args.counterion_templ, args.counterion_mol2, args.counterion_frcmod):
            if not os.path.isfile(need):
                sys.stderr.write(f"Counterion template missing: {need}\n")
                sys.exit(1)
        print(f"Placing {args.counterions} counterions around the cage "
              f"(templ={args.counterion_templ}):")
        current = list(cage_atoms) + list(guest_atoms)
        counterion_lines, last_serial, last_res, placed_counterions = \
            place_counterions_outside(
                cage_atoms, current,
                args.counterion_templ, args.counterion_mol2, args.counterion_frcmod,
                args.counterions, last_serial, last_res, seed=args.seed,
            )

    with open(args.out, 'w') as f:
        for h in cage_headers: f.write(h)
        for a in cage_atoms:
            f.write(write_atom_line(a, a['serial'], a['res_seq'],
                                    a['res_name'], (a['x'], a['y'], a['z'])))
        f.write("TER\n")
        f.writelines(guest_lines)
        f.write("TER\n")
        if counterion_lines:
            f.writelines(counterion_lines)
        f.write("END\n")
    msg = (f"Wrote {args.out} with 1 {args.resname} inside "
           f"(min cage-guest dist={min_d:.2f} A)")
    if args.counterions > 0:
        msg += f" + {placed_counterions}/{args.counterions} counterions"
    print(msg + ".")


if __name__ == "__main__":
    main()
