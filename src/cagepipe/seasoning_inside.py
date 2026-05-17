#!/usr/bin/env python3
"""seasoning_inside.py - place ONE PFOA (BFA template) INSIDE the cage cavity.

Strategy:
  1. Align PFOA's longest principal axis with the cage's longest principal axis
     (PCA on heavy atoms) so the elongated chain lies along the natural channel.
  2. Try random rotations + small jitter at the cage centroid.
  3. Use a tight clash cutoff (default 2.5 A); relax progressively if no fit.

Inputs (same contract as seasoning.py):
  bone.pdb  + BFA.pdb / BFA.mol2 / BFA.frcmod
Output:
  tastybone.pdb   (cage + 1 BFA residue at the cavity)
"""
import argparse, math, os, random, sys
from itertools import combinations

try:
    import numpy as np
except ImportError:
    sys.stderr.write("numpy is required for seasoning_inside.py\n")
    sys.exit(2)


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
                    print(f"  placed: cutoff={cut:.2f} A, min cage-PFOA dist={d:.2f} A, attempt {k+1}")
                return cand.tolist(), d, cut
            if best is None or d > best[1]:
                best = (cand.tolist(), d)
        if verbose:
            print(f"  no fit at cutoff {cut:.2f} A (best min-dist {best[1]:.2f} A); relaxing.")
        if best_overall is None or best[1] > best_overall[1]:
            best_overall = best
    return best_overall[0], best_overall[1], cutoffs[-1]


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--cage',    default='bone.pdb',     help='cage PDB (no guest)')
    ap.add_argument('--out',     default='tastybone.pdb',help='output PDB')
    ap.add_argument('--templ',   default='BFA.pdb',      help='PFOA template PDB')
    ap.add_argument('--resname', default='BFA',          help='residue name for placed PFOA')
    ap.add_argument('--clash',   type=float, default=2.0, help='initial clash cutoff (A)')
    ap.add_argument('--jitter',  type=float, default=1.0, help='max +/- radial offset from centroid (A)')
    ap.add_argument('--axial',   type=float, default=2.0, help='max +/- slide along Pt-Pt axis (A)')
    ap.add_argument('--attempts',type=int,   default=20000)
    ap.add_argument('--seed',    type=int,   default=0)
    args = ap.parse_args()

    cage_atoms, cage_headers = read_pdb_atoms(args.cage)
    tmpl_atoms, _ = read_pdb_atoms(args.templ)
    if not cage_atoms or not tmpl_atoms:
        sys.stderr.write("Missing atoms in cage or template.\n")
        sys.exit(1)

    print(f"Cage atoms: {len(cage_atoms)}   PFOA template atoms: {len(tmpl_atoms)}")
    new_coords, min_d, cut_used = place_inside(
        cage_atoms, tmpl_atoms,
        clash=args.clash, jitter=args.jitter, axial_slide=args.axial,
        attempts=args.attempts, seed=args.seed,
    )
    if min_d < cut_used:
        print(f"WARNING: best-effort placement; min cage-PFOA distance {min_d:.2f} A < target {cut_used:.2f} A")

    last_serial = max((a['serial'] for a in cage_atoms), default=0)
    last_res    = max((a['res_seq'] for a in cage_atoms), default=0)

    new_lines = []
    for i, ta in enumerate(tmpl_atoms):
        last_serial += 1
        a = ta.copy()
        a['alt_loc'] = ' '
        cx, cy, cz = new_coords[i]
        new_lines.append(write_atom_line(a, last_serial, last_res + 1,
                                         args.resname, (cx, cy, cz)))

    with open(args.out, 'w') as f:
        for h in cage_headers: f.write(h)
        for a in cage_atoms:
            f.write(write_atom_line(a, a['serial'], a['res_seq'],
                                    a['res_name'], (a['x'], a['y'], a['z'])))
        f.write("TER\n")
        f.writelines(new_lines)
        f.write("TER\nEND\n")
    print(f"Wrote {args.out} with 1 {args.resname} inside (min cage-PFOA dist={min_d:.2f} A).")


if __name__ == "__main__":
    main()
