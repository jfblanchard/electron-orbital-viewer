/**
 * molecules-lib.js — put hybrid orbitals onto real molecular geometry.
 *
 * Pure JavaScript, no Three.js, so it runs in the browser and under Node.
 *
 * An orbital on a center is a 4-vector of coefficients over (2s, 2px, 2py, 2pz).
 * A hybrid pointing along unit vector d̂ is
 *
 *   h = a·s + b·(d̂x·px + d̂y·py + d̂z·pz),   a² + b² = 1
 *
 * so any direction is available, not just the cartesian axes. The s-fraction
 * a² follows from the geometry (Coulson's theorem): two hybrids at angle θ are
 * orthogonal when a₁a₂ + b₁b₂·cos θ = 0. At 109.47° that gives the textbook
 * a = ½ (sp³); at water's 104.5° the O–H hybrids come out 20% s and the lone
 * pairs pick up the rest, 30% s each.
 *
 * Sections
 *   1. Vectors       tiny 3-vector helpers
 *   2. Elements      CPK colors, masses, effective nuclear charge
 *   3. Hybrids       hybridAlong, Coulson s-fractions, lone-pair completion
 *   4. Molecule      atoms + bonds + centers → orbital list
 *   5. Canonical     every lobe as one cached mesh along +z, plus a rotation
 */

// ============================================================================
// 1. VECTORS
// ============================================================================

export const add = (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
export const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
export const scale = (a, s) => [a[0] * s, a[1] * s, a[2] * s];
export const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
export const cross = (a, b) => [
  a[1] * b[2] - a[2] * b[1],
  a[2] * b[0] - a[0] * b[2],
  a[0] * b[1] - a[1] * b[0],
];
export const norm = (a) => Math.hypot(a[0], a[1], a[2]);
export const unit = (a) => scale(a, 1 / norm(a));
export const angleDeg = (a, b) => (Math.acos(Math.max(-1, Math.min(1, dot(unit(a), unit(b))))) * 180) / Math.PI;

const dot4 = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + a[3] * b[3];
const scale4 = (a, s) => a.map((x) => x * s);
const sub4 = (a, b) => a.map((x, i) => x - b[i]);

// ============================================================================
// 2. ELEMENTS
// ============================================================================
//
// Orbitals are hydrogen-like with an effective nuclear charge Z. For the
// valence 2s/2p of C, N, O that is Slater's Z_eff; for hydrogen in a molecule
// it is 1.24, the standard contracted 1s exponent (the free atom's Z = 1 gives
// a 1s bigger than the C–H bond itself).

export const A0 = 0.529177; // Bohr radius in Å; orbital meshes come out in a0

export const ELEMENTS = {
  H: { n: 1, Z: 1.24, mass: 1.008, color: '#ffffff', radius: 0.2 },
  C: { n: 2, Z: 3.25, mass: 12.011, color: '#808080', radius: 0.3 },
  N: { n: 2, Z: 3.9, mass: 14.007, color: '#3050f8', radius: 0.29 },
  O: { n: 2, Z: 4.55, mass: 15.999, color: '#ff0000', radius: 0.28 },
};

// ============================================================================
// 3. HYBRIDS
// ============================================================================

/** s coefficient of N equivalent hybrids: sp 1/√2, sp² 1/√3, sp³ 1/2. */
export const EQUIVALENT_A = { sp: Math.SQRT1_2, sp2: 1 / Math.sqrt(3), sp3: 0.5 };

/** Coefficient 4-vector (s, px, py, pz) of a hybrid along unit vector d. */
export function hybridVector(dir, a) {
  const b = Math.sqrt(1 - a * a);
  return [a, b * dir[0], b * dir[1], b * dir[2]];
}

/** 4-vector → orbitals.js terms [[coef, n, ang]], dropping zero terms. */
export function termsOf(vec, n = 2) {
  const keys = ['s', 'px', 'py', 'pz'];
  return vec.flatMap((c, i) => (Math.abs(c) > 1e-12 ? [[c, n, keys[i]]] : []));
}

/** orbitals.js terms for a hybrid along unit vector d with s coefficient a. */
export function hybridAlong(dir, a = EQUIVALENT_A.sp3, n = 2) {
  return termsOf(hybridVector(dir, a), n);
}

/**
 * s coefficients for hybrids aimed along the given directions, from
 * Coulson's orthogonality condition. With x = a/b, hybrids i and j are
 * orthogonal when x_i·x_j = −cos θ_ij. Two directions are taken as
 * equivalent (x² = −cos θ); three or more are solved exactly from any
 * triple, x_i² = c_ij·c_ik / c_jk.
 */
export function coulsonA(dirs) {
  const c = (i, j) => -dot(dirs[i], dirs[j]);
  const x = dirs.map((_, i) => {
    if (dirs.length === 2) return Math.sqrt(c(0, 1));
    const [j, k] = dirs.map((__, q) => q).filter((q) => q !== i);
    return Math.sqrt((c(i, j) * c(i, k)) / c(j, k));
  });
  return x.map((xi) => xi / Math.sqrt(1 + xi * xi));
}

/**
 * Lone pairs that complete an sp³ set: the orbitals orthogonal to the given
 * bonding hybrids. One lone pair (NH₃) is the unique leftover direction. Two
 * (H₂O) are split symmetrically above and below the bond plane.
 */
export function completeLonePairs(hybrids, count, bondDirs) {
  const basis = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]];
  const orthogonalize = (v, against) => {
    let r = v;
    for (const u of against) r = sub4(r, scale4(u, dot4(r, u)));
    return r;
  };
  // the leftover vector orthogonal to `against` with the largest remainder
  const leftover = (against) => {
    let best = null;
    for (const e of basis) {
      const r = orthogonalize(e, against);
      if (!best || dot4(r, r) > dot4(best, best)) best = r;
    }
    const v = scale4(best, 1 / Math.sqrt(dot4(best, best)));
    return v[0] < 0 ? scale4(v, -1) : v; // s > 0: big lobe along the p part
  };
  if (count === 1) return [leftover(hybrids)];
  if (count === 2) {
    const nrm = unit(cross(bondDirs[0], bondDirs[1]));
    const perp = [0, ...nrm];
    const inPlane = leftover([...hybrids, perp]);
    const k = Math.SQRT1_2;
    return [
      inPlane.map((x, i) => k * (x + perp[i])),
      inPlane.map((x, i) => k * (x - perp[i])),
    ];
  }
  throw new Error(`unsupported lone pair count ${count}`);
}

/** Split a 4-vector into s coefficient a and unit p direction. */
export function decompose(vec) {
  const p = [vec[1], vec[2], vec[3]];
  const len = norm(p);
  return { a: vec[0], dir: len > 1e-12 ? scale(p, 1 / len) : null };
}

// ============================================================================
// 4. MOLECULE
// ============================================================================
//
// Input: atoms [element, [x, y, z] Å], bonds [i, j, order, perp?], and
// hybridized centers:
//   { atom, sigma: [neighbor atom indices],
//     lonePairs: count (derived by orthogonality) or [directions],
//     pi: [{ dir, partner }] }  leftover pure p orbitals
// Hydrogens get a 1s automatically. Output is centered on the center of mass
// and carries a flat `orbitals` list the viewer and tests both consume.

export function buildMolecule(meta, atomList, bondList, centers) {
  let atoms = atomList.map(([el, pos]) => ({ el, pos }));
  let M = 0;
  let com = [0, 0, 0];
  for (const a of atoms) {
    const m = ELEMENTS[a.el].mass;
    M += m;
    com = add(com, scale(a.pos, m));
  }
  com = scale(com, 1 / M);
  atoms = atoms.map((a) => ({ ...a, pos: sub(a.pos, com) }));

  const bonds = bondList.map(([a, b, order = 1, perp = null]) => ({ a, b, order, perp }));
  const orbitals = [];

  for (const c of centers) {
    const at = atoms[c.atom];
    const { n, Z } = ELEMENTS[at.el];
    const bondDirs = c.sigma.map((j) => unit(sub(atoms[j].pos, at.pos)));
    const lpDirs = Array.isArray(c.lonePairs) ? c.lonePairs.map(unit) : [];
    const domains = [...bondDirs, ...lpDirs];
    const as = coulsonA(domains);
    const vecs = domains.map((d, i) => hybridVector(d, as[i]));
    const base = { atom: c.atom, center: at.pos, n, Z };
    vecs.forEach((vec, i) => {
      const isBond = i < bondDirs.length;
      orbitals.push({
        ...base, vec, ...decompose(vec),
        role: isBond ? 'sigma' : 'lp',
        partner: isBond ? c.sigma[i] : null,
      });
    });
    if (typeof c.lonePairs === 'number' && c.lonePairs > 0) {
      for (const vec of completeLonePairs(vecs, c.lonePairs, bondDirs)) {
        orbitals.push({ ...base, vec, ...decompose(vec), role: 'lp', partner: null });
      }
    }
    for (const p of c.pi || []) {
      const vec = [0, ...unit(p.dir)];
      orbitals.push({ ...base, vec, ...decompose(vec), role: 'pi', partner: p.partner });
    }
  }

  atoms.forEach((at, i) => {
    if (at.el !== 'H') return;
    const { n, Z } = ELEMENTS.H;
    const bond = bonds.find((b) => b.a === i || b.b === i);
    orbitals.push({
      atom: i, center: at.pos, n, Z, vec: [1, 0, 0, 0], a: 1, dir: null,
      role: 'h1s', partner: bond ? (bond.a === i ? bond.b : bond.a) : null,
    });
  });

  for (const o of orbitals) o.terms = termsOf(o.vec, o.n);
  return { ...meta, atoms, bonds, orbitals };
}

// ============================================================================
// 5. CANONICAL LOBES
// ============================================================================
//
// Every orbital here is symmetric about its own axis, so a hybrid along d̂ is
// the same shape as the one along +z, just turned. The viewer meshes each
// distinct shape once, keyed by (n, Z, a), and places copies with a rotation.
// That makes a whole molecule cost a handful of meshes instead of one per lobe.

export function canonical(o) {
  const n = o.n;
  if (!o.dir) return { key: `s|${n}|${o.Z}`, n, Z: o.Z, terms: [[1, n, 's']] };
  const a = Math.abs(o.a) < 1e-9 ? 0 : o.a;
  const b = Math.sqrt(1 - a * a);
  const terms = a === 0 ? [[1, n, 'pz']] : [[a, n, 's'], [b, n, 'pz']];
  return { key: `h|${n}|${o.Z}|${a.toFixed(5)}`, n, Z: o.Z, terms };
}

/** Rotation matrix (rows) taking +z to unit vector d (Rodrigues). */
export function rotationZTo(d) {
  const [x, y, z] = d;
  if (z < -1 + 1e-12) return [[1, 0, 0], [0, -1, 0], [0, 0, -1]];
  const k = 1 / (1 + z);
  return [
    [1 - x * x * k, -x * y * k, x],
    [-x * y * k, 1 - y * y * k, y],
    [-x, -y, z],
  ];
}

export const apply3 = (R, v) => [dot(R[0], v), dot(R[1], v), dot(R[2], v)];
