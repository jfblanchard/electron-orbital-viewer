/**
 * orbitals.js — hydrogen-like orbitals, hybrids, and isosurface extraction.
 *
 * Pure JavaScript with no Three.js dependency, so the same module runs in the
 * browser, in a Web Worker (worker.js), and under Node for tests.
 *
 * Units: distances in Bohr radii (a0 = 0.529 Å). Every wavefunction is
 * normalized so that ∫|ψ|² dV = 1 over all space.
 *
 * Sections
 *   1. Radial part       R_nl(r)            (associated Laguerre polynomials)
 *   2. Angular part      real Y_lm(θ, φ)    (written in Cartesian form)
 *   3. Orbitals          ψ = R · Y, and linear combinations (hybrids)
 *   4. Catalog           every orbital the viewer offers, with descriptions
 *   5. Isosurface        sample ψ on a grid, find the 90% contour, mesh it
 */

// ============================================================================
// 1. RADIAL PART
// ============================================================================
//
// For a one-electron atom with nuclear charge Z, the radial function is
//
//   R_nl(r) = N · e^(−ρ/2) · ρ^l · L_(n−l−1)^(2l+1)(ρ),     ρ = 2Zr / n
//
//   N = sqrt( (2Z/n)³ · (n−l−1)! / (2n · (n+l)!) )
//
// where L_k^α is the generalized (associated) Laguerre polynomial. The
// polynomial has n−l−1 roots, which are the radial nodes: 2s has one
// spherical node, 3s has two, 3p has one, 3d and 4f have none.
//
// Phase convention: the overall sign of any orbital is arbitrary. We multiply
// by (−1)^(n−l−1) so every radial function is POSITIVE at large r (the
// outermost lobe). Two consequences: 3pz's big outer lobe is red on +z just
// like 2pz, and s + p hybrids put their big lobe along the p direction.
// With the raw Laguerre sign, 2s and 3p are negative far out and the
// hybrids' big lobes would point backwards.

function factorial(n) {
  let f = 1;
  for (let i = 2; i <= n; i++) f *= i;
  return f;
}

/** Generalized Laguerre polynomial L_k^alpha(x), by the standard recurrence. */
export function laguerre(k, alpha, x) {
  if (k === 0) return 1;
  let prev = 1;
  let cur = 1 + alpha - x;
  for (let j = 1; j < k; j++) {
    const next = ((2 * j + 1 + alpha - x) * cur - (j + alpha) * prev) / (j + 1);
    prev = cur;
    cur = next;
  }
  return cur;
}

/** Returns R_nl as a function of r (in a0). */
export function radialFn(n, l, Z = 1) {
  const k = n - l - 1;
  const alpha = 2 * l + 1;
  const phase = k % 2 === 0 ? 1 : -1; // positive at large r, see above
  const norm = phase * Math.sqrt((2 * Z / n) ** 3 * factorial(k) / (2 * n * factorial(n + l)));
  return (r) => {
    const rho = 2 * Z * r / n;
    return norm * Math.exp(-rho / 2) * rho ** l * laguerre(k, alpha, rho);
  };
}

/**
 * Radius (a0) that encloses a given fraction of the radial probability
 * ∫ r² R² dr. Used to size the sampling grid so it holds essentially all
 * of the electron.
 */
export function radialExtent(n, l, fraction = 0.998, Z = 1) {
  const R = radialFn(n, l, Z);
  const dr = 0.01 / Z;
  let cum = 0;
  let r = 0;
  while (cum < fraction && r < 400) {
    const rm = r + dr / 2; // midpoint rule
    const Rm = R(rm);
    cum += rm * rm * Rm * Rm * dr;
    r += dr;
  }
  return r;
}

// ============================================================================
// 2. ANGULAR PART — real spherical harmonics
// ============================================================================
//
// The complex harmonics Y_l^m carry a factor e^(imφ). Chemistry uses real
// combinations of +m and −m instead, which is what gives px, py, pz, the five
// d orbitals, and the seven f orbitals their familiar lobes along the axes.
//
// Written in terms of the unit vector (x, y, z) = r̂, each real harmonic is a
// polynomial times a normalization constant chosen so ∫ Y² dΩ = 1. The zeros
// of each polynomial are the angular nodes (planes or cones); there are
// always exactly l of them.

const PI = Math.PI;
const c0 = 0.5 * Math.sqrt(1 / PI);
const c1 = Math.sqrt(3 / (4 * PI));
const c2a = 0.5 * Math.sqrt(15 / PI);         // dxy, dyz, dxz
const c2b = 0.25 * Math.sqrt(5 / PI);         // dz²
const c2c = 0.25 * Math.sqrt(15 / PI);        // dx²−y²
const c3a = 0.25 * Math.sqrt(35 / (2 * PI));  // fx(x²−3y²), fy(3x²−y²)
const c3b = 0.5 * Math.sqrt(105 / PI);        // fxyz
const c3c = 0.25 * Math.sqrt(21 / (2 * PI));  // fxz², fyz²
const c3d = 0.25 * Math.sqrt(7 / PI);         // fz³
const c3e = 0.25 * Math.sqrt(105 / PI);       // fz(x²−y²)

/** Real spherical harmonics, keyed by orbital suffix. f(x, y, z) takes a unit vector. */
export const ANGULAR = {
  s:        { l: 0, f: () => c0 },

  px:       { l: 1, f: (x) => c1 * x },
  py:       { l: 1, f: (x, y) => c1 * y },
  pz:       { l: 1, f: (x, y, z) => c1 * z },

  dxy:      { l: 2, f: (x, y) => c2a * x * y },
  dyz:      { l: 2, f: (x, y, z) => c2a * y * z },
  dxz:      { l: 2, f: (x, y, z) => c2a * x * z },
  dz2:      { l: 2, f: (x, y, z) => c2b * (3 * z * z - 1) },
  dx2y2:    { l: 2, f: (x, y) => c2c * (x * x - y * y) },

  fz3:      { l: 3, f: (x, y, z) => c3d * z * (5 * z * z - 3) },
  fxz2:     { l: 3, f: (x, y, z) => c3c * x * (5 * z * z - 1) },
  fyz2:     { l: 3, f: (x, y, z) => c3c * y * (5 * z * z - 1) },
  fxyz:     { l: 3, f: (x, y, z) => c3b * x * y * z },
  fzx2y2:   { l: 3, f: (x, y, z) => c3e * z * (x * x - y * y) },
  fxx23y2:  { l: 3, f: (x, y) => c3a * x * (x * x - 3 * y * y) },
  fy3x2y2:  { l: 3, f: (x, y) => c3a * y * (3 * x * x - y * y) },
};

// ============================================================================
// 3. ORBITALS — ψ = R(r)·Y(r̂), and linear combinations
// ============================================================================

/** A single hydrogen-like orbital ψ_n,ang(x, y, z). */
export function atomicOrbital(n, ang, Z = 1) {
  const { l, f } = ANGULAR[ang];
  const R = radialFn(n, l, Z);
  const center = l === 0 ? R(0) * c0 : 0; // at r = 0 only s orbitals are nonzero
  return (x, y, z) => {
    const r = Math.sqrt(x * x + y * y + z * z);
    if (r < 1e-12) return center;
    return R(r) * f(x / r, y / r, z / r);
  };
}

/**
 * A normalized linear combination Σ cᵢ ψᵢ, with terms given as
 * [coefficient, n, angularKey]. Because the atomic orbitals are orthonormal,
 * the combination is normalized whenever Σ cᵢ² = 1.
 */
export function combination(terms, Z = 1) {
  const parts = terms.map(([c, n, ang]) => [c, atomicOrbital(n, ang, Z)]);
  if (parts.length === 1 && parts[0][0] === 1) return parts[0][1];
  return (x, y, z) => {
    let sum = 0;
    for (let i = 0; i < parts.length; i++) sum += parts[i][0] * parts[i][1](x, y, z);
    return sum;
  };
}

/**
 * Equivalent spᵏ hybrids pointing along unit vectors d̂. Each is
 *
 *   h = a·s + b·(d̂x·px + d̂y·py + d̂z·pz),   a = 1/√N,  b = √(1 − 1/N)
 *
 * for N equivalent hybrids. The p-part is just a p orbital pointing along d̂.
 * Two such hybrids overlap by a² + b²·cos θ, which is zero exactly at the
 * textbook angles: 180° (sp), 120° (sp²), 109.47° (sp³).
 */
function spHybrids(n, directions) {
  const N = directions.length;
  const a = 1 / Math.sqrt(N);
  const b = Math.sqrt(1 - 1 / N);
  return directions.map(([dx, dy, dz]) => {
    const terms = [[a, n, 's']];
    if (Math.abs(dx) > 1e-12) terms.push([b * dx, n, 'px']);
    if (Math.abs(dy) > 1e-12) terms.push([b * dy, n, 'py']);
    if (Math.abs(dz) > 1e-12) terms.push([b * dz, n, 'pz']);
    return terms;
  });
}

// ============================================================================
// 4. CATALOG
// ============================================================================
//
// Each entry has a list of functions ("set"). Atomic orbitals have one;
// hybrid entries list every hybrid in the set so they can be shown together.

const deg = PI / 180;
const tetraPolar = Math.acos(-1 / 3); // 109.47°, angle between sp³ hybrids

// Azimuths 90°, 210°, 330°: one lobe along +y, and none aimed at the default
// camera (which looks in from roughly +x), so no back lobe hides the center.
const AZIMUTHS = [90, 210, 330];
const TRIGONAL = AZIMUTHS.map((a) => [Math.cos(a * deg), Math.sin(a * deg), 0]);
const TETRAHEDRAL = [
  [0, 0, 1],
  ...AZIMUTHS.map((a) => [
    Math.sin(tetraPolar) * Math.cos(a * deg),
    Math.sin(tetraPolar) * Math.sin(a * deg),
    Math.cos(tetraPolar),
  ]),
];

const r2 = Math.SQRT2;
const r3 = Math.sqrt(3);
const r6 = Math.sqrt(6);

const atomic = (id, group, n, ang, label, desc, extra = {}) => ({
  id, group, label, desc,
  kind: 'atomic',
  n, l: ANGULAR[ang].l,
  set: [{ name: id, terms: [[1, n, ang]] }],
  ...extra,
});

const hybrid = (id, label, desc, termsList, extra) => ({
  id, group: 'hybrid', label, desc,
  kind: 'hybrid',
  set: termsList.map((terms, i) => ({ name: `${id} #${i + 1}`, terms })),
  ...extra,
});

export const ORBITALS = [
  // --- s ---------------------------------------------------------------
  atomic('1s', 's', 1, 's', '1s', 'Sphere centered on the nucleus'),
  atomic('2s', 's', 2, 's', '2s', 'Sphere inside a larger shell, opposite sign; one spherical node'),
  atomic('3s', 's', 3, 's', '3s', 'Nested spheres separated by two spherical nodes'),

  // --- p ---------------------------------------------------------------
  atomic('2px', 'p', 2, 'px', '2p<sub>x</sub>', 'Dumbbell along the x-axis; nodal plane yz'),
  atomic('2py', 'p', 2, 'py', '2p<sub>y</sub>', 'Dumbbell along the y-axis; nodal plane xz'),
  atomic('2pz', 'p', 2, 'pz', '2p<sub>z</sub>', 'Dumbbell along the z-axis; nodal plane xy'),
  atomic('3px', 'p', 3, 'px', '3p<sub>x</sub>', 'Dumbbell along x with an inner lobe pair; one spherical node'),
  atomic('3py', 'p', 3, 'py', '3p<sub>y</sub>', 'Dumbbell along y with an inner lobe pair; one spherical node'),
  atomic('3pz', 'p', 3, 'pz', '3p<sub>z</sub>', 'Dumbbell along z with an inner lobe pair; one spherical node'),

  // --- d ---------------------------------------------------------------
  atomic('3dz2', 'd', 3, 'dz2', '3d<sub>z²</sub>', 'Dumbbell along z with a ring (torus) around the middle; two nodal cones'),
  atomic('3dx2y2', 'd', 3, 'dx2y2', '3d<sub>x²−y²</sub>', 'Four-leaf clover with lobes on the x and y axes'),
  atomic('3dxy', 'd', 3, 'dxy', '3d<sub>xy</sub>', 'Four-leaf clover in the xy plane, lobes between the axes'),
  atomic('3dxz', 'd', 3, 'dxz', '3d<sub>xz</sub>', 'Four-leaf clover in the xz plane, lobes between the axes'),
  atomic('3dyz', 'd', 3, 'dyz', '3d<sub>yz</sub>', 'Four-leaf clover in the yz plane, lobes between the axes'),

  // --- f ---------------------------------------------------------------
  atomic('4fz3', 'f', 4, 'fz3', '4f<sub>z³</sub>', 'Lobes along z with two rings; two nodal cones and the xy plane'),
  atomic('4fxz2', 'f', 4, 'fxz2', '4f<sub>xz²</sub>', 'Six lobes in the xz plane; nodal plane yz and two nodal cones'),
  atomic('4fyz2', 'f', 4, 'fyz2', '4f<sub>yz²</sub>', 'Six lobes in the yz plane; nodal plane xz and two nodal cones'),
  atomic('4fxyz', 'f', 4, 'fxyz', '4f<sub>xyz</sub>', 'Eight lobes toward the corners of a cube; three nodal planes'),
  atomic('4fzx2y2', 'f', 4, 'fzx2y2', '4f<sub>z(x²−y²)</sub>', 'Eight lobes, like fxyz rotated 45° about z'),
  atomic('4fxx23y2', 'f', 4, 'fxx23y2', '4f<sub>x(x²−3y²)</sub>', 'Six lobes around the xy plane, one on +x'),
  atomic('4fy3x2y2', 'f', 4, 'fy3x2y2', '4f<sub>y(3x²−y²)</sub>', 'Six lobes around the xy plane, one on +y'),

  // --- hybrids ---------------------------------------------------------
  hybrid('sp', 'sp', 'Two hybrids pointing opposite ways, 180° apart (linear)',
    spHybrids(2, [[1, 0, 0], [-1, 0, 0]]),
    { formula: 'h = (2s ± 2p<sub>x</sub>) / √2', geometry: 'Linear · 180°', mix: '2s + 2p' }),

  hybrid('sp2', 'sp²', 'Three hybrids in a plane, 120° apart (trigonal planar)',
    spHybrids(2, TRIGONAL),
    { formula: 'h = (1/√3)·2s + √(2/3)·2p, with the 2p aimed along the lobe', geometry: 'Trigonal planar · 120°', mix: '2s + two 2p',
      view: [0.55, 0.25, 1] }), // look down on the plane

  hybrid('sp3', 'sp³', 'Four hybrids toward the corners of a tetrahedron, 109.5° apart',
    spHybrids(2, TETRAHEDRAL),
    { formula: 'h = ½·2s + (√3/2)·2p, with the 2p aimed along the lobe', geometry: 'Tetrahedral · 109.5°', mix: '2s + three 2p' }),

  hybrid('sp3d', 'sp³d', 'Three equatorial hybrids at 120° plus two axial ones (trigonal bipyramidal)',
    [
      ...spHybrids(3, TRIGONAL),
      [[1 / r2, 3, 'dz2'], [1 / r2, 3, 'pz']],
      [[1 / r2, 3, 'dz2'], [-1 / r2, 3, 'pz']],
    ],
    { formula: 'equatorial: (1/√3)·3s + √(2/3)·3p · axial: (3d<sub>z²</sub> ± 3p<sub>z</sub>) / √2',
      geometry: 'Trigonal bipyramidal · 90° / 120°', mix: '3s + three 3p + 3d' }),

  hybrid('sp3d2', 'sp³d²', 'Six hybrids along ±x, ±y, ±z (octahedral)',
    [
      [[1 / r6, 3, 's'], [ 1 / r2, 3, 'px'], [-1 / (2 * r3), 3, 'dz2'], [ 0.5, 3, 'dx2y2']],
      [[1 / r6, 3, 's'], [-1 / r2, 3, 'px'], [-1 / (2 * r3), 3, 'dz2'], [ 0.5, 3, 'dx2y2']],
      [[1 / r6, 3, 's'], [ 1 / r2, 3, 'py'], [-1 / (2 * r3), 3, 'dz2'], [-0.5, 3, 'dx2y2']],
      [[1 / r6, 3, 's'], [-1 / r2, 3, 'py'], [-1 / (2 * r3), 3, 'dz2'], [-0.5, 3, 'dx2y2']],
      [[1 / r6, 3, 's'], [ 1 / r2, 3, 'pz'], [ 1 / r3, 3, 'dz2']],
      [[1 / r6, 3, 's'], [-1 / r2, 3, 'pz'], [ 1 / r3, 3, 'dz2']],
    ],
    { formula: 'h = (1/√6)·3s ± (1/√2)·3p + d-part pointing the same way',
      geometry: 'Octahedral · 90°', mix: '3s + three 3p + two 3d' }),
];

export const ORBITALS_BY_ID = Object.fromEntries(ORBITALS.map((o) => [o.id, o]));

// ============================================================================
// 5. ISOSURFACE
// ============================================================================
//
// Textbook orbital pictures draw the surface that encloses ~90% of the
// probability of finding the electron. We find it numerically:
//
//   a. Sample ψ on an N×N×N grid spanning [−L, L]³, where L encloses
//      essentially all of the radial probability.
//   b. Sort the samples by |ψ|² and take the densest ones until they hold
//      90% of the total. The |ψ|² of the last one taken is the contour level.
//      The region |ψ|² > t is the smallest volume holding 90%.
//   c. |ψ|² = t means ψ = +√t or ψ = −√t. Mesh each sign separately, which
//      gives the red (+) and blue (−) lobes as independent closed surfaces.
//
// Meshing uses "surface nets": one vertex per grid cell the surface passes
// through, joined by quads across every grid edge the surface crosses. Each
// vertex is then pushed exactly onto the true surface with a few Newton steps
// on ψ itself, and its normal comes from ∇ψ. So the shape is exact to
// floating point, not a blocky voxel approximation; the grid only sets how
// many triangles are used.
//
// Known limit: where a lobe gets thinner than one grid cell (the rim of the
// cup-shaped outer lobe of 3p, some n = 3 hybrids) a single vertex can serve
// both walls, leaving a handful of non-manifold edges. They are invisible in
// the translucent render.

/** Sample fn on an N³ grid over [−L, L]³. Index = i + N·(j + N·k). */
export function sampleGrid(fn, N, L) {
  const h = 2 * L / (N - 1);
  const field = new Float32Array(N * N * N);
  let p = 0;
  for (let k = 0; k < N; k++) {
    const z = -L + k * h;
    for (let j = 0; j < N; j++) {
      const y = -L + j * h;
      for (let i = 0; i < N; i++) field[p++] = fn(-L + i * h, y, z);
    }
  }
  return field;
}

/**
 * Contour level c such that the region |ψ| > c holds `fraction` of the
 * probability captured by the grid. Returns { c, captured } where captured is
 * Σ|ψ|²·dV over the whole grid (≈ 1 when the grid is large enough).
 */
export function contourLevel(field, h, fraction = 0.9) {
  const dens = new Float32Array(field.length);
  let total = 0;
  for (let i = 0; i < field.length; i++) {
    const d = field[i] * field[i];
    dens[i] = d;
    total += d;
  }
  dens.sort(); // ascending
  const target = fraction * total;
  let cum = 0;
  let t = 0;
  for (let i = dens.length - 1; i >= 0; i--) {
    cum += dens[i];
    if (cum >= target) { t = dens[i]; break; }
  }
  return { c: Math.sqrt(t), captured: total * h * h * h };
}

// Cube corners (bit b of the index = offset along axis b) and the 12 edges.
const CORNERS = [
  [0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0],
  [0, 0, 1], [1, 0, 1], [0, 1, 1], [1, 1, 1],
];
const EDGES = [
  [0, 1], [2, 3], [4, 5], [6, 7],
  [0, 2], [1, 3], [4, 6], [5, 7],
  [0, 4], [1, 5], [2, 6], [3, 7],
];

/**
 * Mesh the surface sign·ψ = c, with sign = +1 or −1.
 * Inside is where sign·ψ > c; triangles wind counter-clockwise seen from
 * outside and normals point outward.
 * Returns { positions, normals, indices } as typed arrays.
 */
export function extractSurface(field, N, L, c, sign, fn) {
  const h = 2 * L / (N - 1);
  const M = N - 1; // cells per axis
  const g = (i, j, k) => sign * field[i + N * (j + N * k)] - c;

  // --- vertices: one per cell that the surface crosses ---
  const cellVert = new Int32Array(M * M * M).fill(-1);
  const pos = [];
  const v = new Float64Array(8);
  for (let k = 0; k < M; k++) {
    for (let j = 0; j < M; j++) {
      for (let i = 0; i < M; i++) {
        let mask = 0;
        for (let q = 0; q < 8; q++) {
          const [cx, cy, cz] = CORNERS[q];
          v[q] = g(i + cx, j + cy, k + cz);
          if (v[q] > 0) mask |= 1 << q;
        }
        if (mask === 0 || mask === 255) continue;
        // average of the points where the surface crosses the cell's edges
        let sx = 0, sy = 0, sz = 0, n = 0;
        for (const [a, b] of EDGES) {
          if ((v[a] > 0) === (v[b] > 0)) continue;
          const t = v[a] / (v[a] - v[b]);
          const A = CORNERS[a], B = CORNERS[b];
          sx += A[0] + t * (B[0] - A[0]);
          sy += A[1] + t * (B[1] - A[1]);
          sz += A[2] + t * (B[2] - A[2]);
          n++;
        }
        cellVert[i + M * (j + M * k)] = pos.length / 3;
        pos.push(-L + (i + sx / n) * h, -L + (j + sy / n) * h, -L + (k + sz / n) * h);
      }
    }
  }

  // --- faces: one quad per grid edge with a sign change ---
  // For an edge along axis a, the four cells sharing it differ in the other
  // two axes u, v. Taking (a, u, v) as a cyclic permutation makes the quad
  // (u−,v−) → (u+,v−) → (u+,v+) → (u−,v+) face +a, which is outward when
  // the edge starts inside. Otherwise the quad is reversed.
  const idx = [];
  const P = [0, 0, 0];
  const C = [0, 0, 0];
  const cell = (du, dv, a, u, vv) => {
    C[a] = P[a];
    C[u] = P[u] - 1 + du;
    C[vv] = P[vv] - 1 + dv;
    return cellVert[C[0] + M * (C[1] + M * C[2])];
  };
  const dist2 = (p, q) => {
    const dx = pos[3 * p] - pos[3 * q], dy = pos[3 * p + 1] - pos[3 * q + 1], dz = pos[3 * p + 2] - pos[3 * q + 2];
    return dx * dx + dy * dy + dz * dz;
  };
  for (let k = 0; k < N; k++) {
    for (let j = 0; j < N; j++) {
      for (let i = 0; i < N; i++) {
        P[0] = i; P[1] = j; P[2] = k;
        const in0 = g(i, j, k) > 0;
        for (let a = 0; a < 3; a++) {
          const u = (a + 1) % 3, vv = (a + 2) % 3;
          if (P[a] >= M || P[u] < 1 || P[u] > M - 1 || P[vv] < 1 || P[vv] > M - 1) continue;
          const in1 = a === 0 ? g(i + 1, j, k) > 0 : a === 1 ? g(i, j + 1, k) > 0 : g(i, j, k + 1) > 0;
          if (in0 === in1) continue;
          const q00 = cell(0, 0, a, u, vv), q10 = cell(1, 0, a, u, vv);
          const q11 = cell(1, 1, a, u, vv), q01 = cell(0, 1, a, u, vv);
          const quad = in0 ? [q00, q10, q11, q01] : [q00, q01, q11, q10];
          // split along the shorter diagonal for better-shaped triangles
          if (dist2(quad[0], quad[2]) <= dist2(quad[1], quad[3])) {
            idx.push(quad[0], quad[1], quad[2], quad[0], quad[2], quad[3]);
          } else {
            idx.push(quad[0], quad[1], quad[3], quad[1], quad[2], quad[3]);
          }
        }
      }
    }
  }

  // --- snap vertices onto the exact surface and compute normals from ∇ψ ---
  const positions = new Float32Array(pos);
  const normals = new Float32Array(pos.length);
  const eps = h * 1e-3;
  const f = (x, y, z) => sign * fn(x, y, z) - c;
  const grad = (x, y, z, out) => {
    out[0] = (f(x + eps, y, z) - f(x - eps, y, z)) / (2 * eps);
    out[1] = (f(x, y + eps, z) - f(x, y - eps, z)) / (2 * eps);
    out[2] = (f(x, y, z + eps) - f(x, y, z - eps)) / (2 * eps);
  };
  const gr = [0, 0, 0];
  for (let p = 0; p < positions.length; p += 3) {
    let x = positions[p], y = positions[p + 1], z = positions[p + 2];
    for (let it = 0; it < 4; it++) {
      const val = f(x, y, z);
      grad(x, y, z, gr);
      const g2 = gr[0] * gr[0] + gr[1] * gr[1] + gr[2] * gr[2];
      if (g2 < 1e-30) break;
      let s = val / g2;
      const step = Math.abs(s) * Math.sqrt(g2);
      if (step > 0.5 * h) s *= (0.5 * h) / step; // stay near the cell
      x -= s * gr[0]; y -= s * gr[1]; z -= s * gr[2];
      if (step < 1e-7 * h) break;
    }
    positions[p] = x; positions[p + 1] = y; positions[p + 2] = z;
    grad(x, y, z, gr);
    const gl = Math.hypot(gr[0], gr[1], gr[2]) || 1;
    // f decreases going outward, so the outward normal is −∇f
    normals[p] = -gr[0] / gl; normals[p + 1] = -gr[1] / gl; normals[p + 2] = -gr[2] / gl;
  }

  return { positions, normals, indices: new Uint32Array(idx) };
}

/** Grid half-width for a set of term lists: enough to hold ~99.8% of each part. */
export function gridExtent(set, Z = 1) {
  let L = 0;
  for (const { terms } of set) {
    for (const [, n, ang] of terms) L = Math.max(L, radialExtent(n, ANGULAR[ang].l, 0.998, Z));
  }
  return L;
}

/**
 * Compute every surface for a catalog entry.
 * Returns { id, extent, surfaces: [{ member, sign, positions, normals, indices }],
 *           members: [{ name, c, captured }] }.
 */
export function computeOrbital(orbital, { N = 0, fraction = 0.9, Z = 1 } = {}) {
  const res = N || (orbital.set.length > 1 ? 80 : 96);
  const L = gridExtent(orbital.set, Z);
  const h = 2 * L / (res - 1);
  const surfaces = [];
  const members = [];
  orbital.set.forEach((member, m) => {
    const fn = combination(member.terms, Z);
    const field = sampleGrid(fn, res, L);
    const { c, captured } = contourLevel(field, h, fraction);
    members.push({ name: member.name, c, captured });
    for (const sign of [1, -1]) {
      const surf = extractSurface(field, res, L, c, sign, fn);
      if (surf.indices.length) surfaces.push({ member: m, sign, ...surf });
    }
  });
  return { id: orbital.id, extent: L, N: res, surfaces, members };
}
