// Numerical checks for deploy/orbitals.js. Run: node tests/test_orbitals.mjs
import {
  radialFn, ANGULAR, ORBITALS, ORBITALS_BY_ID, combination, computeOrbital,
} from '../deploy/orbitals.js';

let failures = 0;
const check = (ok, msg) => {
  console.log(`${ok ? 'PASS' : 'FAIL'}  ${msg}`);
  if (!ok) failures++;
};
const near = (a, b, tol) => Math.abs(a - b) <= tol;

// --- radial normalization: ∫ r² R² dr = 1 --------------------------------
for (const [n, l] of [[1, 0], [2, 0], [2, 1], [3, 0], [3, 1], [3, 2], [4, 3]]) {
  const R = radialFn(n, l);
  let s = 0;
  const dr = 0.002;
  for (let r = dr / 2; r < 200; r += dr) s += r * r * R(r) ** 2 * dr;
  check(near(s, 1, 1e-6), `radial ${n},${l} normalized (${s.toFixed(8)})`);
}

// --- closed forms --------------------------------------------------------
{
  const r = 1.7;
  check(near(radialFn(1, 0)(r), 2 * Math.exp(-r), 1e-12), 'R10 = 2e^-r');
  check(near(radialFn(2, 1)(r), r * Math.exp(-r / 2) / (2 * Math.sqrt(6)), 1e-12), 'R21 closed form');
  check(near(radialFn(3, 2)(r), 4 / (81 * Math.sqrt(30)) * r * r * Math.exp(-r / 3), 1e-12), 'R32 closed form');
  // 2s has one radial node, so the outer-positive phase flips its textbook sign
  const R20 = (r) => -(1 / (2 * Math.SQRT2)) * (2 - r) * Math.exp(-r / 2);
  check(near(radialFn(2, 0)(r), R20(r), 1e-12), 'R20 closed form (outer-positive phase)');
  check(radialFn(3, 1)(20) > 0 && radialFn(3, 0)(25) > 0 && radialFn(2, 0)(10) > 0, 'every R_nl positive at large r');
}

// --- radial nodes: count sign changes of R_nl, expect n − l − 1 ----------
for (const [n, l] of [[2, 0], [3, 0], [3, 1], [3, 2], [4, 3]]) {
  const R = radialFn(n, l);
  let changes = 0, prev = R(1e-3);
  for (let r = 0.01; r < 60; r += 0.01) { const v = R(r); if (v * prev < 0) changes++; prev = v; }
  check(changes === n - l - 1, `${n},${l} has ${changes} radial nodes`);
}

// --- angular orthonormality over the sphere (Gauss-like quadrature) ------
{
  const keys = Object.keys(ANGULAR);
  const nt = 200, np = 400;
  const pts = [];
  for (let a = 0; a < nt; a++) {
    const th = (a + 0.5) * Math.PI / nt;
    for (let b = 0; b < np; b++) {
      const ph = (b + 0.5) * 2 * Math.PI / np;
      pts.push([Math.sin(th) * Math.cos(ph), Math.sin(th) * Math.sin(ph), Math.cos(th),
        Math.sin(th) * (Math.PI / nt) * (2 * Math.PI / np)]);
    }
  }
  let worst = 0;
  for (let i = 0; i < keys.length; i++) {
    for (let j = i; j < keys.length; j++) {
      let s = 0;
      for (const [x, y, z, w] of pts) s += ANGULAR[keys[i]].f(x, y, z) * ANGULAR[keys[j]].f(x, y, z) * w;
      worst = Math.max(worst, Math.abs(s - (i === j ? 1 : 0)));
    }
  }
  check(worst < 1e-3, `real Y_lm orthonormal, worst error ${worst.toExponential(2)}`);
}

// --- hybrid sets are orthonormal (coefficient vectors) -------------------
for (const o of ORBITALS.filter((o) => o.kind === 'hybrid')) {
  const vec = (terms) => Object.fromEntries(terms.map(([c, n, a]) => [`${n}${a}`, c]));
  const vs = o.set.map((m) => vec(m.terms));
  let worst = 0;
  for (let i = 0; i < vs.length; i++) {
    for (let j = 0; j < vs.length; j++) {
      const keys = new Set([...Object.keys(vs[i]), ...Object.keys(vs[j])]);
      let dot = 0;
      for (const k of keys) dot += (vs[i][k] || 0) * (vs[j][k] || 0);
      worst = Math.max(worst, Math.abs(dot - (i === j ? 1 : 0)));
    }
  }
  check(worst < 1e-12, `${o.id} hybrids orthonormal`);
}

// --- surfaces: 90% contour, closed, outward normals, lobe geometry -------
const centroid = (s) => {
  let x = 0, y = 0, z = 0;
  const n = s.positions.length / 3;
  for (let i = 0; i < s.positions.length; i += 3) { x += s.positions[i]; y += s.positions[i + 1]; z += s.positions[i + 2]; }
  return [x / n, y / n, z / n];
};
// connected components of a mesh (by shared vertex indices)
const components = (s) => {
  const nV = s.positions.length / 3;
  const parent = Int32Array.from({ length: nV }, (_, i) => i);
  const find = (a) => { while (parent[a] !== a) { parent[a] = parent[parent[a]]; a = parent[a]; } return a; };
  for (let t = 0; t < s.indices.length; t += 3) {
    const a = find(s.indices[t]), b = find(s.indices[t + 1]), c = find(s.indices[t + 2]);
    parent[b] = a; parent[find(c)] = a;
  }
  const groups = new Map();
  for (let i = 0; i < nV; i++) {
    const r = find(i);
    if (!groups.has(r)) groups.set(r, []);
    groups.get(r).push(i);
  }
  return [...groups.values()].map((vs) => {
    let x = 0, y = 0, z = 0;
    for (const v of vs) { x += s.positions[3 * v]; y += s.positions[3 * v + 1]; z += s.positions[3 * v + 2]; }
    return { size: vs.length, c: [x / vs.length, y / vs.length, z / vs.length] };
  });
};
// fraction of directed edges without exactly one opposite partner (0 = closed, manifold mesh)
const openFraction = (s) => {
  const edges = new Map();
  for (let t = 0; t < s.indices.length; t += 3) {
    for (let e = 0; e < 3; e++) {
      const a = s.indices[t + e], b = s.indices[t + (e + 1) % 3];
      const key = `${a},${b}`;
      edges.set(key, (edges.get(key) || 0) + 1);
    }
  }
  let bad = 0;
  for (const [key, count] of edges) {
    const [a, b] = key.split(',');
    if (count !== 1 || edges.get(`${b},${a}`) !== 1) bad++;
  }
  return bad / edges.size;
};
// winding agrees with gradient normals
const windingAgrees = (s) => {
  let agree = 0, total = 0;
  const P = s.positions, Nn = s.normals;
  for (let t = 0; t < s.indices.length; t += 3) {
    const [a, b, c] = [s.indices[t], s.indices[t + 1], s.indices[t + 2]];
    const e1 = [P[3 * b] - P[3 * a], P[3 * b + 1] - P[3 * a + 1], P[3 * b + 2] - P[3 * a + 2]];
    const e2 = [P[3 * c] - P[3 * a], P[3 * c + 1] - P[3 * a + 1], P[3 * c + 2] - P[3 * a + 2]];
    const n = [e1[1] * e2[2] - e1[2] * e2[1], e1[2] * e2[0] - e1[0] * e2[2], e1[0] * e2[1] - e1[1] * e2[0]];
    const avg = [0, 1, 2].map((d) => Nn[3 * a + d] + Nn[3 * b + d] + Nn[3 * c + d]);
    if (n[0] * avg[0] + n[1] * avg[1] + n[2] * avg[2] > 0) agree++;
    total++;
  }
  return agree / total;
};
const onSurface = (orb, res) => {
  let worst = 0;
  for (const s of res.surfaces) {
    const fn = combination(orb.set[s.member].terms);
    const c = res.members[s.member].c;
    for (let i = 0; i < s.positions.length; i += 3) {
      worst = Math.max(worst, Math.abs(s.sign * fn(s.positions[i], s.positions[i + 1], s.positions[i + 2]) - c) / c);
    }
  }
  return worst;
};
const maxRadius = (res) => {
  let m = 0;
  for (const s of res.surfaces) for (let i = 0; i < s.positions.length; i += 3) {
    m = Math.max(m, Math.hypot(s.positions[i], s.positions[i + 1], s.positions[i + 2]));
  }
  return m;
};
const fmt = (v) => `(${v.map((x) => x.toFixed(2)).join(', ')})`;
const unit = (v) => { const l = Math.hypot(...v); return v.map((x) => x / l); };

const results = {};
for (const orb of ORBITALS) {
  const t0 = performance.now();
  const res = computeOrbital(orb);
  const ms = performance.now() - t0;
  results[orb.id] = res;
  const tris = res.surfaces.reduce((a, s) => a + s.indices.length / 3, 0);
  const open = Math.max(...res.surfaces.map(openFraction));
  const wind = Math.min(...res.surfaces.map(windingAgrees));
  const err = onSurface(orb, res);
  const capt = Math.min(...res.members.map((m) => m.captured));
  const rmax = maxRadius(res);
  const ok = open < 2e-3 && wind > 0.995 && err < 1e-3 && capt > 0.995 && rmax < 0.9 * res.extent * Math.SQRT2;
  check(ok, `${orb.id.padEnd(9)} ${ms.toFixed(0).padStart(5)} ms  tris ${String(tris).padStart(6)}  non-manifold ${(open * 100).toFixed(3)}%  winding ${(wind * 100).toFixed(2)}%  on-surface err ${err.toExponential(1)}  grid holds ${(capt * 100).toFixed(2)}%  r90max ${rmax.toFixed(1)} / L ${res.extent.toFixed(1)}`);
}

// Geometry checks for the orbitals the spec calls out.
{
  // 1s: one + sphere, no − surface; radius ≈ 2.66 a0 (analytic 90% radius of 1s density contour)
  const r = results['1s'];
  const s = r.surfaces;
  const radii = [];
  for (let i = 0; i < s[0].positions.length; i += 3) radii.push(Math.hypot(s[0].positions[i], s[0].positions[i + 1], s[0].positions[i + 2]));
  const mn = Math.min(...radii), mx = Math.max(...radii);
  // For a sphere, the 90% region is r < r90 where 1 − e^(−2r)(1 + 2r + 2r²) = 0.9
  let r90 = 0;
  for (let x = 0; x < 10; x += 1e-5) { if (1 - Math.exp(-2 * x) * (1 + 2 * x + 2 * x * x) >= 0.9) { r90 = x; break; } }
  check(s.length === 1 && s[0].sign === 1 && mx - mn < 0.02 && near(mn, r90, 0.05),
    `1s: single + sphere, radius ${mn.toFixed(3)}–${mx.toFixed(3)} vs analytic ${r90.toFixed(3)}`);
}
{
  // 2s: − inner sphere; + outer shell, whose boundary is two concentric spheres
  const r = results['2s'];
  const pos = r.surfaces.filter((s) => s.sign === 1).flatMap(components);
  const neg = r.surfaces.filter((s) => s.sign === -1).flatMap(components);
  check(neg.length === 1 && pos.length === 2 && [...pos, ...neg].every((c) => Math.hypot(...c.c) < 0.05),
    `2s: ${neg.length} − sphere, ${pos.length} + shell walls, all centered`);
}
{
  const r = results['2pz'];
  const pos = r.surfaces.find((s) => s.sign === 1), neg = r.surfaces.find((s) => s.sign === -1);
  const cp = centroid(pos), cn = centroid(neg);
  check(components(pos).length === 1 && components(neg).length === 1 && cp[2] > 3 && cn[2] < -3 && Math.hypot(cp[0], cp[1]) < 0.05,
    `2pz: + lobe at ${fmt(cp)}, − lobe at ${fmt(cn)}`);
}
{
  const r = results['3dz2'];
  const pos = components(r.surfaces.find((s) => s.sign === 1));
  const neg = components(r.surfaces.find((s) => s.sign === -1));
  const ringOk = neg.length === 1 && Math.hypot(...neg[0].c) < 0.05; // torus centered at the origin
  const lobesOk = pos.length === 2 && pos.every((p) => Math.abs(p.c[2]) > 3 && Math.hypot(p.c[0], p.c[1]) < 0.05);
  check(ringOk && lobesOk, `3dz2: + lobes at ${pos.map((p) => fmt(p.c)).join(' ')}; − torus at ${neg.map((p) => fmt(p.c)).join(' ')}`);
}
{
  const r = results['3dx2y2'];
  const pos = components(r.surfaces.find((s) => s.sign === 1));
  const neg = components(r.surfaces.find((s) => s.sign === -1));
  check(pos.length === 2 && neg.length === 2 && pos.every((p) => Math.abs(p.c[0]) > 3) && neg.every((p) => Math.abs(p.c[1]) > 3),
    `3dx2-y2: + lobes on x, − lobes on y`);
}
{
  // lobe counts for all 4f: fz3 and fxz2/fyz2 have rings; the rest have 8 or 6 lobes
  const expect = { '4fxyz': 8, '4fzx2y2': 8, '4fxx23y2': 6, '4fy3x2y2': 6, '4fz3': 4, '4fxz2': 6, '4fyz2': 6 };
  for (const [id, n] of Object.entries(expect)) {
    const r = results[id];
    const lobes = r.surfaces.flatMap(components).filter((c) => c.size > 20).length;
    check(lobes === n, `${id}: ${lobes} lobes (expected ${n})`);
  }
  // fxyz lobes point toward cube corners: sign(x·y·z) matches lobe sign
  const r = results['4fxyz'];
  const ok = r.surfaces.every((s) => components(s).every(({ c }) => Math.sign(c[0] * c[1] * c[2]) === s.sign && Math.abs(Math.abs(c[0]) - Math.abs(c[2])) < 0.1));
  check(ok, '4fxyz: lobes on cube diagonals with sign = sign(xyz)');
}
{
  // sp3: four big + lobes, pairwise 109.47° apart
  const r = results.sp3;
  const dirs = r.surfaces.filter((s) => s.sign === 1).map((s) => {
    const big = components(s).sort((a, b) => b.size - a.size)[0];
    return unit(big.c);
  });
  const angles = [];
  for (let i = 0; i < dirs.length; i++) for (let j = i + 1; j < dirs.length; j++) {
    angles.push(Math.acos(dirs[i].reduce((a, x, k) => a + x * dirs[j][k], 0)) * 180 / Math.PI);
  }
  check(dirs.length === 4 && angles.every((a) => near(a, 109.47, 0.5)),
    `sp3: big-lobe angles ${angles.map((a) => a.toFixed(2)).join(', ')}°`);
}
// every hybrid's big + lobe points along its intended direction
{
  const t = Math.acos(-1 / 3), a = [90, 210, 330].map((d) => d * Math.PI / 180);
  const trig = a.map((p) => [Math.cos(p), Math.sin(p), 0]);
  const tetra = [[0, 0, 1], ...a.map((p) => [Math.sin(t) * Math.cos(p), Math.sin(t) * Math.sin(p), Math.cos(t)])];
  for (const [id, want] of [
    ['sp', [[1, 0, 0], [-1, 0, 0]]],
    ['sp2', trig],
    ['sp3', tetra],
    ['sp3d', [...trig, [0, 0, 1], [0, 0, -1]]],
    ['sp3d2', [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]]],
  ]) {
    const got = results[id].surfaces.filter((s) => s.sign === 1)
      .map((s) => ({ m: s.member, d: unit(components(s).sort((p, q) => q.size - p.size)[0].c) }));
    const cosines = got.map(({ m, d }) => d.reduce((acc, x, k) => acc + x * want[m][k], 0));
    check(got.length === want.length && cosines.every((c) => c > 0.999),
      `${id}: big + lobes along intended directions (min cos ${Math.min(...cosines).toFixed(4)})`);
  }
}

console.log(failures ? `\n${failures} FAILED` : '\nall passed');
process.exit(failures ? 1 : 0);
