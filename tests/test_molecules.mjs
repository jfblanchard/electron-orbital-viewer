// Geometry and orbital checks for the molecule viewer.
// Run from deploy/:  node --test test_molecules.mjs
import { test } from 'node:test';
import assert from 'node:assert/strict';
import { MOLECULES_BY_ID as M } from '../deploy/molecules.js';
import {
  ELEMENTS, sub, dot, norm, unit, cross, angleDeg, canonical, rotationZTo, apply3,
  hybridAlong, EQUIVALENT_A,
} from '../deploy/molecules-lib.js';
import { combination, computeOrbital } from '../deploy/orbitals.js';

const near = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol, `${msg}: ${a} vs ${b} (±${tol})`);
const pos = (mol, i) => mol.atoms[i].pos;
const dist = (mol, i, j) => norm(sub(pos(mol, i), pos(mol, j)));
const angle = (mol, i, c, j) => angleDeg(sub(pos(mol, i), pos(mol, c)), sub(pos(mol, j), pos(mol, c)));
const on = (mol, atom, role) => mol.orbitals.filter((o) => o.atom === atom && (!role || o.role === role));
const TETRA = (Math.acos(-1 / 3) * 180) / Math.PI;

// Deterministic pseudo-random points for field comparisons.
function* points(count, R, seed = 7) {
  let s = seed;
  const rnd = () => ((s = (s * 16807) % 2147483647) / 2147483647) * 2 - 1;
  for (let i = 0; i < count; i++) yield [rnd() * R, rnd() * R, rnd() * R];
}

// --- geometry -----------------------------------------------------------

test('bond lengths match experiment', () => {
  const cases = [
    ['ch4', [[0, 1, 1.09], [0, 2, 1.09], [0, 3, 1.09], [0, 4, 1.09]]],
    ['nh3', [[0, 1, 1.01], [0, 2, 1.01], [0, 3, 1.01]]],
    ['h2o', [[0, 1, 0.96], [0, 2, 0.96]]],
    ['c2h4', [[0, 1, 1.34], [0, 2, 1.09], [0, 3, 1.09], [1, 4, 1.09], [1, 5, 1.09]]],
    ['c2h2', [[0, 1, 1.2], [0, 2, 1.06], [1, 3, 1.06]]],
    ['co2', [[0, 1, 1.16], [0, 2, 1.16]]],
  ];
  for (const [id, bonds] of cases) {
    for (const [i, j, d] of bonds) near(dist(M[id], i, j), d, 1e-9, `${id} ${i}-${j}`);
  }
});

test('bond angles match experiment', () => {
  const ch4 = M.ch4;
  for (let i = 1; i <= 4; i++) for (let j = i + 1; j <= 4; j++) near(angle(ch4, i, 0, j), TETRA, 1e-9, `CH4 H${i}-C-H${j}`);

  const nh3 = M.nh3;
  for (const [i, j] of [[1, 2], [1, 3], [2, 3]]) near(angle(nh3, i, 0, j), 107, 1e-9, `NH3 H${i}-N-H${j}`);

  near(angle(M.h2o, 1, 0, 2), 104.5, 1e-9, 'H2O H-O-H');

  const e = M.c2h4;
  near(angle(e, 2, 0, 3), 117, 1e-9, 'C2H4 H-C-H (C1)');
  near(angle(e, 4, 1, 5), 117, 1e-9, 'C2H4 H-C-H (C2)');
  for (const [h, c, o] of [[2, 0, 1], [3, 0, 1], [4, 1, 0], [5, 1, 0]]) near(angle(e, h, c, o), 121.5, 1e-9, `C2H4 C-C-H${h}`);

  near(angle(M.co2, 1, 0, 2), 180, 1e-6, 'CO2 O-C-O');
  near(angle(M.c2h2, 2, 0, 1), 180, 1e-6, 'C2H2 H-C-C');
  near(angle(M.c2h2, 3, 1, 0), 180, 1e-6, 'C2H2 C-C-H');
});

test('ethylene is planar and every molecule sits on its center of mass', () => {
  for (const at of M.c2h4.atoms) near(at.pos[2], 0, 1e-12, 'C2H4 atom off the xy plane');
  for (const mol of Object.values(M)) {
    let com = [0, 0, 0], total = 0;
    for (const at of mol.atoms) {
      const m = ELEMENTS[at.el].mass;
      total += m;
      com = com.map((c, k) => c + m * at.pos[k]);
    }
    for (const c of com) near(c / total, 0, 1e-12, `${mol.id} center of mass`);
  }
});

// --- hybrids ------------------------------------------------------------

test('each atom\'s hybrids and π orbitals are orthonormal', () => {
  for (const mol of Object.values(M)) {
    const atoms = new Set(mol.orbitals.filter((o) => o.role !== 'h1s').map((o) => o.atom));
    for (const a of atoms) {
      const set = on(mol, a);
      assert.ok(set.length <= 4, `${mol.id} atom ${a} uses at most 2s + three 2p`);
      for (let i = 0; i < set.length; i++) {
        for (let j = i; j < set.length; j++) {
          const d = set[i].vec.reduce((s, x, k) => s + x * set[j].vec[k], 0);
          near(d, i === j ? 1 : 0, 1e-9, `${mol.id} atom ${a} <h${i}|h${j}>`);
        }
      }
    }
  }
});

test('valence s is fully used on every sp³ and sp² center', () => {
  const sTotal = (mol, a) => on(mol, a).reduce((s, o) => s + o.a * o.a, 0);
  for (const [id, a] of [['ch4', 0], ['nh3', 0], ['h2o', 0], ['c2h4', 0], ['c2h4', 1], ['c2h2', 0], ['co2', 0], ['co2', 1], ['co2', 2]]) {
    near(sTotal(M[id], a), 1, 1e-9, `${id} atom ${a} s-character`);
  }
});

test('s-character follows the geometry (Coulson)', () => {
  for (const o of on(M.ch4, 0)) near(o.a, EQUIVALENT_A.sp3, 1e-12, 'CH4 sp3');
  for (const o of on(M.c2h2, 0, 'sigma')) near(o.a, EQUIVALENT_A.sp, 1e-12, 'C2H2 sp');
  for (const o of on(M.co2, 0, 'sigma')) near(o.a, EQUIVALENT_A.sp, 1e-12, 'CO2 carbon sp');
  for (const o of on(M.co2, 1).filter((x) => x.role !== 'pi')) near(o.a, EQUIVALENT_A.sp2, 1e-12, 'CO2 oxygen sp2');
  // water: bonds cos θ/(cos θ − 1) = 20.0% s, lone pairs take the rest
  const c = Math.cos((104.5 * Math.PI) / 180);
  for (const o of on(M.h2o, 0, 'sigma')) near(o.a * o.a, c / (c - 1), 1e-12, 'H2O bond s-fraction');
  for (const o of on(M.h2o, 0, 'lp')) near(o.a * o.a, (1 - 2 * c / (c - 1)) / 2, 1e-12, 'H2O lone pair s-fraction');
  // ethylene: C–C hybrid has more s than C–H (Bent's rule)
  const [cc] = on(M.c2h4, 0, 'sigma').filter((o) => o.partner === 1);
  const ch = on(M.c2h4, 0, 'sigma').filter((o) => o.partner !== 1);
  assert.ok(cc.a > ch[0].a, 'C2H4 C–C hybrid richer in s than C–H');
});

test('σ hybrid directions match bond directions', () => {
  for (const mol of Object.values(M)) {
    for (const o of mol.orbitals.filter((x) => x.role === 'sigma')) {
      const bond = unit(sub(pos(mol, o.partner), o.center));
      near(angleDeg(o.dir, bond), 0, 1e-6, `${mol.id} atom ${o.atom}→${o.partner}`);
    }
  }
});

test('σ hybrid wavefunction actually peaks along its bond', () => {
  // scan a sphere of directions for the largest ψ, a few a0/Z out
  const dirs = [];
  const K = 20000;
  const golden = Math.PI * (3 - Math.sqrt(5));
  for (let i = 0; i < K; i++) {
    const z = 1 - (2 * (i + 0.5)) / K;
    const r = Math.sqrt(1 - z * z);
    dirs.push([r * Math.cos(golden * i), r * Math.sin(golden * i), z]);
  }
  for (const id of ['ch4', 'c2h4', 'co2', 'h2o']) {
    const mol = M[id];
    for (const o of mol.orbitals.filter((x) => x.role === 'sigma')) {
      const f = combination(o.terms, o.Z);
      const r = 3 / o.Z;
      let best = null, bestV = -Infinity;
      for (const d of dirs) {
        const v = f(d[0] * r, d[1] * r, d[2] * r);
        if (v > bestV) { bestV = v; best = d; }
      }
      const bond = sub(pos(mol, o.partner), o.center);
      assert.ok(angleDeg(best, bond) < 2, `${id} ${o.atom}→${o.partner} peaks ${angleDeg(best, bond).toFixed(2)}° off the bond`);
    }
  }
});

test('hybridAlong builds the textbook sp³ hybrid', () => {
  const t = hybridAlong([0, 0, 1]);
  assert.deepEqual(t.map(([, n, k]) => [n, k]), [[2, 's'], [2, 'pz']]);
  near(t[0][0], 0.5, 1e-15, 'a');
  near(t[1][0], Math.sqrt(3) / 2, 1e-15, 'b');
});

test('each lobe equals its canonical +z shape, rotated onto the bond', () => {
  for (const mol of Object.values(M)) {
    for (const o of mol.orbitals) {
      const c = canonical(o);
      const world = combination(o.terms, o.Z);
      const local = combination(c.terms, c.Z);
      const R = rotationZTo(o.dir || [0, 0, 1]);
      for (const p of points(40, 4)) {
        const q = apply3(R, p);
        near(world(...q), local(...p), 1e-10, `${mol.id} ${o.role} on atom ${o.atom}`);
      }
    }
  }
});

// --- lone pairs ---------------------------------------------------------

test('NH3 lone pair points straight up the three-fold axis', () => {
  const [lp] = on(M.nh3, 0, 'lp');
  const sumBonds = on(M.nh3, 0, 'sigma').reduce((s, o) => s.map((x, k) => x + o.dir[k]), [0, 0, 0]);
  near(angleDeg(lp.dir, sumBonds), 180, 1e-6, 'lone pair opposite the bonds');
});

test('H2O lone pairs are mirror images across the molecular plane', () => {
  const nrm = unit(cross(sub(pos(M.h2o, 1), pos(M.h2o, 0)), sub(pos(M.h2o, 2), pos(M.h2o, 0))));
  const [a, b] = on(M.h2o, 0, 'lp');
  near(dot(a.dir, nrm), -dot(b.dir, nrm), 1e-12, 'opposite sides');
  assert.ok(Math.abs(dot(a.dir, nrm)) > 0.5, 'clearly out of plane');
  near(a.a, b.a, 1e-12, 'equal s-character');
  // the lone-pair angle opens wider than the bond angle
  assert.ok(angleDeg(a.dir, b.dir) > 104.5, `lone pairs ${angleDeg(a.dir, b.dir).toFixed(1)}° apart`);
});

// --- π systems ----------------------------------------------------------

test('C2H4 π orbitals stand perpendicular to the molecular plane', () => {
  for (const c of [0, 1]) {
    const [pi] = on(M.c2h4, c, 'pi');
    near(Math.abs(pi.dir[2]), 1, 1e-12, `C${c} π along the plane normal`);
    for (const s of on(M.c2h4, c, 'sigma')) near(s.dir[2], 0, 1e-12, 'σ hybrid in plane');
  }
  const [p0] = on(M.c2h4, 0, 'pi');
  const [p1] = on(M.c2h4, 1, 'pi');
  near(dot(p0.dir, p1.dir), 1, 1e-12, 'both p orbitals in phase (bonding π)');
});

test('C2H4 π lobes are symmetric above and below the plane', () => {
  for (const c of [0, 1]) {
    const [pi] = on(M.c2h4, c, 'pi');
    const f = combination(pi.terms, pi.Z);
    for (const [x, y, z] of points(50, 3)) near(f(x, y, z), -f(x, y, -z), 1e-12, 'ψ(z) = −ψ(−z)');

    const res = computeOrbital({ set: [{ name: 'pi', terms: pi.terms }] }, { N: 48, fraction: 0.9, Z: pi.Z });
    const centroid = (s) => {
      let z = 0;
      for (let i = 2; i < s.positions.length; i += 3) z += s.positions[i];
      return z / (s.positions.length / 3);
    };
    const up = res.surfaces.find((s) => s.sign > 0);
    const down = res.surfaces.find((s) => s.sign < 0);
    assert.ok(centroid(up) > 0.5 && centroid(down) < -0.5, 'one lobe above, one below');
    near(centroid(up), -centroid(down), 1e-3, 'lobes mirror each other');
    near(up.positions.length, down.positions.length, 0, 'lobes have identical meshes');
  }
});

test('CO2: the two π systems are orthogonal', () => {
  const co2 = M.co2;
  const axis = unit(sub(pos(co2, 2), pos(co2, 1)));
  const cPi = on(co2, 0, 'pi');
  assert.equal(cPi.length, 2, 'carbon keeps two p orbitals');
  near(dot(cPi[0].dir, cPi[1].dir), 0, 1e-12, 'carbon p orbitals at 90°');
  for (const p of cPi) near(dot(p.dir, axis), 0, 1e-12, 'carbon p ⊥ O=C=O axis');

  const [o1] = on(co2, 1, 'pi');
  const [o2] = on(co2, 2, 'pi');
  near(dot(o1.dir, o2.dir), 0, 1e-12, 'oxygen p orbitals rotated 90° from each other');
  // each oxygen's p lines up with the carbon p it bonds to
  for (const o of [o1, o2]) {
    const partner = cPi.find((p) => p.partner === o.atom);
    near(Math.abs(dot(o.dir, partner.dir)), 1, 1e-12, `O${o.atom} p parallel to its carbon partner`);
  }
});

test('CO2: each oxygen\'s lone pairs lie in the plane perpendicular to its own π', () => {
  const co2 = M.co2;
  const lpNormal = (a) => {
    const [l1, l2] = on(co2, a, 'lp');
    return unit(cross(l1.dir, l2.dir));
  };
  for (const a of [1, 2]) {
    const [pi] = on(co2, a, 'pi');
    for (const lp of on(co2, a, 'lp')) near(dot(lp.dir, pi.dir), 0, 1e-12, `O${a} lone pair ⊥ its π`);
    for (const lp of on(co2, a, 'lp')) near(angleDeg(lp.dir, on(co2, a, 'sigma')[0].dir), 120, 1e-9, `O${a} lone pair 120° from σ`);
  }
  near(dot(lpNormal(1), lpNormal(2)), 0, 1e-12, 'the two lone-pair planes are perpendicular');
});

test('C2H2: two orthogonal π bonds wrap the axis', () => {
  const axis = unit(sub(pos(M.c2h2, 0), pos(M.c2h2, 1)));
  for (const c of [0, 1]) {
    const pis = on(M.c2h2, c, 'pi');
    assert.equal(pis.length, 2);
    near(dot(pis[0].dir, pis[1].dir), 0, 1e-12, 'π orbitals at 90°');
    for (const p of pis) near(dot(p.dir, axis), 0, 1e-12, 'π ⊥ axis');
  }
});

test('bond orders: σ/π counts', () => {
  const count = (mol) => ({ sigma: mol.bonds.length, pi: mol.bonds.reduce((s, b) => s + b.order - 1, 0) });
  assert.deepEqual(count(M.ch4), { sigma: 4, pi: 0 });
  assert.deepEqual(count(M.c2h4), { sigma: 5, pi: 1 });
  assert.deepEqual(count(M.c2h2), { sigma: 3, pi: 2 });
  assert.deepEqual(count(M.co2), { sigma: 2, pi: 2 });
  // every π bond has one p orbital on each end
  for (const mol of Object.values(M)) {
    const pOrbitals = mol.orbitals.filter((o) => o.role === 'pi').length;
    assert.equal(pOrbitals, 2 * count(mol).pi, `${mol.id} p orbitals per π bond`);
  }
});

// --- meshes -------------------------------------------------------------

test('canonical lobes mesh cleanly, big lobe forward', () => {
  const seen = new Map();
  for (const mol of Object.values(M)) for (const o of mol.orbitals) seen.set(canonical(o).key, canonical(o));
  assert.ok(seen.size <= 14, `distinct shapes stay few (${seen.size})`);
  for (const c of seen.values()) {
    const res = computeOrbital({ set: [{ name: c.key, terms: c.terms }] }, { N: 48, fraction: 0.5, Z: c.Z });
    assert.ok(res.surfaces.length >= 1, `${c.key} has a surface`);
    for (const s of res.surfaces) {
      assert.ok(s.indices.length > 0 && s.indices.length % 3 === 0, `${c.key} triangles`);
      assert.ok(s.positions.every(Number.isFinite), `${c.key} finite vertices`);
    }
    if (c.terms.length === 2) {
      // hybrid: the + lobe reaches farther along +z than anything reaches along −z
      let zmax = -Infinity, zmin = Infinity;
      for (const s of res.surfaces) {
        for (let i = 2; i < s.positions.length; i += 3) { zmax = Math.max(zmax, s.positions[i]); zmin = Math.min(zmin, s.positions[i]); }
      }
      assert.ok(zmax > -zmin, `${c.key} big lobe forward (${zmax.toFixed(2)} vs ${zmin.toFixed(2)})`);
    }
  }
});
