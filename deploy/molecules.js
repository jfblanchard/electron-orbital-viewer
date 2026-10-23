/**
 * molecules.js — geometry and bonding for each molecule in the viewer.
 *
 * Experimental bond lengths (Å) and angles. z is "up", matching the orbital
 * viewer; the default camera looks in from roughly +x, so each molecule is
 * laid out across y with its most telling feature (π lobes, lone pairs)
 * pointing up or toward the viewer (an optional `view` overrides the camera
 * direction). buildMolecule() centers everything on the
 * center of mass and derives every hybrid from the actual bond directions.
 */
import { buildMolecule, unit, scale, add, cross } from './molecules-lib.js';

const deg = Math.PI / 180;
const X = [1, 0, 0];
const Y = [0, 1, 0];
const Z = [0, 0, 1];

// Directions at polar angle `polar` from `axis`, spread evenly in azimuth
// starting at +y (so no lobe aims straight at the camera).
const cone = (polar, count, flip = false) =>
  Array.from({ length: count }, (_, i) => {
    const az = (90 + (360 / count) * i) * deg;
    const s = Math.sin(polar * deg);
    const zc = Math.cos(polar * deg);
    return [s * Math.cos(az), s * Math.sin(az), flip ? -zc : zc];
  });

// Two sp² lone pairs on a terminal atom: 120° from the σ bond, in the plane
// perpendicular to that atom's own π orbital.
const sp2LonePairs = (sigmaDir, piDir) => {
  const q = unit(cross(piDir, sigmaDir));
  const back = scale(unit(sigmaDir), -0.5);
  const side = Math.sqrt(3) / 2;
  return [add(back, scale(q, side)), add(back, scale(q, -side))];
};

// --- CH₄ -------------------------------------------------------------------
const CH = 1.09;
const tetra = [Z, ...cone(Math.acos(-1 / 3) / deg, 3)];
const methane = buildMolecule(
  {
    id: 'ch4', group: 'single', formula: 'CH<sub>4</sub>', name: 'Methane',
    hybridization: 'C sp³', shape: 'Tetrahedral', angle: 'H–C–H 109.5°',
    desc: 'Carbon mixes 2s with all three 2p to make four identical sp³ hybrids aimed at the corners of a tetrahedron. Each overlaps a hydrogen 1s: four C–H σ bonds.',
  },
  [['C', [0, 0, 0]], ...tetra.map((d) => ['H', scale(d, CH)])],
  [[0, 1], [0, 2], [0, 3], [0, 4]],
  [{ atom: 0, sigma: [1, 2, 3, 4], lonePairs: 0 }],
);

// --- NH₃ -------------------------------------------------------------------
// Three bonds at polar angle α below the axis meet at 107° when
// cos(107°) = 1.5·cos²α − 0.5.
const NH = 1.01;
const nhPolar = Math.acos(Math.sqrt((Math.cos(107 * deg) + 0.5) / 1.5)) / deg;
const ammonia = buildMolecule(
  {
    id: 'nh3', group: 'single', formula: 'NH<sub>3</sub>', name: 'Ammonia',
    hybridization: 'N sp³', shape: 'Trigonal pyramidal', angle: 'H–N–H 107°',
    desc: 'Nitrogen is sp³, but only three hybrids bond to hydrogen. The fourth holds a lone pair on top, which crowds the bonds down to 107°.',
  },
  [['N', [0, 0, 0]], ...cone(nhPolar, 3, true).map((d) => ['H', scale(d, NH)])],
  [[0, 1], [0, 2], [0, 3]],
  [{ atom: 0, sigma: [1, 2, 3], lonePairs: 1 }],
);

// --- H₂O -------------------------------------------------------------------
// In the yz plane, facing the camera; the two lone pairs point up and out
// toward and away from the viewer ("rabbit ears").
const OH = 0.96;
const half = (104.5 / 2) * deg;
const water = buildMolecule(
  {
    id: 'h2o', group: 'single', formula: 'H<sub>2</sub>O', name: 'Water',
    hybridization: 'O sp³', shape: 'Bent', angle: 'H–O–H 104.5°',
    desc: 'Oxygen is sp³ with two bonding hybrids and two lone pairs. The lone pairs spread out above the molecule and squeeze the H–O–H angle to 104.5°.',
    view: [1, 0.9, 0.45], // oblique, so neither lone pair hides the oxygen
  },
  [
    ['O', [0, 0, 0]],
    ['H', [0, OH * Math.sin(half), -OH * Math.cos(half)]],
    ['H', [0, -OH * Math.sin(half), -OH * Math.cos(half)]],
  ],
  [[0, 1], [0, 2]],
  [{ atom: 0, sigma: [1, 2], lonePairs: 2 }],
);

// --- C₂H₄ ------------------------------------------------------------------
// Flat in the xy plane with C=C along y; the π lobes stand up along z.
const CC2 = 1.34;
const hOut = (180 - 121.5) * deg; // C–H angle off the C=C axis
const ethyleneH = (sy, sx) => [sx * CH * Math.sin(hOut), sy * (CC2 / 2 + CH * Math.cos(hOut)), 0];
const ethylene = buildMolecule(
  {
    id: 'c2h4', group: 'pi', formula: 'C<sub>2</sub>H<sub>4</sub>', name: 'Ethylene',
    hybridization: 'C sp²', shape: 'Trigonal planar at each C', angle: 'H–C–H 117° · C–C–H 121.5°',
    desc: 'Each carbon mixes 2s with two 2p: three sp² hybrids in a plane. One overlaps the other carbon (σ), two reach hydrogens. The unused 2p on each carbon stands straight up; side by side they form the π bond, one lobe above the plane and one below.',
  },
  [
    ['C', [0, CC2 / 2, 0]],
    ['C', [0, -CC2 / 2, 0]],
    ['H', ethyleneH(1, 1)],
    ['H', ethyleneH(1, -1)],
    ['H', ethyleneH(-1, 1)],
    ['H', ethyleneH(-1, -1)],
  ],
  [[0, 1, 2, Z], [0, 2], [0, 3], [1, 4], [1, 5]],
  [
    { atom: 0, sigma: [1, 2, 3], lonePairs: 0, pi: [{ dir: Z, partner: 1 }] },
    { atom: 1, sigma: [0, 4, 5], lonePairs: 0, pi: [{ dir: Z, partner: 0 }] },
  ],
);

// --- C₂H₂ ------------------------------------------------------------------
const CC3 = 1.2;
const CH3 = 1.06;
const acetylene = buildMolecule(
  {
    id: 'c2h2', group: 'pi', formula: 'C<sub>2</sub>H<sub>2</sub>', name: 'Acetylene',
    hybridization: 'C sp', shape: 'Linear', angle: 'H–C–C 180°',
    desc: 'Each carbon is sp: two hybrids 180° apart, one to the other carbon and one to hydrogen. That leaves two 2p orbitals on each carbon, at right angles. They pair up into two π bonds, one vertical and one horizontal, wrapping the C≡C axis.',
  },
  [
    ['C', [0, CC3 / 2, 0]],
    ['C', [0, -CC3 / 2, 0]],
    ['H', [0, CC3 / 2 + CH3, 0]],
    ['H', [0, -CC3 / 2 - CH3, 0]],
  ],
  [[0, 1, 3, Z], [0, 2], [1, 3]],
  [
    { atom: 0, sigma: [1, 2], lonePairs: 0, pi: [{ dir: Z, partner: 1 }, { dir: X, partner: 1 }] },
    { atom: 1, sigma: [0, 3], lonePairs: 0, pi: [{ dir: Z, partner: 0 }, { dir: X, partner: 0 }] },
  ],
);

// --- CO₂ -------------------------------------------------------------------
// Carbon is sp with 2pz and 2px left over. The left C=O π bond uses pz, the
// right one px, so the two π systems are rotated 90° about the O=C=O axis.
// Each oxygen is sp²; its lone pairs lie in the plane perpendicular to its
// own π orbital, so they turn 90° too.
const CO = 1.16;
const co2 = buildMolecule(
  {
    id: 'co2', group: 'pi', formula: 'CO<sub>2</sub>', name: 'Carbon dioxide',
    hybridization: 'C sp · O sp²', shape: 'Linear', angle: 'O–C–O 180°',
    desc: 'Carbon is sp, leaving two 2p orbitals at right angles. Each C=O π bond uses a different one, so the two π systems are turned 90° from each other: one vertical, one horizontal. Each oxygen\'s two lone pairs lie in the plane perpendicular to its own π bond, so they turn with it.',
  },
  [
    ['C', [0, 0, 0]],
    ['O', [0, -CO, 0]],
    ['O', [0, CO, 0]],
  ],
  [[0, 1, 2, Z], [0, 2, 2, X]],
  [
    { atom: 0, sigma: [1, 2], lonePairs: 0, pi: [{ dir: Z, partner: 1 }, { dir: X, partner: 2 }] },
    { atom: 1, sigma: [0], lonePairs: sp2LonePairs(Y, Z), pi: [{ dir: Z, partner: 0 }] },
    { atom: 2, sigma: [0], lonePairs: sp2LonePairs(scale(Y, -1), X), pi: [{ dir: X, partner: 0 }] },
  ],
);

export const MOLECULES = [methane, ammonia, water, ethylene, acetylene, co2];
export const MOLECULES_BY_ID = Object.fromEntries(MOLECULES.map((m) => [m.id, m]));
