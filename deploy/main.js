/**
 * main.js — Three.js scene and UI for the orbital viewer.
 *
 * The physics lives in orbitals.js (wavefunctions → 90% isosurface meshes).
 * worker.js runs that off the main thread; this file only turns the returned
 * vertex/normal/index arrays into translucent meshes and drives the sidebar.
 *
 * Sections
 *   1. Scene       renderer, camera, lights, axes, nucleus
 *   2. Meshes      translucent sign-colored materials, orbital → Group
 *   3. Worker      compute queue with caching and background prefetch
 *   4. UI          sidebar, info card, controls, selection
 */
import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { CSS2DRenderer, CSS2DObject } from 'three/addons/renderers/CSS2DRenderer.js';
import { ORBITALS, ORBITALS_BY_ID } from './orbitals.js';
import { MOLECULES, MOLECULES_BY_ID } from './molecules.js';
import { A0, ELEMENTS, canonical } from './molecules-lib.js';

const A0_ANGSTROM = 0.529177; // Bohr radius in Å
const PLUS_COLOR = '#e8453c';
const MINUS_COLOR = '#3b82f6';
const TRUE_SCALE_RADIUS = 34; // a0; frames the largest orbital (4f) in true-size mode

// Molecular lobes are colored by role (σ / π / lone pair / H 1s), with the
// small ψ < 0 lobe a darker shade; "color by sign" switches to red +/blue −.
const ROLES = {
  sigma: { color: '#38bdf8', label: 'σ hybrid', chip: 'σ' },
  pi: { color: '#4ade80', label: 'π (unhybridized p)', chip: 'π' },
  lp: { color: '#f5a524', label: 'Lone pair', chip: 'lone pairs' },
  h1s: { color: '#e5e7eb', label: 'Hydrogen 1s', chip: 'H 1s' },
};
const ROLE_ORDER = ['sigma', 'pi', 'lp', 'h1s'];

// ============================================================================
// 1. SCENE
// ============================================================================

const viewport = document.getElementById('viewport');

const renderer = new THREE.WebGLRenderer({ antialias: true });
renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
viewport.prepend(renderer.domElement);

const labelRenderer = new CSS2DRenderer({ element: document.getElementById('labels') });

const scene = new THREE.Scene();
scene.background = new THREE.Color(0x111111);

// z is "up", as in chemistry texts. The camera sits off the +x axis so that
// x points out of the screen toward the lower left and y points right.
const camera = new THREE.PerspectiveCamera(35, 1, 0.05, 2000);
camera.up.set(0, 0, 1);
const VIEW_DIR = new THREE.Vector3(1, 0.42, 0.55).normalize();
camera.position.copy(VIEW_DIR).multiplyScalar(40);
scene.add(camera);

const controls = new OrbitControls(camera, renderer.domElement);
controls.enableDamping = true;
controls.dampingFactor = 0.08;
controls.autoRotate = true;
controls.autoRotateSpeed = 0.8;

// Soft fill plus a key light that rides with the camera, so the lit side
// stays consistent however the orbital is turned.
scene.add(new THREE.HemisphereLight(0xffffff, 0x303040, 0.9));
const key = new THREE.DirectionalLight(0xffffff, 2.6);
key.position.set(1, 1.5, 2);
camera.add(key);

// Axes: thin lines through the nucleus with italic x/y/z labels.
const axes = new THREE.Group();
const axisMat = new THREE.LineBasicMaterial({ color: 0x5a606a, transparent: true, opacity: 0.8 });
const axisLabels = [];
for (const [name, dir] of [['x', [1, 0, 0]], ['y', [0, 1, 0]], ['z', [0, 0, 1]]]) {
  const geo = new THREE.BufferGeometry().setFromPoints([
    new THREE.Vector3(...dir).multiplyScalar(-1),
    new THREE.Vector3(...dir),
  ]);
  axes.add(new THREE.Line(geo, axisMat));
  const el = document.createElement('div');
  el.className = 'axis-label';
  el.textContent = name;
  const label = new CSS2DObject(el);
  label.userData.dir = new THREE.Vector3(...dir);
  axes.add(label);
  axisLabels.push(label);
}
scene.add(axes);

function setAxisLength(len) {
  axes.children.forEach((c) => {
    if (c.isLine) c.scale.setScalar(len);
  });
  axisLabels.forEach((l) => l.position.copy(l.userData.dir).multiplyScalar(len * 1.05));
}

// Nucleus: a small opaque point at the origin, sized relative to the view.
const nucleus = new THREE.Mesh(
  new THREE.SphereGeometry(1, 24, 16),
  new THREE.MeshStandardMaterial({ color: 0xffe6a8, emissive: 0x7a5a20, roughness: 0.4 }),
);
scene.add(nucleus);

// CSS2D labels ignore their parent's visibility, so hide them explicitly too.
function showAxes(on) {
  axes.visible = on;
  axisLabels.forEach((l) => { l.visible = on; });
}

// --- molecular framework ----------------------------------------------------
// Ball-and-stick, opaque, drawn before the translucent lobes so the lobes
// tint whatever sits behind them. Balls are kept small so the orbitals read.

const sphereGeo = new THREE.SphereGeometry(1, 32, 20);
const stickGeo = new THREE.CylinderGeometry(1, 1, 1, 16, 1);
const atomMats = Object.fromEntries(Object.entries(ELEMENTS).map(([el, e]) => [
  el, new THREE.MeshStandardMaterial({ color: e.color, roughness: 0.45, metalness: 0.05 }),
]));
const stickMat = new THREE.MeshStandardMaterial({ color: 0x9aa0a6, roughness: 0.5 });
const Y_AXIS = new THREE.Vector3(0, 1, 0);
const Z_AXIS = new THREE.Vector3(0, 0, 1);

function buildFramework(mol) {
  const group = new THREE.Group();
  for (const at of mol.atoms) {
    const m = new THREE.Mesh(sphereGeo, atomMats[at.el]);
    m.position.set(...at.pos);
    m.scale.setScalar(ELEMENTS[at.el].radius);
    group.add(m);
  }
  // Multiple bonds are drawn as parallel sticks spread along the π direction,
  // like the bent bonds of a model kit, so CO₂'s two C=O visibly turn 90°.
  for (const b of mol.bonds) {
    const p = new THREE.Vector3(...mol.atoms[b.a].pos);
    const q = new THREE.Vector3(...mol.atoms[b.b].pos);
    const axis = q.clone().sub(p);
    const len = axis.length();
    axis.normalize();
    const perp = b.perp ? new THREE.Vector3(...b.perp) : new THREE.Vector3(1, 0, 0);
    const offsets = b.order === 1 ? [0] : b.order === 2 ? [-1, 1] : [0, 120, 240];
    const r = b.order === 1 ? 0.07 : 0.05;
    for (const o of offsets) {
      const shift = b.order === 3
        ? perp.clone().applyAxisAngle(axis, THREE.MathUtils.degToRad(o)).multiplyScalar(0.085)
        : perp.clone().multiplyScalar(0.075 * o);
      const s = new THREE.Mesh(stickGeo, stickMat);
      s.position.copy(p).add(q).multiplyScalar(0.5).add(shift);
      s.quaternion.setFromUnitVectors(Y_AXIS, axis);
      s.scale.set(r, len, r);
      group.add(s);
    }
  }
  return group;
}

function resize() {
  const w = viewport.clientWidth, h = viewport.clientHeight;
  renderer.setSize(w, h);
  labelRenderer.setSize(w, h);
  camera.aspect = w / h;
  camera.updateProjectionMatrix();
}
window.addEventListener('resize', resize);
resize();

// Smoothly move the camera so a sphere of radius R fills the view, keeping
// the current viewing direction unless the orbital asks for its own.
let fit = null;
function fitCamera(R, animate = true, view = null) {
  const fov = THREE.MathUtils.degToRad(camera.fov);
  const halfAngle = Math.min(fov, 2 * Math.atan(Math.tan(fov / 2) * camera.aspect)) / 2;
  const dist = (R / Math.sin(halfAngle)) * 1.3;
  const dir = view
    ? new THREE.Vector3(...view).normalize()
    : camera.position.clone().sub(controls.target).normalize();
  const to = dir.multiplyScalar(dist);
  if (!animate) {
    camera.position.copy(to);
    controls.target.set(0, 0, 0);
    fit = null;
    return;
  }
  fit = { from: camera.position.clone(), to, target: controls.target.clone(), t0: performance.now() };
}

function stepFit(now) {
  if (!fit) return;
  const t = Math.min(1, (now - fit.t0) / 450);
  const e = t * t * (3 - 2 * t); // smoothstep
  // interpolate length and direction separately so the camera arcs, not cuts through
  const len = THREE.MathUtils.lerp(fit.from.length(), fit.to.length(), e);
  const dir = fit.from.clone().normalize().lerp(fit.to.clone().normalize(), e).normalize();
  camera.position.copy(dir.multiplyScalar(len));
  controls.target.copy(fit.target).multiplyScalar(1 - e);
  if (t === 1) fit = null;
}

// ============================================================================
// 2. MESHES
// ============================================================================
//
// The "textbook" look: each lobe is a translucent shell that gets denser
// toward its silhouette, like looking through a soap bubble. A rim term
// injected into the standard physical material raises opacity where the
// surface turns edge-on. Back faces draw first, then front faces, so the
// far wall of each lobe shows through the near wall. depthWrite is off so
// overlapping lobes (hybrids, nested 2s/3s shells) blend instead of clipping.

const opacity = { value: 0.45 };
const materials = [];

function makeMaterial(color, side) {
  const mat = new THREE.MeshPhysicalMaterial({
    color,
    side,
    transparent: true,
    opacity: opacity.value,
    depthWrite: false,
    roughness: 0.35,
    metalness: 0,
    clearcoat: 0.4,
    clearcoatRoughness: 0.3,
  });
  mat.onBeforeCompile = (shader) => {
    shader.fragmentShader = shader.fragmentShader.replace(
      '#include <normal_fragment_maps>',
      `#include <normal_fragment_maps>
      float rimF = 1.0 - abs(dot(normalize(vViewPosition), normal));
      diffuseColor.a = clamp(diffuseColor.a * (0.5 + 1.4 * rimF * rimF), 0.0, 1.0);`,
    );
  };
  materials.push(mat);
  return mat;
}

// Materials are cached per color. Hybrid sets can switch from sign colors to
// one hue per hybrid, with the small ψ < 0 lobe a darker shade of that hue.
const HYBRID_COLORS = ['#e8453c', '#f5a524', '#3ecf6e', '#38bdf8', '#a78bfa', '#f472b6'];
const matCache = new Map();
function materialPair(color) {
  if (!matCache.has(color)) {
    matCache.set(color, { back: makeMaterial(color, THREE.BackSide), front: makeMaterial(color, THREE.FrontSide) });
  }
  return matCache.get(color);
}
const darker = (hex) => `#${new THREE.Color(hex).multiplyScalar(0.22).getHexString()}`;

function paint(group, byHybrid) {
  group.children.forEach((m) => {
    const { sign, member, side } = m.userData;
    const hue = HYBRID_COLORS[member % HYBRID_COLORS.length];
    const color = byHybrid ? (sign > 0 ? hue : darker(hue)) : (sign > 0 ? PLUS_COLOR : MINUS_COLOR);
    m.material = materialPair(color)[side];
  });
}

function setOpacity(v) {
  opacity.value = v;
  materials.forEach((m) => { m.opacity = v; });
}

/** Build a Group from a worker result. userData.radius = farthest vertex. */
function buildGroup(res) {
  const group = new THREE.Group();
  let radius = 0;
  for (const s of res.surfaces) {
    const geo = new THREE.BufferGeometry();
    geo.setAttribute('position', new THREE.BufferAttribute(s.positions, 3));
    geo.setAttribute('normal', new THREE.BufferAttribute(s.normals, 3));
    geo.setIndex(new THREE.BufferAttribute(s.indices, 1));
    geo.computeBoundingSphere();
    for (let i = 0; i < s.positions.length; i += 3) {
      radius = Math.max(radius, Math.hypot(s.positions[i], s.positions[i + 1], s.positions[i + 2]));
    }
    const back = new THREE.Mesh(geo);
    const front = new THREE.Mesh(geo);
    back.renderOrder = 1;
    front.renderOrder = 2;
    back.userData = { member: s.member, sign: s.sign, side: 'back' };
    front.userData = { member: s.member, sign: s.sign, side: 'front' };
    group.add(back, front);
  }
  group.userData.radius = radius;
  paint(group, false);
  return group;
}

// ============================================================================
// 3. WORKER — compute queue
// ============================================================================
//
// One orbital at a time, keyed by "id@fraction". A user click jumps the
// queue; idle time prefetches the rest (at 90%) so later clicks are instant.

const worker = new Worker(new URL('./worker.js', import.meta.url), { type: 'module' });
const cache = new Map();   // key → THREE.Group (orbitals) or [{ sign, geometry }] (molecular lobes)
const specs = new Map();   // key → worker request, for molecular lobes (no catalog id)
const waiting = new Map(); // key → [resolve]
const queue = [];
let busy = null;

worker.onmessage = (e) => {
  const res = e.data;
  let out;
  if (specs.has(res.key)) {
    out = res.surfaces.map((s) => {
      const geo = new THREE.BufferGeometry();
      geo.setAttribute('position', new THREE.BufferAttribute(s.positions, 3));
      geo.setAttribute('normal', new THREE.BufferAttribute(s.normals, 3));
      geo.setIndex(new THREE.BufferAttribute(s.indices, 1));
      geo.computeBoundingSphere();
      return { sign: s.sign, geometry: geo };
    });
  } else {
    out = buildGroup(res);
    out.userData.members = res.members;
  }
  cache.set(res.key, out);
  (waiting.get(res.key) || []).forEach((resolve) => resolve(out));
  waiting.delete(res.key);
  busy = null;
  pump();
};

function pump() {
  if (busy) return;
  while (queue.length && cache.has(queue[0])) queue.shift();
  if (!queue.length) return;
  busy = queue.shift();
  const msg = specs.get(busy) || (() => {
    const [id, fraction] = busy.split('@');
    return { key: busy, id, fraction: +fraction };
  })();
  worker.postMessage(msg);
}

const keyOf = (id, fraction) => `${id}@${fraction}`;

function request(key, urgent = false) {
  if (cache.has(key)) return Promise.resolve(cache.get(key));
  const p = new Promise((resolve) => {
    if (!waiting.has(key)) waiting.set(key, []);
    waiting.get(key).push(resolve);
  });
  const i = queue.indexOf(key);
  if (i >= 0) queue.splice(i, 1);
  if (busy !== key) urgent ? queue.unshift(key) : queue.push(key);
  pump();
  return p;
}

function prefetchAll() {
  ORBITALS.forEach((o) => {
    const key = keyOf(o.id, 0.9);
    if (!cache.has(key) && !queue.includes(key) && busy !== key) queue.push(key);
  });
  pump();
}

// --- molecular lobes --------------------------------------------------------
//
// Keys are "<canonical key>@<fraction>" and resolve to a list of
// { sign, geometry } in a0, shared by every lobe of that shape.

function shapeKey(o, fraction) {
  const c = canonical(o);
  const k = `${c.key}@${fraction}`;
  // spheres need fewer cells than directional lobes
  if (!specs.has(k)) specs.set(k, { key: k, Z: c.Z, terms: c.terms, fraction, N: c.n === 1 ? 64 : 80 });
  return k;
}

/** All lobes of a molecule, each a rotated, scaled copy of its canonical shape. */
async function buildOrbitals(mol, fraction) {
  const keys = mol.orbitals.map((o) => shapeKey(o, fraction));
  // urgent requests jump the queue; reverse so the first ends up in front
  const unique = [...new Set(keys)];
  const loaded = new Map();
  await Promise.all(unique.reverse().map(async (k) => loaded.set(k, await request(k, true))));

  const group = new THREE.Group();
  mol.orbitals.forEach((o, idx) => {
    const lobe = new THREE.Group();
    lobe.position.set(...o.center);
    if (o.dir) lobe.quaternion.setFromUnitVectors(Z_AXIS, new THREE.Vector3(...o.dir));
    lobe.scale.setScalar(A0);
    lobe.userData = { role: o.role };
    for (const { sign, geometry } of loaded.get(keys[idx])) {
      const back = new THREE.Mesh(geometry);
      const front = new THREE.Mesh(geometry);
      back.renderOrder = 1;
      front.renderOrder = 2;
      // the small ψ < 0 lobe behind each hybrid (π lobes are equal pairs, not "small")
      const small = sign < 0 && o.role !== 'pi';
      back.userData = { sign, role: o.role, side: 'back', small };
      front.userData = { sign, role: o.role, side: 'front', small };
      lobe.add(back, front);
    }
    group.add(lobe);
  });
  return group;
}

function prefetchMoleculeShapes(fraction) {
  for (const mol of MOLECULES) {
    for (const o of mol.orbitals) {
      const k = shapeKey(o, fraction);
      if (!cache.has(k) && !queue.includes(k) && busy !== k) queue.push(k);
    }
  }
  pump();
}

// ============================================================================
// 4. UI
// ============================================================================

const $ = (id) => document.getElementById(id);
const GROUPS = [
  ['s', 's orbitals', ''],
  ['p', 'p orbitals', ''],
  ['d', 'd orbitals', ''],
  ['f', 'f orbitals', 'wide'],
  ['hybrid', 'Hybrids', 'wide'],
];

// Sidebar
const buttons = new Map();
for (const [key, title, cls] of GROUPS) {
  const section = document.createElement('section');
  section.className = `group ${cls}`;
  section.innerHTML = `<h2>${title}</h2><div class="items"></div>`;
  const items = section.querySelector('.items');
  ORBITALS.filter((o) => o.group === key).forEach((o) => {
    const b = document.createElement('button');
    b.className = 'orb';
    b.innerHTML = o.label;
    b.title = o.desc;
    b.addEventListener('click', () => {
      select(o.id);
      document.body.classList.remove('menu-open');
    });
    items.append(b);
    buttons.set(o.id, b);
  });
  $('list').append(section);
}

// Molecule sections — the six whose hybrid/π orbitals this page can compute.
const MOL_GROUPS = [
  ['single', 'Molecules — σ only'],
  ['pi', 'Molecules — σ + π'],
];
for (const [key, title] of MOL_GROUPS) {
  const section = document.createElement('section');
  section.className = 'group wide';
  section.innerHTML = `<h2>${title}</h2><div class="items"></div>`;
  const items = section.querySelector('.items');
  MOLECULES.filter((m) => m.group === key).forEach((m) => {
    const b = document.createElement('button');
    b.className = 'orb';
    b.innerHTML = `${m.formula} <span class="molname">${m.name}</span>`;
    b.title = m.desc;
    b.addEventListener('click', () => {
      select(m.id);
      document.body.classList.remove('menu-open');
    });
    items.append(b);
    buttons.set(m.id, b);
  });
  $('list').append(section);
}

const plural = (n, word) => `${n} ${word}${n === 1 ? '' : 's'}`;

function renderInfo(o) {
  $('title').innerHTML = o.label;
  $('desc').textContent = o.desc;
  const facts = [];
  if (o.kind === 'atomic') {
    const letters = 'spdf';
    facts.push(`n = ${o.n} · l = ${o.l} (${letters[o.l]}) · ${plural(o.l, 'angular node')} · ${plural(o.n - o.l - 1, 'radial node')}`);
  } else {
    facts.push(`${o.geometry} · ${o.set.length} hybrids from ${o.mix}`);
    facts.push(o.formula);
  }
  facts.push('<span id="size"></span>');
  $('facts').innerHTML = facts.join('<br>');
}

function renderSize(group) {
  const r = group.userData.radius;
  $('size').textContent = `Reaches ${r.toFixed(1)} a₀ (${(r * A0_ANGSTROM).toFixed(1)} Å) from the nucleus`;
}

// Legend: sign colors, or one color per hybrid.
function renderLegend(byHybrid) {
  $('legend-colors').innerHTML = byHybrid
    ? '<div><span class="swatch" style="background: conic-gradient(#e8453c, #f5a524, #3ecf6e, #38bdf8, #a78bfa, #f472b6, #e8453c)"></span>One color per hybrid</div>'
      + '<div><span class="swatch" style="background: #333"></span>Darker shade: ψ &lt; 0</div>'
    : '<div><span class="swatch" style="background: var(--plus)"></span>ψ &gt; 0</div>'
      + '<div><span class="swatch" style="background: var(--minus)"></span>ψ &lt; 0</div>';
}

// Hybrid chips: toggle individual hybrids in a set, and the coloring mode.
let colorByHybrid = false;
function renderChips(o, group) {
  const chips = $('chips');
  chips.innerHTML = '';
  const byHybrid = o.set.length > 1 && colorByHybrid;
  paint(group, byHybrid);
  renderLegend(byHybrid);
  if (o.set.length < 2) return;
  const visible = new Array(o.set.length).fill(true);
  const apply = () => {
    group.children.forEach((m) => { m.visible = visible[m.userData.member]; });
    chips.querySelectorAll('[data-member]').forEach((c) => c.classList.toggle('on', visible[+c.dataset.member]));
  };
  o.set.forEach((_, i) => {
    const c = document.createElement('button');
    c.className = 'chip on';
    c.dataset.member = i;
    c.textContent = `h${String.fromCharCode(0x2081 + i)}`; // h₁, h₂, …
    c.title = 'Show or hide this hybrid';
    c.addEventListener('click', () => { visible[i] = !visible[i]; apply(); });
    chips.append(c);
  });
  const only = document.createElement('button');
  only.className = 'chip';
  only.textContent = 'one / all';
  only.title = 'Show just the first hybrid, or all of them';
  only.addEventListener('click', () => {
    const all = visible.every(Boolean);
    visible.fill(!all);
    visible[0] = true;
    apply();
  });
  chips.append(only);
  const color = document.createElement('button');
  color.className = `chip${colorByHybrid ? ' on' : ''}`;
  color.textContent = 'color by hybrid';
  color.title = 'Give each hybrid its own color instead of coloring by sign';
  color.addEventListener('click', () => {
    colorByHybrid = !colorByHybrid;
    color.classList.toggle('on', colorByHybrid);
    paint(group, colorByHybrid);
    renderLegend(colorByHybrid);
  });
  chips.append(color);
  apply();
}

// --- molecule info & controls ------------------------------------------------

const visibleRoles = { sigma: true, pi: true, lp: true, h1s: false };
let showSmall = true;      // the little ψ < 0 lobe behind each hybrid
let molColorBySign = false;
let shownMol = null;       // Group (framework + lobes) currently in the scene
let molPrefetched = false;

function applyMolVisibility() {
  if (!shownMol) return;
  const lobes = shownMol.userData.lobes;
  lobes.visible = $('orbitals').checked;
  lobes.children.forEach((l) => {
    const role = l.userData.role;
    l.children.forEach((m) => {
      m.visible = visibleRoles[role] && (showSmall || !m.userData.small);
    });
  });
}

function paintMol(group, bySign) {
  group.traverse((m) => {
    if (!m.isMesh) return;
    const { sign, role, side } = m.userData;
    const hue = ROLES[role].color;
    const color = bySign ? (sign > 0 ? PLUS_COLOR : MINUS_COLOR) : (sign > 0 ? hue : darker(hue));
    m.material = materialPair(color)[side];
  });
}

function renderMolInfo(mol) {
  $('title').innerHTML = `${mol.formula} <span class="name">${mol.name}</span>`;
  $('desc').textContent = mol.desc;
  const pi = mol.bonds.reduce((s, b) => s + b.order - 1, 0);
  const lp = mol.orbitals.filter((o) => o.role === 'lp').length;
  $('facts').innerHTML = [
    `${mol.hybridization} · ${mol.shape} · ${mol.angle}`,
    `${plural(mol.bonds.length, 'σ bond')} · ${plural(pi, 'π bond')} · ${plural(lp, 'lone pair')} · ${plural(mol.orbitals.length, 'orbital')}`,
  ].join('<br>');
}

// Role chips: toggle σ / π / lone pairs / H 1s, the small back lobes, and
// sign coloring. Same .chip style as the hybrid chips.
function renderMolChips(mol) {
  const chips = $('chips');
  chips.innerHTML = '';
  const roles = new Set(mol.orbitals.map((o) => o.role));
  const toggle = (label, title, on, fn) => {
    const c = document.createElement('button');
    c.className = `chip${on ? ' on' : ''}`;
    c.textContent = label;
    c.title = title;
    c.addEventListener('click', () => { fn(); c.classList.toggle('on'); });
    chips.append(c);
  };
  for (const role of ROLE_ORDER) {
    if (!roles.has(role)) continue;
    toggle(ROLES[role].chip, `Show or hide ${ROLES[role].label}`, visibleRoles[role], () => {
      visibleRoles[role] = !visibleRoles[role];
      applyMolVisibility();
      renderMolLegend(mol);
    });
  }
  if (mol.orbitals.some((o) => o.role !== 'pi')) {
    toggle('back lobes', 'Show or hide the small ψ < 0 lobe behind each hybrid', showSmall, () => {
      showSmall = !showSmall;
      applyMolVisibility();
    });
  }
  toggle('color by sign', 'Color lobes by wavefunction sign instead of role', molColorBySign, () => {
    molColorBySign = !molColorBySign;
    paintMol(shownMol, molColorBySign);
    renderMolLegend(mol);
  });
}

function renderMolLegend(mol) {
  const roles = new Set(mol.orbitals.map((o) => o.role));
  $('legend-colors').innerHTML = molColorBySign
    ? '<div><span class="swatch" style="background: var(--plus)"></span>ψ &gt; 0</div>'
      + '<div><span class="swatch" style="background: var(--minus)"></span>ψ &lt; 0</div>'
    : ROLE_ORDER.filter((r) => roles.has(r) && visibleRoles[r])
      .map((r) => `<div><span class="swatch" style="background:${ROLES[r].color}"></span>${ROLES[r].label}</div>`)
      .join('');
}

function frameMol(mol, animate = true) {
  // Radius of the framework plus the molecular lobes — they hang well
  // outside the bond ends. Sample vertices of the placed meshes; cheap.
  let R = Math.max(...mol.atoms.map((a) => Math.hypot(...a.pos) + ELEMENTS[a.el].radius));
  const lobes = shownMol.userData.lobes;
  const v = new THREE.Vector3();
  shownMol.updateMatrixWorld(true);
  lobes.children.forEach((lobe) => {
    if (lobe.userData.role === 'h1s') return;
    lobe.children.forEach((m) => {
      const pos = m.geometry.attributes.position;
      for (let i = 0; i < pos.count; i += 7) {
        R = Math.max(R, v.fromBufferAttribute(pos, i).applyMatrix4(m.matrixWorld).length());
      }
    });
  });
  setAxisLength(R);
  fitCamera(R, animate, mol.view || null);
}

async function selectMolecule(id) {
  const mol = MOLECULES_BY_ID[id];
  current = mol.id;
  buttons.forEach((b, k) => b.classList.toggle('active', k === mol.id));
  buttons.get(mol.id).scrollIntoView({ block: 'nearest' });
  if (location.hash.slice(1) !== mol.id) history.replaceState(null, '', `#${mol.id}`);
  renderMolInfo(mol);

  if (shown) { scene.remove(shown); shown = null; }
  if (shownMol) scene.remove(shownMol);
  nucleus.visible = false;
  showAxes(false);

  const outer = new THREE.Group();
  outer.add(buildFramework(mol));

  const loading = setTimeout(() => $('loading').classList.add('show'), 120);
  const fraction = +$('fraction').value;
  const lobes = await buildOrbitals(mol, fraction);
  clearTimeout(loading);
  if (current !== mol.id || fraction !== +$('fraction').value) return;
  $('loading').classList.remove('show');

  paintMol(lobes, molColorBySign);
  outer.userData.lobes = lobes;
  outer.add(lobes);
  scene.add(outer);
  shownMol = outer;
  document.body.dataset.shown = `${mol.id}@${fraction}`;

  renderMolChips(mol);
  applyMolVisibility();
  renderMolLegend(mol);
  $('pct').textContent = `${Math.round(fraction * 100)}%`;

  frameMol(mol, !firstLoad || !keepView);
  keepView = false;
  firstLoad = false;
  if (!molPrefetched) {
    molPrefetched = true;
    prefetchMoleculeShapes(fraction);
  }
}

let current = null;
let shown = null;
let firstLoad = true;
let keepView = false; // set when only the contour level changes

async function select(id) {
  if (MOLECULES_BY_ID[id]) return selectMolecule(id);
  const o = ORBITALS_BY_ID[id] || ORBITALS[0];
  current = o.id;
  buttons.forEach((b, k) => b.classList.toggle('active', k === o.id));
  buttons.get(o.id).scrollIntoView({ block: 'nearest' });
  if (location.hash.slice(1) !== o.id) history.replaceState(null, '', `#${o.id}`);
  renderInfo(o);
  $('chips').innerHTML = '';

  const loading = setTimeout(() => $('loading').classList.add('show'), 120);
  const fraction = +$('fraction').value;
  const group = await request(keyOf(o.id, fraction), true);
  clearTimeout(loading);
  // user moved on while this was computing
  if (current !== o.id || fraction !== +$('fraction').value) return;
  $('loading').classList.remove('show');

  const sameOrbital = shown && shown.userData.id === o.id;
  if (shown) scene.remove(shown);
  if (shownMol) { scene.remove(shownMol); shownMol = null; }
  nucleus.visible = true;
  showAxes($('axes').checked);
  group.children.forEach((m) => { m.visible = true; });
  group.visible = $('orbitals').checked;
  scene.add(group);
  shown = group;
  document.body.dataset.shown = `${o.id}@${fraction}`;
  renderSize(group);
  renderChips(o, group);
  $('pct').textContent = `${Math.round(fraction * 100)}%`;
  group.userData.id = o.id;
  frame(!firstLoad, !(keepView && sameOrbital));
  keepView = false;
  if (firstLoad) {
    firstLoad = false;
    prefetchAll();
  }
}

// Frame the view, size the axes and nucleus for the current orbital.
// Contour changes skip the camera move so the view doesn't jump.
function frame(animate = true, moveCamera = true) {
  if (!shown) return;
  const trueScale = $('truescale').checked;
  const R = trueScale ? TRUE_SCALE_RADIUS : shown.userData.radius;
  setAxisLength(R * 1.12);
  nucleus.scale.setScalar(R * 0.018);
  if (moveCamera) fitCamera(R, animate, ORBITALS_BY_ID[current].view);
}

// Controls
$('opacity').addEventListener('input', (e) => setOpacity(+e.target.value));
$('spin').addEventListener('change', (e) => { controls.autoRotate = e.target.checked; });
$('axes').addEventListener('change', (e) => showAxes(e.target.checked));
$('orbitals').addEventListener('change', (e) => {
  if (shownMol) shownMol.userData.lobes.visible = e.target.checked;
  if (shown) shown.visible = e.target.checked;
});
$('truescale').addEventListener('change', () => frame(true));
$('fraction').addEventListener('change', () => {
  keepView = true;
  select(current);
});
$('menu').addEventListener('click', () => {
  document.body.classList.remove('settings-open');
  document.body.classList.toggle('menu-open');
});
$('settings').addEventListener('click', () => {
  document.body.classList.remove('menu-open');
  document.body.classList.toggle('settings-open');
});
$('settings-close').addEventListener('click', () => document.body.classList.remove('settings-open'));
renderer.domElement.addEventListener('pointerdown', () => {
  document.body.classList.remove('menu-open');
  document.body.classList.remove('settings-open');
});

// Arrow keys step through the whole list, orbitals then molecules.
const NAV = [...ORBITALS.map((o) => o.id), ...MOLECULES.map((m) => m.id)];
window.addEventListener('keydown', (e) => {
  if (['INPUT', 'SELECT'].includes(e.target.tagName) || !['ArrowUp', 'ArrowDown'].includes(e.key)) return;
  e.preventDefault();
  const i = NAV.indexOf(current);
  const next = NAV[(i + (e.key === 'ArrowDown' ? 1 : -1) + NAV.length) % NAV.length];
  select(next);
});
window.addEventListener('hashchange', () => {
  const id = location.hash.slice(1);
  if (id !== current && (ORBITALS_BY_ID[id] || MOLECULES_BY_ID[id])) select(id);
});

// Render loop
renderer.setAnimationLoop((now) => {
  stepFit(now);
  controls.update();
  renderer.render(scene, camera);
  labelRenderer.render(scene, camera);
});

select(location.hash.slice(1) || '2pz');
