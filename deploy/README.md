# Electron Orbital Viewer — deploy unit

Static site: `index.html`, `main.js`, `orbitals.js`, `worker.js`. Three.js
0.160 loads from unpkg via an import map. No build step, no data files.

## Run locally

```bash
cd deploy
python3 -m http.server 8000
# open http://localhost:8000
```

It needs a server (not `file://`) because it uses ES modules and a module
Web Worker. Deep links work: `http://localhost:8000/#3dz2`.

## Deploy

```bash
cd deploy
wrangler pages deploy . --project-name electron-orbital-viewer
```

## Files

| File | Role |
|------|------|
| `orbitals.js` | The physics. Radial functions, real spherical harmonics, hybrids, the orbital catalog, and 90% isosurface meshing. No Three.js dependency, so it runs in the browser, the worker, and Node. |
| `worker.js` | Runs `computeOrbital()` off the main thread. |
| `main.js` | Three.js scene, translucent materials, compute queue with prefetch, sidebar and controls. |
| `index.html` | Layout and styles. |

## Adding an orbital

Add an entry to `ORBITALS` in `orbitals.js`. An atomic orbital is one line:

```js
atomic('4dz2', 'd', 4, 'dz2', '4d<sub>z²</sub>', 'Like 3dz² with one spherical node'),
```

A hybrid is any orthonormal set of linear combinations, given as
`[coefficient, n, angularKey]` terms. It shows up in the sidebar automatically.
Run `node tests/test_orbitals.mjs` from the project root afterwards.
