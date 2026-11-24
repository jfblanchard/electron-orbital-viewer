# Electron Orbital Viewer

An interactive, textbook-style 3D viewer for **atomic orbitals** (s, p and d),
built as a visual aid for chemistry: shells, nodes, and the shapes behind
bonding.

Semi-transparent probability-density isosurfaces, sign-colored lobes, and free
three.js rotation / zoom / pan let you explore the shape of the wavefunctions
directly, in a way a printed diagram can't.

## What's included

- **Atomic orbitals** — 1s–3s, 2p/3p, all five 3d for hydrogen-like atoms
  (radial parts are hydrogenic; real spherical harmonics give the angular
  shapes).
- **Customizable display** — surface contour (50–95%), opacity, auto-rotate,
  axes, and true relative size.
- **Physics-first math** — wavefunctions are normalized so ∫|ψ|² dV = 1, and
  surfaces enclose exactly the requested fraction of the electron; sign coloring
  follows the wavefunction phase.

## Screenshots

| 2p<sub>z</sub> | 3d<sub>xy</sub> | 3d<sub>z²</sub> |
|:---:|:---:|:---:|
| ![2pz orbital](assets/README-cap-2pz.png) | ![3dxy orbital](assets/README-cap-3dxy.png) | ![3dz2 orbital](assets/README-cap-3dz2.png) |

## Try it live

**<https://orbital-viewer.pages.dev>**

A deluxe edition adds the seven f orbitals, hybrid orbitals (sp to sp³d²), and
hybrid and π orbitals on real molecules (CH₄, NH₃, H₂O, C₂H₄, C₂H₂, CO₂).


## Run locally

```bash
cd deploy
python3 -m http.server 8000
# open http://localhost:8000   (a server is required — ES modules + a Web Worker)
```

To deploy your own copy to Cloudflare Pages:

```bash
cd deploy
npx wrangler pages deploy . --project-name <your-project-name>
```

## Tests

```bash
node tests/test_orbitals.mjs
```

## Layout

```
deploy/   self-contained web app — the deploy unit (see deploy/README.md)
tests/    Node checks for the orbital math and surface meshes
assets/   README screenshots
```

## License / source

Source: <https://github.com/jfblanchard/electron-orbital-viewer>

MIT License (see LICENSE).