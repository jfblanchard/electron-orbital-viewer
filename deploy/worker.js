// Meshes orbitals off the main thread so the UI never stalls.
// Request: { key, id, fraction }   Response: computeOrbital() result plus key, buffers transferred.
import { ORBITALS_BY_ID, computeOrbital } from './orbitals.js';

self.onmessage = (e) => {
  const { key, id, fraction } = e.data;
  const res = computeOrbital(ORBITALS_BY_ID[id], { fraction });
  const transfer = res.surfaces.flatMap((s) => [s.positions.buffer, s.normals.buffer, s.indices.buffer]);
  self.postMessage({ key, ...res }, transfer);
};
