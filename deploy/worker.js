// Meshes orbitals off the main thread so the UI never stalls.
// Request: { key, id, fraction }              — catalog entry (atomic / hybrid set)
//      or: { key, Z, terms, fraction, N }     — canonical molecular lobe along +z
// Response: computeOrbital() result plus key, buffers transferred.
import { ORBITALS_BY_ID, computeOrbital } from './orbitals.js';

self.onmessage = (e) => {
  const { key, id, terms, Z, fraction, N } = e.data;
  const res = id
    ? computeOrbital(ORBITALS_BY_ID[id], { fraction })
    : computeOrbital({ id: key, set: [{ name: key, terms }] }, { N, fraction, Z });
  const transfer = res.surfaces.flatMap((s) => [s.positions.buffer, s.normals.buffer, s.indices.buffer]);
  self.postMessage({ key, ...res }, transfer);
};
