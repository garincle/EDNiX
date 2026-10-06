import * as THREE from 'three';
import { OrbitControls } from './vendor/OrbitControls.js';

// ================================================================== API & helpers
const TOKEN = (() => {
  const m = location.hash.match(/token=([^&]+)/);
  if (m) { sessionStorage.setItem('ednixToken', m[1]); history.replaceState(null, '', location.pathname); return m[1]; }
  return sessionStorage.getItem('ednixToken') || '';
})();
async function api(path, body) {
  const opt = { headers: { 'X-EDNiX-Token': TOKEN } };
  if (body !== undefined) { opt.method = 'POST'; opt.headers['Content-Type'] = 'application/json'; opt.body = JSON.stringify(body); }
  const r = await fetch(path, opt);
  const j = await r.json().catch(() => ({ error: `invalid response (${r.status})` }));
  if (!r.ok || j.error) throw new Error(j.error || `error ${r.status}`);
  return j;
}
const $ = (id) => document.getElementById(id);
const esc = (s) => String(s).replace(/[&<>"]/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
function setOptions(sel, items, keep) {
  const prev = keep ?? sel.value;
  sel.innerHTML = items.map(([v, l]) => `<option value="${esc(v)}">${esc(l)}</option>`).join('');
  if (items.some(([v]) => String(v) === String(prev))) sel.value = prev;
}
function banner(text, err = false) { const b = $('banner'); b.hidden = !text; b.textContent = text || ''; b.classList.toggle('err', err); }
const plural = (n, w) => `${n} ${w}${n === 1 ? '' : 's'}`;

// ================================================================== navigation
let currentView = 'atlas';
let currentPipe = 'func';
document.querySelectorAll('.nav').forEach((b) => b.addEventListener('click', () => showView(b.dataset.view, b.dataset.pipe)));
function showView(v, pipe) {
  if (v === 'launcher' && pipe) switchPipe(pipe);
  currentView = v;
  document.querySelectorAll('.nav').forEach((b) => b.classList.toggle('on', b.dataset.view === v && (v !== 'launcher' || b.dataset.pipe === currentPipe)));
  document.querySelectorAll('.view').forEach((s) => s.classList.toggle('on', s.id === `view-${v}`));
  if (v === 'atlas') onResize();
  if (v === 'runs') refreshRuns();
  if (v === 'qc') loadQcCandidates();
}

// ================================================================== startup
let INFO = null;
async function init() {
  if (!TOKEN) { banner('Open the interface with the EDNiX launcher (double-click): its address contains an access key.', true); return; }
  try { INFO = await api('/api/info'); } catch (e) { banner(`EDNiX server unreachable: ${e.message}`, true); return; }
  $('env').textContent = `Python ${INFO.python_version}`;
  const probs = [];
  if (!INFO.repo) probs.push('EDNiX repository not found');
  if (INFO.missing.length) probs.push(`missing Python modules: ${INFO.missing.join(', ')}`);
  if (!INFO.atlas_library_ok) probs.push(`atlas library not found (${INFO.atlas_library})`);
  if (probs.length) banner(probs.join('; '), !INFO.repo);
  switchPipe('func');
  await loadSpecies();
  refreshRuns();
  setInterval(() => { if (currentView === 'runs') refreshRuns(); else refreshBadge(); }, 4000);
}

// ================================================================== 3D scene
const stage = $('stage');
const renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true });
renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
stage.prepend(renderer.domElement);
const scene = new THREE.Scene();
const camera = new THREE.PerspectiveCamera(35, 1, 0.1, 10000);
const controls = new OrbitControls(camera, renderer.domElement);
controls.enableDamping = true;
scene.add(new THREE.HemisphereLight(0xffffff, 0x8a9a95, 1.1));
const key = new THREE.DirectionalLight(0xffffff, 1.4); key.position.set(1, 2, 1.5); camera.add(key); scene.add(camera);
const fill = new THREE.DirectionalLight(0xffffff, 0.5); fill.position.set(-2, -1, -1); scene.add(fill);
// RAS millimetres -> three.js (x = right, y = superior, -z = anterior); user orientation fixes are applied here too
const world = new THREE.Group(); scene.add(world);
const BASE_Q = new THREE.Quaternion().setFromAxisAngle(new THREE.Vector3(1, 0, 0), -Math.PI / 2);
let userQ = new THREE.Quaternion();
function applyWorldRotation() { world.quaternion.copy(userQ).multiply(BASE_Q); }
applyWorldRotation();
const EOSIN = new THREE.Color('#C8466F');
function onResize() {
  const w = stage.clientWidth, h = stage.clientHeight;
  if (!w || !h) return;
  renderer.setSize(w, h); camera.aspect = w / h; camera.updateProjectionMatrix();
}
window.addEventListener('resize', onResize);
(function loop() { requestAnimationFrame(loop); controls.update(); renderer.render(scene, camera); })();

function b64ToArray(b64, Type) {
  const bin = atob(b64); const buf = new Uint8Array(bin.length);
  for (let i = 0; i < bin.length; i++) buf[i] = bin.charCodeAt(i);
  return new Type(buf.buffer);
}
function hashColor(id) { return new THREE.Color(`hsl(${(id * 137.508) % 360}, 55%, 58%)`); }

// ================================================================== layers
const LAYERS = [];          // {id, kind, label, meta, group, inner, regions: Map, selected: Set, visible, opacity, center, radius, ...}
let ACTIVE = null;
let LAYOUT = 'overlay';
let nextLayerId = 1;
const layerById = (id) => LAYERS.find((l) => l.id === id);

async function loadSpecies() {
  const sp = await api('/api/atlas/species');
  setOptions($('sel-species'), sp.map((s) => [s, s]));
  if (sp.length) await loadRefs();
}
async function loadRefs() {
  const refs = await api(`/api/atlas/references?species=${encodeURIComponent($('sel-species').value)}`);
  setOptions($('sel-ref'), refs.map((r) => [r, r]));
  if (refs.length) await loadAtlasList();
}
let CATALOG = { atlases: [], templates: [], surfaces: { surfaces: [], atlases: [] } };
async function loadAtlasList() {
  CATALOG = await api(`/api/atlas/list?species=${encodeURIComponent($('sel-species').value)}&reference=${encodeURIComponent($('sel-ref').value)}`);
  const hasSurf = CATALOG.surfaces.surfaces.length && CATALOG.surfaces.atlases.length;
  $('sel-source').querySelector('option[value=surface]').disabled = !hasSurf;
  if (!hasSurf) $('sel-source').value = 'volume';
  setOptions($('sel-template'), [['', 'none'], ...CATALOG.templates.map((t) => [t, t])], CATALOG.templates[0] || '');
  setOptions($('sel-surface'), CATALOG.surfaces.surfaces.map((s) => [s, s.replace('_', ' ')]), 'pial');
  onSource();
}
function onSource() {
  const surf = $('sel-source').value === 'surface';
  document.querySelectorAll('.vol-only').forEach((e) => { e.hidden = surf; });
  document.querySelectorAll('.surf-only').forEach((e) => { e.hidden = !surf; });
  const list = surf ? CATALOG.surfaces.atlases.map((a) => [a, a])
    : CATALOG.atlases.map((a) => [a.key, a.has_labels ? a.key : `${a.key} (no region names)`]);
  setOptions($('sel-atlas'), list);
  updateLevels();
}
function updateLevels() {
  const a = CATALOG.atlases.find((x) => x.key === $('sel-atlas').value);
  const n = a ? a.levels : 1;
  const names = ['coarse', 'intermediate', 'fine', 'finest'];
  setOptions($('sel-level'), Array.from({ length: n }, (_, i) => [i, n > 1 ? `${i} (${names[Math.min(i, 3)]})` : '0']), n - 1);
}
$('sel-species').addEventListener('change', loadRefs);
$('sel-ref').addEventListener('change', loadAtlasList);
$('sel-source').addEventListener('change', onSource);
$('sel-atlas').addEventListener('change', updateLevels);
$('btn-add').addEventListener('click', addLayer);

async function addLayer() {
  const surf = $('sel-source').value === 'surface';
  const meta = { species: $('sel-species').value, reference: $('sel-ref').value, atlas: $('sel-atlas').value,
                 level: surf ? null : Number($('sel-level').value), surface: surf ? $('sel-surface').value : null,
                 template: surf ? null : $('sel-template').value };
  if (!meta.atlas) { banner('No atlas available for this species and reference.'); return; }
  $('stage-empty').hidden = true; $('stage-loading').hidden = false;
  $('loading-text').textContent = surf ? 'Loading the cortical surfaces' : 'Building the region surfaces';
  let data;
  try {
    const q = new URLSearchParams(Object.fromEntries(Object.entries(meta).filter(([, v]) => v !== null && v !== '')));
    data = await api(surf ? `/api/atlas/surface?${q}` : `/api/atlas/mesh?${q}`);
  } catch (e) {
    $('stage-loading').hidden = true; if (!LAYERS.length) $('stage-empty').hidden = false;
    banner(`Could not display this atlas: ${e.message}`, true); return;
  }
  $('stage-loading').hidden = true;
  banner(data.label_warning || '');
  stopQuiz();
  const layer = surf ? buildSurfaceLayer(data, meta) : buildVolumeLayer(data, meta);
  LAYERS.push(layer); world.add(layer.group);
  setActive(layer.id); arrange(); fitView('left'); renderLayers(); onResize();
}

function newLayer(kind, meta, label) {
  const group = new THREE.Group(); const inner = new THREE.Group(); group.add(inner);
  return { id: nextLayerId++, kind, meta, label, group, inner, regions: new Map(), selected: new Set(), visible: true,
           opacity: Number($('rng-opacity').value), explode: 0, center: new THREE.Vector3(), radius: 10, shell: null, hover: null };
}

function buildVolumeLayer(data, meta) {
  const L = newLayer('volume', meta, `${meta.species} ${data.atlas_name} level ${meta.level}`);
  L.meta.atlasName = data.atlas_name;
  const box = new THREE.Box3();
  for (const r of data.regions) {
    const g = new THREE.BufferGeometry();
    g.setAttribute('position', new THREE.BufferAttribute(b64ToArray(r.v, Float32Array), 3));
    g.setIndex(new THREE.BufferAttribute(b64ToArray(r.f, Uint32Array), 1));
    g.computeVertexNormals(); g.computeBoundingBox(); box.union(g.boundingBox);
    const color = r.rgb ? new THREE.Color(`rgb(${r.rgb[0]},${r.rgb[1]},${r.rgb[2]})`) : hashColor(r.id);
    const mat = new THREE.MeshStandardMaterial({ color, roughness: 0.55, metalness: 0.02 });
    const m = new THREE.Mesh(g, mat); m.userData = { layer: L.id, id: r.id }; L.inner.add(m);
    L.regions.set(r.id, { id: r.id, name: r.name, hemi: r.hemi, color, size: r.mm3, sizeLabel: `${r.mm3 < 10 ? r.mm3.toFixed(2) : r.mm3.toFixed(0)} mm³`,
                          voxels: r.voxels, centroid: new THREE.Vector3(...r.centroid), mesh: m, hidden: false });
  }
  if (data.shell) {
    const g = new THREE.BufferGeometry();
    g.setAttribute('position', new THREE.BufferAttribute(b64ToArray(data.shell.v, Float32Array), 3));
    g.setIndex(new THREE.BufferAttribute(b64ToArray(data.shell.f, Uint32Array), 1));
    g.computeVertexNormals();
    L.shell = new THREE.Mesh(g, new THREE.MeshStandardMaterial({ color: 0x8FA19B, transparent: true, opacity: 0.12, depthWrite: false, roughness: 0.9 }));
    L.shell.renderOrder = 10; L.shell.visible = $('chk-shell').checked; L.inner.add(L.shell);
  }
  box.getCenter(L.center); L.radius = box.getSize(new THREE.Vector3()).length() / 2 || 10;
  L.inner.position.copy(L.center).multiplyScalar(-1);
  return L;
}

function buildSurfaceLayer(data, meta) {
  const L = newLayer('surface', meta, `${meta.species} ${data.atlas_name} (${meta.surface.replace('_', ' ')} surface)`);
  L.meta.atlasName = data.atlas_name;
  const labels = new Map(Object.entries(data.labels).map(([k, v]) => [Number(k), v]));
  const box = new THREE.Box3();
  L.hemis = [];
  const counts = new Map();
  for (const h of data.hemis) {
    const pos = b64ToArray(h.v, Float32Array), idx = b64ToArray(h.f, Uint32Array), lab = b64ToArray(h.lab, Int32Array);
    const g = new THREE.BufferGeometry();
    g.setAttribute('position', new THREE.BufferAttribute(pos, 3));
    g.setIndex(new THREE.BufferAttribute(idx, 1));
    const base = new Float32Array(lab.length * 3);
    for (let i = 0; i < lab.length; i++) {
      const info = labels.get(lab[i]);
      const c = info ? (info.rgb ? new THREE.Color(`rgb(${info.rgb[0]},${info.rgb[1]},${info.rgb[2]})`) : hashColor(lab[i])) : new THREE.Color(0xB9C2BE);
      base[3 * i] = c.r; base[3 * i + 1] = c.g; base[3 * i + 2] = c.b;
      if (info) counts.set(lab[i], (counts.get(lab[i]) || 0) + 1);
    }
    g.setAttribute('color', new THREE.BufferAttribute(base.slice(), 3));
    g.computeVertexNormals(); g.computeBoundingBox(); box.union(g.boundingBox);
    const m = new THREE.Mesh(g, new THREE.MeshStandardMaterial({ vertexColors: true, roughness: 0.6, metalness: 0.0, side: THREE.DoubleSide }));
    m.userData = { layer: L.id, surface: true, hemi: h.hemi };
    L.inner.add(m); L.hemis.push({ mesh: m, lab, base });
  }
  for (const [id, info] of labels) {
    const col = info.rgb ? new THREE.Color(`rgb(${info.rgb[0]},${info.rgb[1]},${info.rgb[2]})`) : hashColor(id);
    L.regions.set(id, { id, name: info.name, hemi: info.hemi, color: col, size: counts.get(id) || 0,
                        sizeLabel: plural(counts.get(id) || 0, 'vertex').replace('vertexs', 'vertices'), voxels: counts.get(id) || 0, hidden: false });
  }
  box.getCenter(L.center); L.radius = box.getSize(new THREE.Vector3()).length() / 2 || 10;
  L.inner.position.copy(L.center).multiplyScalar(-1);
  return L;
}

function removeLayer(id) {
  const L = layerById(id); if (!L) return;
  world.remove(L.group);
  L.group.traverse((o) => { if (o.isMesh) { o.geometry.dispose(); o.material.dispose(); } });
  LAYERS.splice(LAYERS.indexOf(L), 1);
  if (ACTIVE === id) ACTIVE = LAYERS.length ? LAYERS[LAYERS.length - 1].id : null;
  if (!LAYERS.length) $('stage-empty').hidden = false;
  arrange(); renderLayers(); renderList(); renderBasket();
}

function setActive(id) {
  ACTIVE = id; const L = layerById(id);
  $('active-layer').innerHTML = L ? `Regions of <b>${esc(L.label)}</b>` : 'No layer yet';
  if (L) { $('rng-opacity').value = L.opacity; $('rng-explode').value = L.explode; }
  document.querySelectorAll('.vol-tool').forEach((e) => { e.hidden = !L || L.kind !== 'volume'; });
  $('btn-to-sba').disabled = $('btn-to-matrix').disabled = !L || L.kind !== 'volume';
  $('btn-to-sba').title = L && L.kind === 'surface' ? 'Seeds are defined on volume atlases: add the volume atlas to send seeds.' : '';
  renderLayers(); renderList(); renderBasket();
}

function renderLayers() {
  $('layers-box').hidden = LAYERS.length === 0;
  $('layer-list').innerHTML = LAYERS.map((L) => `<li data-id="${L.id}" class="${L.id === ACTIVE ? 'on' : ''}" title="Click to make it the active layer">
    <input type="checkbox" data-vis="${L.id}" ${L.visible ? 'checked' : ''} aria-label="show ${esc(L.label)}">
    <span class="lname" title="${esc(L.label)}">${esc(L.label)}</span><button data-rm="${L.id}" aria-label="remove">Remove</button></li>`).join('');
}
$('layer-list').addEventListener('click', (e) => {
  const rm = e.target.closest('[data-rm]'); if (rm) { removeLayer(Number(rm.dataset.rm)); return; }
  const vis = e.target.closest('[data-vis]');
  if (vis) { const L = layerById(Number(vis.dataset.vis)); L.visible = vis.checked; L.group.visible = L.visible; arrange(); return; }
  const li = e.target.closest('li'); if (li) setActive(Number(li.dataset.id));
});
document.querySelectorAll('.segb').forEach((b) => b.addEventListener('click', () => {
  document.querySelectorAll('.segb').forEach((x) => x.classList.toggle('on', x === b)); LAYOUT = b.dataset.layout; arrange(); fitView();
}));
$('chk-samesize').addEventListener('change', () => { arrange(); fitView(); });

function screenAxis() {
  // camera's horizontal screen direction, expressed in the RAS frame of the world group
  camera.updateMatrixWorld();
  const right = new THREE.Vector3().setFromMatrixColumn(camera.matrixWorld, 0);
  return right.applyQuaternion(world.quaternion.clone().invert()).normalize();
}
function arrange() {
  const vis = LAYERS.filter((L) => L.visible);
  const same = $('chk-samesize').checked;
  const ref = vis.length ? Math.max(...vis.map((L) => L.radius)) : 1;
  const axis = screenAxis();
  let x = 0;
  const widths = vis.map((L) => 2 * L.radius * (same ? ref / L.radius : 1));
  const total = widths.reduce((a, b) => a + b, 0) + Math.max(0, vis.length - 1) * ref * 0.25;
  vis.forEach((L, i) => {
    const s = same ? ref / L.radius : 1;
    L.group.scale.setScalar(s);
    if (LAYOUT === 'side') { L.group.position.copy(axis).multiplyScalar(-total / 2 + x + widths[i] / 2); x += widths[i] + ref * 0.25; }
    else L.group.position.set(0, 0, 0);
  });
  refreshAllStyles();
}

// ---- views & orientation
function sceneRadius() {
  const vis = LAYERS.filter((L) => L.visible); if (!vis.length) return 10;
  const box = new THREE.Box3(); vis.forEach((L) => box.expandByObject(L.group));
  return box.getSize(new THREE.Vector3()).length() / 2 || 10;
}
let lastView = 'left';
function fitView(v) {
  v = v || lastView; lastView = v;
  placeCamera(v);
  if (LAYOUT === 'side') { arrange(); placeCamera(v); }   // side by side follows the screen's horizontal axis
}
function placeCamera(v) {
  const vis = LAYERS.filter((L) => L.visible);
  const box = new THREE.Box3(); vis.forEach((L) => box.expandByObject(L.group));
  const c = vis.length ? box.getCenter(new THREE.Vector3()) : new THREE.Vector3();
  const d = sceneRadius() / Math.tan((camera.fov * Math.PI) / 360) * 1.05;
  const dir = { left: [-1, 0, 0.0001], right: [1, 0, 0.0001], top: [0, 1, 0.0001], bottom: [0, -1, 0.0001],
                front: [0, 0, -1], back: [0, 0, 1] }[v];
  camera.position.set(c.x + dir[0] * d, c.y + dir[1] * d, c.z + dir[2] * d);
  camera.near = d / 100; camera.far = d * 20; camera.updateProjectionMatrix();
  controls.target.copy(c); controls.update();
}
document.querySelectorAll('[data-view3d]').forEach((b) => b.addEventListener('click', () => fitView(b.dataset.view3d)));
document.querySelectorAll('[data-rot]').forEach((b) => b.addEventListener('click', () => {
  const axis = { x: [1, 0, 0], y: [0, 1, 0], z: [0, 0, 1] }[b.dataset.rot];
  userQ = new THREE.Quaternion().setFromAxisAngle(new THREE.Vector3(...axis), Math.PI / 2).multiply(userQ);
  applyWorldRotation(); fitView();
}));
$('btn-rot-reset').addEventListener('click', () => { userQ = new THREE.Quaternion(); applyWorldRotation(); fitView(); });

// ---- styles
function refreshStyles(L) {
  const iso = $('chk-isolate').checked && L.selected.size > 0;
  const dim = L.selected.size > 0 && !(QUIZ.active && QUIZ.layer === L.id);
  if (L.kind === 'volume') {
    for (const r of L.regions.values()) {
      const sel = L.selected.has(r.id);
      r.mesh.visible = !r.hidden && (!iso || sel);
      const op = dim && !sel ? Math.min(L.opacity, 0.18) : L.opacity;
      r.mesh.material.opacity = op; r.mesh.material.transparent = op < 1; r.mesh.material.depthWrite = op >= 0.5;
      r.mesh.material.emissive.copy(sel ? EOSIN : new THREE.Color(0x000000));
      r.mesh.material.emissiveIntensity = sel ? 0.55 : 0;
      if (L.hover === r.id && !sel) { r.mesh.material.emissive.set(0x3a3a3a); r.mesh.material.emissiveIntensity = 1; }
    }
    for (const r of L.regions.values()) r.mesh.position.copy(r.centroid.clone().sub(L.center).multiplyScalar(L.explode));
    if (L.shell) L.shell.visible = $('chk-shell').checked;
  } else {
    const grey = new THREE.Color(0xC9D1CD);
    for (const h of L.hemis) {
      const col = h.mesh.geometry.getAttribute('color');
      for (let i = 0; i < h.lab.length; i++) {
        const id = h.lab[i]; const sel = L.selected.has(id); const hov = L.hover === id;
        let r = h.base[3 * i], g = h.base[3 * i + 1], b = h.base[3 * i + 2];
        if (sel) { r = r * 0.35 + EOSIN.r * 0.65; g = g * 0.35 + EOSIN.g * 0.65; b = b * 0.35 + EOSIN.b * 0.65; }
        else if (iso || dim) { r = r * 0.25 + grey.r * 0.75; g = g * 0.25 + grey.g * 0.75; b = b * 0.25 + grey.b * 0.75; }
        if (hov && !sel) { r = Math.min(1, r + 0.18); g = Math.min(1, g + 0.18); b = Math.min(1, b + 0.18); }
        col.array[3 * i] = r; col.array[3 * i + 1] = g; col.array[3 * i + 2] = b;
      }
      col.needsUpdate = true;
      h.mesh.material.opacity = L.opacity; h.mesh.material.transparent = L.opacity < 1;
    }
  }
}
function refreshAllStyles() { LAYERS.forEach(refreshStyles); }
$('rng-opacity').addEventListener('input', () => { const L = layerById(ACTIVE); if (L) { L.opacity = Number($('rng-opacity').value); refreshStyles(L); } });
$('rng-explode').addEventListener('input', () => { const L = layerById(ACTIVE); if (L) { L.explode = Number($('rng-explode').value); refreshStyles(L); } });
$('chk-shell').addEventListener('change', refreshAllStyles);
$('chk-isolate').addEventListener('change', refreshAllStyles);

// ---- picking
const ray = new THREE.Raycaster(); const ptr = new THREE.Vector2();
function pick(ev) {
  const rect = renderer.domElement.getBoundingClientRect();
  ptr.set(((ev.clientX - rect.left) / rect.width) * 2 - 1, -((ev.clientY - rect.top) / rect.height) * 2 + 1);
  ray.setFromCamera(ptr, camera);
  const objs = [];
  for (const L of LAYERS) {
    if (!L.visible) continue;
    if (L.kind === 'volume') for (const r of L.regions.values()) { if (r.mesh.visible) objs.push(r.mesh); }
    else L.hemis.forEach((h) => objs.push(h.mesh));
  }
  const hit = ray.intersectObjects(objs, false)[0];
  if (!hit) return null;
  const L = layerById(hit.object.userData.layer);
  if (hit.object.userData.surface) {
    const h = L.hemis.find((x) => x.mesh === hit.object);
    const id = h.lab[hit.face.a];
    return L.regions.has(id) ? { layer: L, id } : null;
  }
  return { layer: L, id: hit.object.userData.id };
}
let moveRaf = 0;
renderer.domElement.addEventListener('pointermove', (ev) => {
  if (moveRaf) return;
  moveRaf = requestAnimationFrame(() => {
    moveRaf = 0;
    const p = pick(ev); const tip = $('tip');
    for (const L of LAYERS) {
      const nh = p && p.layer === L ? p.id : null;
      if (L.hover !== nh) { L.hover = nh; refreshStyles(L); }
    }
    if (!p || QUIZ.active) { tip.hidden = true; return; }
    const r = p.layer.regions.get(p.id); const rect = stage.getBoundingClientRect();
    tip.innerHTML = `<b>${esc(r.name)}</b>ID ${r.id}; ${esc(r.sizeLabel)}${LAYERS.length > 1 ? `<br><span style="color:#5B6B68">${esc(p.layer.label)}</span>` : ''}`;
    tip.style.left = `${ev.clientX - rect.left + 14}px`; tip.style.top = `${ev.clientY - rect.top + 14}px`; tip.hidden = false;
  });
});
renderer.domElement.addEventListener('pointerleave', () => { $('tip').hidden = true; LAYERS.forEach((L) => { if (L.hover !== null) { L.hover = null; refreshStyles(L); } }); });
let downAt = null;
renderer.domElement.addEventListener('pointerdown', (ev) => { downAt = [ev.clientX, ev.clientY]; });
renderer.domElement.addEventListener('pointerup', (ev) => {
  if (!downAt || Math.hypot(ev.clientX - downAt[0], ev.clientY - downAt[1]) > 5) return;   // it was a rotation
  const p = pick(ev); if (!p) return;
  if (QUIZ.active) { if (p.layer.id === QUIZ.layer) quizAnswer(p.id); return; }
  if (p.layer.id !== ACTIVE) setActive(p.layer.id);
  toggle(p.id);
});
renderer.domElement.addEventListener('dblclick', (ev) => { const p = pick(ev); if (p) focusRegion(p.layer, p.id); });

function regionWorldPoint(L, id) {
  const r = L.regions.get(id);
  if (L.kind === 'volume') return r.mesh.localToWorld(r.centroid.clone());
  const acc = new THREE.Vector3(); let n = 0; const v = new THREE.Vector3();
  for (const h of L.hemis) {
    const pos = h.mesh.geometry.getAttribute('position');
    for (let i = 0; i < h.lab.length; i += 3) if (h.lab[i] === id) { v.fromBufferAttribute(pos, i); acc.add(h.mesh.localToWorld(v.clone())); n++; }
  }
  return n ? acc.divideScalar(n) : null;
}
function focusRegion(L, id) { const p = regionWorldPoint(L, id); if (p) { controls.target.copy(p); controls.update(); } }

// ================================================================== region list & selection
let HEMI = 'all'; let QUERY = '';
function normName(n) {
  // hemisphere-agnostic region name: 'L_Somatosensory', 'ctx-rh-somatosensory', 'Somatosensory_R' -> 'somatosensory'
  return String(n).toLowerCase().replace(/^(ctx[-_])?(lh|rh|l|r|left|right)[-_.\s]+/, '').replace(/[-_.\s]+(lh|rh|l|r|left|right)$/, '')
    .replace(/[^a-z0-9]+/g, '_').replace(/^_|_$/g, '');
}
function propagate(L, id, on) {
  const mode = $('sel-link').value; if (mode === 'off') return 0;
  const src = L.regions.get(id); const key = normName(src.name); let n = 0;
  for (const M of LAYERS) {
    if (M === L) continue;
    for (const r of M.regions.values()) {
      if (mode === 'id' ? r.id === id : normName(r.name) === key) { if (on) M.selected.add(r.id); else M.selected.delete(r.id); n++; }
    }
    refreshStyles(M);
  }
  return n;
}
function toggle(id) {
  const L = layerById(ACTIVE); if (!L) return;
  if (L.selected.has(id)) L.selected.delete(id); else L.selected.add(id);
  const linked = propagate(L, id, L.selected.has(id));
  if ($('sel-link').value !== 'off') $('basket-msg').textContent = linked ? `Linked: ${plural(linked, 'region')} in the other layers.` : 'No matching region in the other layers.';
  refreshStyles(L); renderBasket();
  const li = document.querySelector(`#region-list li[data-id="${id}"]`);
  if (li) { const on = L.selected.has(id); li.classList.toggle('sel', on); li.querySelector('input').checked = on; }
}
document.querySelectorAll('.chip').forEach((c) => c.addEventListener('click', () => {
  document.querySelectorAll('.chip').forEach((x) => x.classList.toggle('on', x === c)); HEMI = c.dataset.hemi; renderList();
}));
$('region-search').addEventListener('input', (e) => { QUERY = e.target.value.trim().toLowerCase(); renderList(); });
function renderList() {
  const L = layerById(ACTIVE);
  if (!L) { $('region-list').innerHTML = ''; $('region-count').textContent = ''; return; }
  const list = [...L.regions.values()].filter((r) => (HEMI === 'all' || r.hemi === HEMI) &&
      (!QUERY || r.name.toLowerCase().includes(QUERY) || String(r.id) === QUERY))
    .sort((a, b) => a.name.localeCompare(b.name, 'en', { numeric: true }));
  $('region-count').textContent = `${list.length} of ${plural(L.regions.size, 'region')} shown`;
  $('region-list').innerHTML = list.map((r) => `<li data-id="${r.id}" class="${L.selected.has(r.id) ? 'sel' : ''}" title="Click: select. Double-click: centre the view on it.">
    <input type="checkbox" ${L.selected.has(r.id) ? 'checked' : ''} tabindex="-1" aria-label="select ${esc(r.name)}">
    <span class="sw" style="background:#${r.color.getHexString()}"></span><span class="nm">${esc(r.name)}</span>
    <span class="vol">${r.id}; ${esc(r.sizeLabel)}</span></li>`).join('');
}
$('region-list').addEventListener('click', (e) => { const li = e.target.closest('li'); if (li) toggle(Number(li.dataset.id)); });
$('region-list').addEventListener('dblclick', (e) => { const li = e.target.closest('li'); if (li) focusRegion(layerById(ACTIVE), Number(li.dataset.id)); });
function renderBasket() {
  const L = layerById(ACTIVE); const ids = L ? [...L.selected] : [];
  $('basket-n').textContent = `${plural(ids.length, 'region')} selected`;
  $('basket-list').textContent = ids.map((i) => `${L.regions.get(i)?.name ?? i} (${i})`).join(', ');
}
$('btn-clear').addEventListener('click', () => {
  const L = layerById(ACTIVE); if (!L) return;
  const targets = $('sel-link').value === 'off' ? [L] : LAYERS;
  targets.forEach((M) => { M.selected.clear(); refreshStyles(M); });
  renderList(); renderBasket();
});
$('btn-copy').addEventListener('click', async () => {
  const L = layerById(ACTIVE); if (!L) return;
  const txt = `${L.label}\n` + [...L.selected].map((i) => `${i}\t${L.regions.get(i)?.name}`).join('\n');
  try { await navigator.clipboard.writeText(txt); $('basket-msg').textContent = 'List copied to the clipboard.'; }
  catch { $('basket-msg').textContent = txt; }
});

// ---- send the selection to the functional launcher
function fieldByName(S, name) { return S.model ? S.model.fields.filter((f) => f.name === name).pop() : null; }
function pyList(text) {
  return JSON.parse(text.replace(/'/g, '"').replace(/\bTrue\b/g, 'true').replace(/\bFalse\b/g, 'false').replace(/\bNone\b/g, 'null').replace(/\(/g, '[').replace(/\)/g, ']'));
}
function toPy(v) { return JSON.stringify(v).replace(/"/g, "'").replace(/,/g, ', '); }
$('btn-to-sba').addEventListener('click', () => {
  const m = $('basket-msg'); const L = layerById(ACTIVE); const S = PIPES.func;
  if (!L || !L.selected.size) { m.textContent = 'Select at least one region (click it in the 3D view or in the list).'; return; }
  const fa = fieldByName(S, 'selected_atlases_SBA'), fb = fieldByName(S, 'panda_files_SBA');
  if (!fa || !fb) { m.textContent = 'Open a functional launcher first (Launchers, Functional).'; return; }
  try {
    const a = pyList(S.values[fa.key]); const b = pyList(S.values[fb.key]);
    a.push([L.meta.atlasName, L.meta.level]); b.push([...L.selected]);
    setField(S, fa, toPy(a)); setField(S, fb, toPy(b));
    m.textContent = `Seeds added to ${S.model.path.split('/').pop()}; remember to save it.`;
  } catch (e) { m.textContent = `Could not read the current launcher value: ${e.message}`; }
});
$('btn-to-matrix').addEventListener('click', () => {
  const m = $('basket-msg'); const L = layerById(ACTIVE); const S = PIPES.func; const f = fieldByName(S, 'selected_atlases_matrix');
  if (!f || !L) { m.textContent = 'Open a functional launcher first, then add a volume atlas layer.'; return; }
  let cur = S.values[f.key];
  cur = (String(cur).trim().startsWith('[')) ? pyList(cur) : [[], []];
  cur[0].push(L.meta.atlasName); cur[1].push(L.meta.level);
  setField(S, f, toPy(cur));
  m.textContent = `${L.meta.atlasName} level ${L.meta.level} added to the matrices; remember to save the launcher.`;
});

// ================================================================== challenge (quiz)
const QUIZ = { active: false, layer: null, round: 0, score: 0, streak: 0, target: null, pool: [], busy: false, kept: [] };
$('btn-quiz').addEventListener('click', startQuiz);
$('quiz-stop').addEventListener('click', stopQuiz);
function quizKey(L) { return `ednixQuizBest:${L.label}`; }
function startQuiz() {
  const L = layerById(ACTIVE);
  if (!L) { banner('Add an atlas layer before starting the challenge.'); return; }
  const regs = [...L.regions.values()].filter((r) => !r.hidden && r.size > 0);
  const cut = regs.map((r) => r.size).sort((a, b) => a - b)[Math.floor(regs.length * 0.15)] || 0;
  QUIZ.pool = regs.filter((r) => r.size >= cut).map((r) => r.id);
  Object.assign(QUIZ, { active: true, layer: L.id, round: 0, score: 0, streak: 0, kept: [...L.selected] });
  L.selected.clear(); refreshStyles(L); renderList(); renderBasket();
  $('quiz').hidden = false; $('btn-quiz').hidden = true; document.querySelector('.views').hidden = true; nextQuestion();
}
function stopQuiz() {
  const L = layerById(QUIZ.layer);
  const was = QUIZ.active; QUIZ.active = false;
  $('quiz').hidden = true; $('btn-quiz').hidden = false; $('quiz').className = 'quiz'; document.querySelector('.views').hidden = false;
  if (was && L) { QUIZ.kept.forEach((i) => L.selected.add(i)); QUIZ.kept = []; refreshStyles(L); renderList(); renderBasket(); }
}
function nextQuestion() {
  QUIZ.round += 1;
  if (QUIZ.round > 10) return endQuiz();
  const L = layerById(QUIZ.layer);
  QUIZ.target = QUIZ.pool[Math.floor(Math.random() * QUIZ.pool.length)];
  $('quiz').className = 'quiz';
  $('quiz-q').textContent = `Find: ${L.regions.get(QUIZ.target).name}`;
  $('quiz-round').textContent = `question ${QUIZ.round} of 10`;
  $('quiz-score').textContent = `${QUIZ.score} points`;
  $('quiz-streak').textContent = QUIZ.streak > 1 ? `streak of ${QUIZ.streak}` : '';
}
function flash(L, id, color, ms) {
  L.selected.add(id); const keep = EOSIN.getHex(); EOSIN.set(color); refreshStyles(L);
  setTimeout(() => { L.selected.delete(id); EOSIN.set(keep); refreshStyles(L); }, ms);
}
function quizAnswer(id) {
  if (QUIZ.busy) return; QUIZ.busy = true;
  const L = layerById(QUIZ.layer);
  if (id === QUIZ.target) {
    QUIZ.streak += 1; QUIZ.score += 10 + 5 * (QUIZ.streak - 1);
    $('quiz').className = 'quiz good'; $('quiz-q').textContent = 'Correct.'; flash(L, id, 0x2F7A55, 700);
  } else {
    QUIZ.streak = 0; $('quiz').className = 'quiz bad';
    $('quiz-q').textContent = `That was ${L.regions.get(id).name}. The right region is highlighted.`;
    flash(L, QUIZ.target, 0xC8466F, 1600);
  }
  $('quiz-score').textContent = `${QUIZ.score} points`;
  setTimeout(() => { QUIZ.busy = false; if (QUIZ.active) nextQuestion(); }, id === QUIZ.target ? 800 : 1700);
}
function endQuiz() {
  const L = layerById(QUIZ.layer);
  const best = Math.max(QUIZ.score, Number(localStorage.getItem(quizKey(L)) || 0));
  try { localStorage.setItem(quizKey(L), String(best)); } catch { /* storage unavailable */ }
  $('quiz').className = 'quiz good';
  $('quiz-q').textContent = `Finished: ${QUIZ.score} points`;
  $('quiz-round').textContent = `best score on this layer: ${best}`;
  $('quiz-score').textContent = ''; $('quiz-streak').textContent = '';
  QUIZ.active = false; $('btn-quiz').hidden = false; $('btn-quiz').textContent = 'Play the challenge again';
  QUIZ.kept.forEach((i) => L.selected.add(i)); QUIZ.kept = []; refreshStyles(L); renderList(); renderBasket();
}

// ================================================================== launchers (one state per pipeline)
const PIPE_NAMES = { anat: 'Anatomical pipeline', func: 'Functional pipeline', pet: 'PET pipeline' };
const PIPES = { anat: newPipeState(), func: newPipeState(), pet: newPipeState() };
function newPipeState() { return { model: null, values: {}, original: {}, saved: null, out: '', path: '' }; }
function switchPipe(p) {
  if (PIPES[currentPipe]) { PIPES[currentPipe].out = $('out-path').value; PIPES[currentPipe].path = $('launcher-path').value; }
  currentPipe = p;
  $('pipe-title').textContent = PIPE_NAMES[p];
  const list = (INFO?.launchers || []).filter((l) => l.pipeline === p);
  setOptions($('sel-launcher'), list.length ? list.map((l) => [l.path, l.name]) : [['', 'no example launcher']]);
  $('out-path').value = PIPES[p].out; $('launcher-path').value = PIPES[p].path;
  renderForm(); msg([]);
}
$('btn-open').addEventListener('click', openLauncher);
$('sel-launcher').addEventListener('change', () => { $('launcher-path').value = ''; });
async function openLauncher() {
  const S = PIPES[currentPipe];
  const path = $('launcher-path').value.trim() || $('sel-launcher').value;
  if (!path) return;
  let model;
  try { model = await api(`/api/launcher?path=${encodeURIComponent(path)}`); }
  catch (e) { msg([['bad', `Could not open this launcher: ${e.message}`]]); return; }
  const expect = { anat: ['preprocess_anat'], func: ['preprocess_fmri'], pet: ['preprocess_pet', 'preprocess_data'] }[currentPipe];
  Object.assign(S, { model, values: {}, original: {}, saved: null });
  for (const f of model.fields) { S.values[f.key] = f.value; S.original[f.key] = f.value; }
  $('out-path').value = path.replace(/\.py$/, '_custom.py');
  renderForm();
  msg(model.pipeline && !expect.includes(model.pipeline) ? [['wn', `This launcher calls ${model.pipeline}(), not ${expect[0]}().`]]
    : !model.pipeline ? [['wn', 'This script does not call a pipeline function: only its variables are shown.']] : []);
}
function setField(S, f, value) {
  S.values[f.key] = value;
  if (S !== PIPES[currentPipe]) return;
  const el = document.querySelector(`[data-key="${CSS.escape(f.key)}"]`);
  if (el) { if (el.type === 'checkbox') el.checked = value; else el.value = value; markChanged(f.key); }
}
function markChanged(key) {
  const S = PIPES[currentPipe];
  const row = document.querySelector(`.field[data-row="${CSS.escape(key)}"]`);
  if (row) row.classList.toggle('changed', S.values[key] !== S.original[key]);
}
function stepsWidget(f, steps) {
  let cur = [];
  try { cur = pyList(PIPES[currentPipe].values[f.key]); } catch { cur = []; }
  const known = new Set(steps.map((s) => String(s.id)));
  const extra = cur.filter((c) => !known.has(String(c)));
  return `<div class="steps" data-steps="${esc(f.key)}">` + steps.map((s) => `<label><input type="checkbox" data-step="${esc(JSON.stringify(s.id))}"
      ${cur.some((c) => String(c) === String(s.id)) ? 'checked' : ''}><code>${esc(s.id)}</code><span>${esc(s.label)}</span></label>`).join('') +
    (extra.length ? `<div class="note">Other values kept as they are: ${esc(JSON.stringify(extra))}</div>` : '') + '</div>';
}
function inputFor(f, steps) {
  if (f.name === 'Skip_step' && f.kind === 'list' && steps.length) return stepsWidget(f, steps);
  const k = `data-key="${esc(f.key)}"`;
  const v = PIPES[currentPipe].values[f.key];
  if (f.kind === 'bool') return `<label class="check"><input type="checkbox" ${k} ${v ? 'checked' : ''}> ${v ? 'yes' : 'no'}</label>`;
  if (f.kind === 'list') return `<textarea ${k} rows="${String(v).length > 90 ? 3 : 1}" spellcheck="false">${esc(v)}</textarea>`;
  const dl = f.options && f.options.length ? `list="dl-${esc(f.key)}"` : '';
  const opts = dl ? `<datalist id="dl-${esc(f.key)}">${f.options.map((o) => `<option value="${esc(o)}">`).join('')}</datalist>` : '';
  return `<input type="text" ${k} ${dl} value="${esc(v)}" spellcheck="false">${opts}`;
}
function renderForm() {
  const S = PIPES[currentPipe];
  if (!S.model) {
    $('form').innerHTML = `<div class="empty"><p class="big">Open a launcher to edit its parameters without touching the code.</p>
      <p>Pick an example above or give the path of your own launcher. Every parameter comes with its description.</p></div>`;
    return;
  }
  const secs = new Map();
  const stepField = S.model.fields.filter((f) => f.name === 'Skip_step' && f.kind === 'list').pop();
  if (stepField && S.model.steps.length) secs.set('Pipeline steps to skip', [stepField]);
  for (const f of S.model.fields) {
    if (f === stepField && secs.has('Pipeline steps to skip')) continue;
    if (!secs.has(f.section)) secs.set(f.section, []); secs.get(f.section).push(f);
  }
  const warn = S.model.warnings.map((w) => `<div class="warn">${esc(w)}</div>`).join('');
  let i = 0;
  $('form').innerHTML = `<h2 style="margin:4px 0 10px;font-size:18px">${esc(S.model.path.split('/').pop())}</h2>${warn}` +
    [...secs.entries()].map(([s, fs]) => `<details ${i++ < 2 ? 'open' : ''}><summary>${esc(s)}</summary><div class="fields">` +
      fs.map((f) => `<div class="field ${S.values[f.key] !== S.original[f.key] ? 'changed' : ''}" data-row="${esc(f.key)}">
        <div class="name" title="${esc(f.where)}">${esc(f.name)}</div><div>${inputFor(f, S.model.steps)}</div>` +
        (f.doc ? `<div class="doc">${esc(f.doc)}</div>` : '') +
        (f.note && f.note !== f.doc ? `<div class="note">Launcher comment: ${esc(f.note)}</div>` : '') + '</div>').join('') + '</div></details>').join('');
  $('form').querySelectorAll('[data-key]').forEach((el) => {
    el.addEventListener(el.type === 'checkbox' ? 'change' : 'input', () => {
      const key = el.dataset.key;
      S.values[key] = el.type === 'checkbox' ? el.checked : el.value;
      if (el.type === 'checkbox') el.parentElement.lastChild.textContent = el.checked ? ' yes' : ' no';
      markChanged(key);
    });
  });
  $('form').querySelectorAll('[data-steps]').forEach((box) => {
    box.addEventListener('change', () => {
      const key = box.dataset.steps;
      let cur = [];
      try { cur = pyList(S.values[key]); } catch { cur = []; }
      const known = new Set([...box.querySelectorAll('[data-step]')].map((c) => String(JSON.parse(c.dataset.step))));
      const extra = cur.filter((c) => !known.has(String(c)));
      const ticked = [...box.querySelectorAll('[data-step]:checked')].map((c) => JSON.parse(c.dataset.step));
      S.values[key] = toPy([...ticked, ...extra]);
      markChanged(key);
    });
  });
}
function msg(items) { $('launcher-msg').innerHTML = items.map(([c, t]) => `<p class="${c}">${esc(t)}</p>`).join(''); }
async function checkOrSave(save) {
  const S = PIPES[currentPipe];
  if (!S.model) { msg([['bad', 'Open a launcher first.']]); return null; }
  const body = { path: S.model.path, values: S.values, out: $('out-path').value.trim() };
  try {
    const r = await api(save ? '/api/launcher/save' : '/api/launcher/check', body);
    if (!r.ok) { msg(r.errors.map((e) => ['bad', e])); return null; }
    const lines = [['ok', save ? `Saved: ${r.written} (${plural(r.changed, 'value')} changed).` : `No syntax error (${plural(r.changed, 'value')} changed).`]];
    if (r.backup) lines.push(['wn', `Previous version kept as ${r.backup}`]);
    r.issues.forEach((i) => lines.push(['wn', i]));
    if (!r.issues.length) lines.push(['ok', 'Every argument is known to the pipeline.']);
    msg(lines);
    if (save) S.saved = r.written;
    return r;
  } catch (e) { msg([['bad', e.message]]); return null; }
}
$('btn-mainpath').addEventListener('click', () => {
  const S = PIPES[currentPipe]; const f = fieldByName(S, 'MAIN_PATH');
  if (!f || !INFO?.repo) { msg([['bad', 'Open a launcher that defines MAIN_PATH first.']]); return; }
  setField(S, f, INFO.repo.endsWith('/') ? INFO.repo : `${INFO.repo}/`);
  msg([['ok', `MAIN_PATH = ${INFO.repo} (remember to save).`]]);
});
$('btn-check').addEventListener('click', () => checkOrSave(false));
$('btn-save').addEventListener('click', () => checkOrSave(true));
$('btn-run').addEventListener('click', async () => {
  const S = PIPES[currentPipe];
  if (!S.model) { msg([['bad', 'Open a launcher first.']]); return; }
  const dirty = S.model.fields.some((f) => S.values[f.key] !== S.original[f.key]);
  if (dirty || !S.saved) { const r = await checkOrSave(true); if (!r) return; }
  try {
    const r = await api('/api/run', { path: S.saved || S.model.path, strict: $('chk-strict').checked });
    msg([['ok', `Run started (${r.id}). Follow it in Runs.`]]);
    RUNS.selected = r.id; showView('runs');
  } catch (e) { msg([['bad', `The run did not start: ${e.message}`]]); }
});

// ================================================================== runs
const RUNS = { list: [], selected: null };
function stClass(s) { return s === 'running' ? 'st-running' : s === 'finished' ? 'st-ok' : s === 'failed' ? 'st-fail' : 'st-other'; }
async function refreshRuns() {
  try { RUNS.list = await api('/api/runs'); } catch { return; }
  refreshBadge(true);
  $('runs-empty').hidden = RUNS.list.length > 0;
  $('runs-tbody').innerHTML = RUNS.list.map((r) => `<tr data-id="${r.id}" class="${r.id === RUNS.selected ? 'on' : ''}">
    <td>${esc(r.name)}</td><td>${esc(r.started.replace('T', ' '))}</td>
    <td class="st ${stClass(r.status)}">${esc(r.status)}${r.code !== null && r.code !== undefined && r.code !== 0 ? ` (code ${r.code})` : ''}</td>
    <td class="prog">${r.progress ? `step ${esc(r.progress.step)}: ${esc(r.progress.function)}` : ''}</td></tr>`).join('');
  if (!RUNS.selected && RUNS.list.length) RUNS.selected = RUNS.list[0].id;
  if (RUNS.selected) showLog(RUNS.selected);
}
function refreshBadge(fromList = false) {
  const setBadge = (l) => { const n = l.filter((r) => r.status === 'running').length; $('run-badge').hidden = n === 0; $('run-badge').textContent = n; };
  if (fromList) setBadge(RUNS.list);
  else if (currentView !== 'runs') api('/api/runs').then((l) => { RUNS.list = l; setBadge(l); }).catch(() => {});
}
$('runs-tbody').addEventListener('click', (e) => { const tr = e.target.closest('tr'); if (!tr) return; RUNS.selected = tr.dataset.id; refreshRuns(); });
let LOG_TEXT = '';
async function showLog(id) {
  const run = RUNS.list.find((r) => r.id === id);
  $('log-title').textContent = run ? `Log of ${run.name}` : 'Log';
  $('btn-stop').disabled = !run || run.status !== 'running';
  try {
    const { text } = await api(`/api/run/log?id=${encodeURIComponent(id)}`);
    if (text === LOG_TEXT && $('log').innerHTML) return;
    LOG_TEXT = text || '';
    paintRunLog();
  } catch { /* keep previous text */ }
}
function paintRunLog() {
  const pre = $('log'); const atEnd = pre.scrollTop + pre.clientHeight >= pre.scrollHeight - 30;
  const lines = LOG_TEXT ? LOG_TEXT.split('\n').map((t, i) => [i + 1, t]) : [];
  const r = renderLines(lines, $('chk-log-issues').checked);
  pre.innerHTML = r.html || '(empty log)';
  $('log-counts').innerHTML = LOG_TEXT ? countsHtml(r.errors, r.warnings) : '';
  if (atEnd) pre.scrollTop = pre.scrollHeight;
}
$('chk-log-issues').addEventListener('change', paintRunLog);
$('btn-next-error').addEventListener('click', () => scrollToNextError($('log')));
$('btn-stop').addEventListener('click', async () => {
  if (!RUNS.selected) return;
  await api('/api/run/stop', { id: RUNS.selected }).catch(() => {});
  setTimeout(refreshRuns, 800);
});
$('btn-quit').addEventListener('click', async () => {
  await api('/api/shutdown', {}).catch(() => {});
  document.body.innerHTML = '<p style="margin:40px;font-size:18px">Interface closed. Runs already started keep going in the background; double-click EDNiX again to follow them.</p>';
});

// ================================================================== quality control
// ------------------------------------------------------------------ log rendering (runs and processing logs)
const ERR_RE = /\bERROR\b|Traceback \(most recent call last\)|exited with code [1-9]|^\s*(?:\*\* )?error:|^\s*\w*(?:Error|Exception)\b(?::|$)|\bFATAL\b|Segmentation fault/;
const WARN_RE = /\bWARNING\b|^\s*(?:\*\* )?warning:/;
const ANSI_CLASS = { 91: 'l-err', 93: 'l-warn', 92: 'l-ok', 96: 'l-info', 94: 'l-info', 95: 'l-head', 1: 'l-bold' };
function lineKind(plain, ansiCodes) {
  if (ansiCodes.includes(91) || ERR_RE.test(plain)) return 'error';
  if (ansiCodes.includes(93) || WARN_RE.test(plain)) return 'warning';
  return '';
}
function ansiLine(line) {
  // pipeline messages carry ANSI colours (run_cmd.printcolor): keep them, as HTML
  let out = '', open = false; const codes = [];
  line.split(/\x1b\[([0-9;]*)m/).forEach((part, i) => {
    if (i % 2 === 0) { out += esc(part); return; }
    if (open) { out += '</span>'; open = false; }
    const c = part.split(';').map(Number); codes.push(...c);
    const cls = c.map((x) => ANSI_CLASS[x]).filter(Boolean);
    if (cls.length) { out += `<span class="${cls.join(' ')}">`; open = true; }
  });
  if (open) out += '</span>';
  return { html: out, plain: line.replace(/\x1b\[[0-9;]*m/g, ''), codes };
}
// lines: [[number, text, kind?], ...]; onlyIssues keeps errors / warnings with 2 lines of context
function renderLines(lines, onlyIssues) {
  const parsed = lines.map(([n, text, kind]) => { const a = ansiLine(text); return { n, html: a.html, kind: kind !== undefined ? kind : lineKind(a.plain, a.codes) }; });
  let keep = parsed.map(() => true);
  if (onlyIssues) {
    keep = parsed.map(() => false);
    parsed.forEach((l, i) => { if (l.kind) for (let j = Math.max(0, i - 2); j <= Math.min(parsed.length - 1, i + 2); j++) keep[j] = true; });
  }
  let html = '', prev = -1;
  parsed.forEach((l, i) => {
    if (!keep[i]) return;
    if (onlyIssues && prev >= 0 && i !== prev + 1) html += '<span class="gap">   ...</span>';
    html += `<span class="ll ${l.kind ? 'k-' + l.kind : ''}"><span class="ln">${l.n}</span>${l.html || ' '}</span>`;
    prev = i;
  });
  return { html, errors: parsed.filter((l) => l.kind === 'error').length, warnings: parsed.filter((l) => l.kind === 'warning').length };
}
function countsHtml(e, w) {
  return `<span class="${e ? 'e' : ''}">${plural(e, 'error')}</span>, <span class="${w ? 'w' : ''}">${plural(w, 'warning')}</span>`;
}
function scrollToNextError(container) {
  const errs = [...container.querySelectorAll('.k-error')]; if (!errs.length) return;
  const top = container.scrollTop + 5;
  const next = errs.find((el) => el.offsetTop > top + 1) || errs[0];
  container.scrollTop = next.offsetTop - 40;
}

// ------------------------------------------------------------------ quality control
const QC = { bids: '', scan: null, mode: 'images', step: null, table: null, tableData: null, sort: null, logs: null, log: null, loaded: false };
async function loadQcCandidates() {
  if (QC.loaded) return; QC.loaded = true;
  let c = [];
  try { c = await api('/api/qc/candidates'); } catch { c = []; }
  setOptions($('sel-qc-bids'), c.length ? c.map((x) => [x.path, `${x.path}${x.has_qc ? '' : ' (no QC folder yet)'}`]) : [['', 'no dataset found on this machine']]);
}
const REPORT_GROUPS = ['Functional reports', 'Anatomical reports'];
$('btn-qc-scan').addEventListener('click', async () => {
  QC.bids = $('qc-path').value.trim() || $('sel-qc-bids').value;
  if (!QC.bids) return;
  try { QC.scan = await api(`/api/qc/scan?bids=${encodeURIComponent(QC.bids)}`); }
  catch (e) { banner(`Could not read the dataset: ${e.message}`, true); return; }
  QC.logs = null;
  if (!QC.scan.exists) { ['qc-steps', 'qc-tables', 'qc-logs'].forEach((id) => { $(id).innerHTML = ''; }); $('qc-grid').innerHTML = `<div class="empty" style="position:static"><p>Nothing to show in ${esc(QC.bids)}: no QC folder, no report and no processing log.</p></div>`; return; }
  QC.scan.steps.sort((a, b) => (REPORT_GROUPS.includes(b.name) - REPORT_GROUPS.includes(a.name)) || REPORT_GROUPS.indexOf(a.name) - REPORT_GROUPS.indexOf(b.name));
  QC.mode = 'images'; QC.step = QC.scan.steps.length ? QC.scan.steps[0].name : null; QC.table = null;
  renderQcSide(); renderQcMain();
  if (QC.scan.n_logs) { try { QC.logs = await api(`/api/qc/logs?bids=${encodeURIComponent(QC.bids)}`); } catch { QC.logs = []; } renderQcSide(); }
});
function renderQcSide() {
  $('qc-steps').innerHTML = QC.scan.steps.map((st) => `<li data-step="${esc(st.name)}" class="${REPORT_GROUPS.includes(st.name) ? 'group ' : ''}${QC.mode === 'images' && st.name === QC.step ? 'on' : ''}"><span>${esc(st.name.replace(/_/g, ' '))}</span><span class="n">${st.images.length}</span></li>`).join('') || '<li>none</li>';
  $('qc-tables').innerHTML = QC.scan.tables.map((t) => `<li data-table="${esc(t.rel)}" class="${t.rel.startsWith('@') ? 'group ' : ''}${QC.mode === 'table' && t.rel === QC.table ? 'on' : ''}"><span>${esc(t.name)}</span></li>`).join('') || '<li>none</li>';
  if (!QC.scan.n_logs) { $('qc-logs').innerHTML = '<li>none</li>'; return; }
  const withErr = QC.logs ? QC.logs.filter((l) => l.errors).length : null;
  $('qc-logs').innerHTML = `<li data-logs="1" class="${QC.mode === 'logs' || QC.mode === 'log' ? 'on' : ''}"><span>Errors and warnings</span>` +
    (withErr === null ? `<span class="n">${QC.scan.n_logs} logs</span>` : `<span class="n ${withErr ? 'err' : ''}">${withErr ? plural(withErr, 'log') + ' with errors' : 'no error'}</span>`) + '</li>';
}
$('qc-steps').addEventListener('click', (e) => { const li = e.target.closest('[data-step]'); if (!li) return; QC.mode = 'images'; QC.step = li.dataset.step; renderQcSide(); renderQcMain(); });
$('qc-tables').addEventListener('click', async (e) => { const li = e.target.closest('[data-table]'); if (!li) return; QC.mode = 'table'; QC.table = li.dataset.table; QC.tableData = null; QC.sort = null; renderQcSide(); await renderQcMain(); });
$('qc-logs').addEventListener('click', async (e) => { const li = e.target.closest('[data-logs]'); if (!li) return; QC.mode = 'logs'; renderQcSide(); await renderQcMain(); });
$('qc-filter').addEventListener('input', () => { renderQcMain(); });
async function renderQcMain() {
  if (QC.mode === 'images') return renderQcGrid();
  if (QC.mode === 'table') return renderQcTable();
  if (QC.mode === 'logs') return renderQcLogs();
  if (QC.mode === 'log') return renderQcLog();
}
const qcUrl = (rel) => `/api/qc/file?bids=${encodeURIComponent(QC.bids)}&rel=${encodeURIComponent(rel)}&token=${encodeURIComponent(TOKEN)}`;
let QC_VISIBLE = [];
function renderQcGrid() {
  const st = QC.scan.steps.find((x) => x.name === QC.step);
  const f = $('qc-filter').value.trim().toLowerCase();
  QC_VISIBLE = st ? st.images.filter((r) => !f || r.toLowerCase().includes(f)) : [];
  $('qc-count').textContent = st ? `${QC_VISIBLE.length} of ${plural(st.images.length, 'image')}` : '';
  const caption = (r) => { const sub = (r.match(/sub-[^/]+/) || [''])[0]; const ses = (r.match(/ses-[^/]+/) || [''])[0]; const name = r.split('/').pop(); return (sub && !name.includes(sub) ? `${sub} ${ses} · ` : '') + name; };
  $('qc-grid').innerHTML = QC_VISIBLE.slice(0, 400).map((r, i) => `<div class="qc-card" data-i="${i}"><img loading="lazy" src="${qcUrl(r)}" alt="${esc(r)}"><div class="cap">${esc(caption(r))}</div></div>`).join('')
    + (QC_VISIBLE.length > 400 ? '<p class="count">Only the first 400 images are shown: use the filter.</p>' : '');
}
async function renderQcTable() {
  if (!QC.tableData) {
    try { QC.tableData = await api(`/api/qc/table?bids=${encodeURIComponent(QC.bids)}&rel=${encodeURIComponent(QC.table)}`); }
    catch (e) { $('qc-grid').innerHTML = `<p class="count">Could not read ${esc(QC.table)}: ${esc(e.message)}</p>`; return; }
  }
  const t = QC.tableData;
  const f = $('qc-filter').value.trim().toLowerCase();
  let rows = t.rows.map((r, i) => [r, t.flags[i]]).filter(([r]) => !f || r.join(' ').toLowerCase().includes(f));
  if (QC.sort) {
    const { col, dir } = QC.sort;
    rows.sort((a, b) => { const x = a[0][col], y = b[0][col]; if (x === y) return 0; if (x === null) return 1; if (y === null) return -1; return (typeof x === 'number' && typeof y === 'number' ? x - y : String(x).localeCompare(String(y))) * dir; });
  }
  const subCol = t.columns.indexOf('subject'); const runCol = t.columns.indexOf('run');
  const linkable = t.link && subCol >= 0;
  $('qc-count').textContent = `${plural(rows.length, 'row')}; ${plural(t.n_flagged_rows, 'row')} with a value beyond ${t.z_thresh} robust SD` + (linkable ? '; click a row to see its images' : '');
  $('qc-grid').innerHTML = `<div class="qc-table"><table><thead><tr>${t.columns.map((c, j) => `<th class="sortable ${QC.sort && QC.sort.col === j ? 'sorted' : ''}" data-col="${j}">${esc(c)}${QC.sort && QC.sort.col === j ? (QC.sort.dir > 0 ? ' ▲' : ' ▼') : ''}</th>`).join('')}</tr></thead><tbody>` +
    rows.map(([r, fl]) => `<tr class="${fl.some(Boolean) ? 'flagrow ' : ''}${linkable ? 'link' : ''}" ${linkable ? `data-sub="${esc(r[subCol])}" data-run="${runCol >= 0 ? esc(r[runCol]) : ''}"` : ''}>${r.map((v, j) => `<td class="${fl[j] ? 'flag' : ''}">${v === null ? '' : esc(v)}</td>`).join('')}</tr>`).join('') + '</tbody></table></div>';
}
$('qc-grid').addEventListener('click', (e) => {
  const th = e.target.closest('th.sortable');
  if (th && QC.mode === 'table') { const col = Number(th.dataset.col); QC.sort = QC.sort && QC.sort.col === col ? { col, dir: -QC.sort.dir } : { col, dir: 1 }; renderQcTable(); return; }
  const tr = e.target.closest('tr.link');
  if (tr && QC.mode === 'table') {
    const group = QC.tableData.link; if (!QC.scan.steps.find((s) => s.name === group)) { banner(`No ${group.toLowerCase()} image in this dataset.`, true); return; }
    $('qc-filter').value = tr.dataset.run || `sub-${tr.dataset.sub}`;
    QC.mode = 'images'; QC.step = group; renderQcSide(); renderQcGrid(); return;
  }
  const lr = e.target.closest('tr[data-log]');
  if (lr) { QC.log = { rel: lr.dataset.log, onlyIssues: true }; QC.mode = 'log'; renderQcSide(); renderQcLog(); return; }
  if (e.target.id === 'qc-log-back') { QC.mode = 'logs'; renderQcLogs(); return; }
  if (e.target.id === 'qc-log-next') { scrollToNextError($('qc-log-text')); return; }
  const c = e.target.closest('.qc-card'); if (c) showLightbox(Number(c.dataset.i));
});
$('qc-grid').addEventListener('change', (e) => { if (e.target.id === 'qc-log-issues') { QC.log.onlyIssues = e.target.checked; renderQcLog(); } });
async function renderQcLogs() {
  if (!QC.logs) { try { QC.logs = await api(`/api/qc/logs?bids=${encodeURIComponent(QC.bids)}`); } catch (e) { $('qc-grid').innerHTML = `<p class="count">${esc(e.message)}</p>`; return; } renderQcSide(); }
  const f = $('qc-filter').value.trim().toLowerCase();
  const rows = QC.logs.filter((l) => !f || `${l.subject} ${l.session} ${l.modality} ${l.rel} ${l.last_error}`.toLowerCase().includes(f));
  const nErr = rows.filter((l) => l.errors).length;
  $('qc-count').textContent = `${plural(rows.length, 'processing log')}, ${nErr} with errors; click a row to read it`;
  $('qc-grid').innerHTML = `<div class="qc-table"><table><thead><tr><th>subject</th><th>session</th><th>modality</th><th>errors</th><th>warnings</th><th style="text-align:left">last error</th><th style="text-align:left">log</th></tr></thead><tbody>` +
    rows.map((l) => `<tr class="link ${l.errors ? 'flagrow' : ''}" data-log="${esc(l.rel)}"><td>${esc(l.subject)}</td><td>${esc(l.session)}</td><td>${esc(l.modality)}</td><td class="${l.errors ? 'flag' : ''}">${l.errors}</td><td>${l.warnings}</td><td class="err">${esc(l.last_error)}</td><td style="text-align:left">${esc(l.rel.split('/').pop())}</td></tr>`).join('') + '</tbody></table></div>';
}
async function renderQcLog() {
  let d;
  try { d = await api(`/api/qc/logtext?bids=${encodeURIComponent(QC.bids)}&rel=${encodeURIComponent(QC.log.rel)}&issues=${QC.log.onlyIssues ? 1 : 0}`); }
  catch (e) { $('qc-grid').innerHTML = `<p class="count">${esc(e.message)}</p>`; return; }
  const r = renderLines(d.lines.map(([n, t, k]) => [n, t, k]), false);
  $('qc-count').textContent = '';
  $('qc-grid').innerHTML = `<div class="qc-loghead"><button id="qc-log-back">Back to the logs</button><b>${esc(QC.log.rel)}</b><span class="log-counts">${countsHtml(d.errors, d.warnings)}</span>` +
    `<label class="check"><input type="checkbox" id="qc-log-issues" ${QC.log.onlyIssues ? 'checked' : ''}> errors and warnings only (2 lines of context)</label><button id="qc-log-next">Next error</button></div>` +
    `<pre class="qc-log log" id="qc-log-text">${r.html || '(nothing to show)'}</pre>`;
}
let LB = -1;
function showLightbox(i) { LB = i; const r = QC_VISIBLE[i]; if (!r) return; $('lb-img').src = qcUrl(r); $('lb-cap').textContent = `${r}   (${i + 1} of ${QC_VISIBLE.length}; arrow keys: previous / next)`; $('lightbox').hidden = false; }
$('qc-grid').addEventListener('click', (e) => { const c = e.target.closest('.qc-card'); if (c) showLightbox(Number(c.dataset.i)); });
$('lightbox').addEventListener('click', () => { $('lightbox').hidden = true; LB = -1; });
document.addEventListener('keydown', (e) => {
  if ($('lightbox').hidden) return;
  if (e.key === 'Escape') { $('lightbox').hidden = true; LB = -1; }
  if (e.key === 'ArrowRight') showLightbox(Math.min(LB + 1, QC_VISIBLE.length - 1));
  if (e.key === 'ArrowLeft') showLightbox(Math.max(LB - 1, 0));
});

init();
