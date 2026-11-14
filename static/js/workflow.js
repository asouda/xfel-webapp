(() => {
  const form = document.getElementById("job-form");
  if (!form) return;
  const errorBox = document.getElementById("form-error");
  const runPanel = document.getElementById("run-panel");
  const runState = document.getElementById("run-state");
  const logBox = document.getElementById("log");
  const result = document.getElementById("result");
  const figure = document.getElementById("figure");
  const viewer = document.getElementById("viewer");
  const viewerStatus = document.getElementById("viewer-status");
  const viewerReset = document.getElementById("viewer-reset");
  const viewerSave = document.getElementById("viewer-save");
  const viewerPymol = document.getElementById("viewer-pymol");
  const viewerOrientation = document.getElementById("viewer-orientation");
  const editBond = document.getElementById("edit-bond");
  const editDistance = document.getElementById("edit-distance");
  const editDelete = document.getElementById("edit-delete");
  const editLabel = document.getElementById("edit-label");
  const editDistanceLabels = document.getElementById("edit-distance-labels");
  const surfaceSolid = document.getElementById("surface-solid");
  const surfaceMesh = document.getElementById("surface-mesh");
  const surfaceTransparency = document.getElementById("surface-transparency");
  const surfaceTransparencyValue = document.getElementById("surface-transparency-value");
  const surfaceSigma = document.getElementById("surface-sigma");
  const surfaceSigmaValue = document.getElementById("surface-sigma-value");
  const surfaceSigmaMin = document.getElementById("surface-sigma-min");
  const surfaceSigmaMax = document.getElementById("surface-sigma-max");
  const surfaceSigmaLabel = document.getElementById("surface-sigma-label");
  const surfaceSigmaUnit = document.getElementById("surface-sigma-unit");
  const surfaceSigmaScale = document.getElementById("surface-sigma-scale");
  const surfaceSigmaSwatch = document.getElementById("surface-sigma-swatch");
  const surfaceSigmaNote = document.getElementById("viewer-sigma-note");
  const colorSigned = document.getElementById("color-signed");
  const colorTemperature = document.getElementById("color-temperature");
  const viewerLight = document.getElementById("viewer-light");
  const lightSiteButton = document.getElementById("light-site");
  const lightColorButtons = Array.from(document.querySelectorAll("[data-light-color]"));
  const lightTransparencyInput = document.getElementById("light-transparency");
  const lightTransparencyValue = document.getElementById("light-transparency-value");
  const viewerTitle = document.getElementById("viewer-title");
  const viewerTitleDefault = viewerTitle ? viewerTitle.textContent : "Interactive 3D scene";
  const caption = document.getElementById("caption");
  const fileList = document.getElementById("file-list");
  const runButton = document.getElementById("run-btn");
  const continueButton = document.getElementById("continue-btn");
  const lightHint = document.getElementById("light-hint");
  const lightHintDefault = lightHint ? lightHint.textContent : "";
  const scanButton = document.getElementById("scan-btn");
  const scanPanel = document.getElementById("scan-panel");
  const scanBody = document.getElementById("scan-body");
  const scanSummary = document.getElementById("scan-summary");
  let stage = null;
  let jobId = null;
  let viewerToken = 0;
  let resultShown = false;
  let shownView = null;
  let lightSiteName = "";
  let structureComponent = null;
  let volumeComponent = null;
  let lightComponent = null;
  let lightRepresentation = null;
  let lightAtomCount = 0;
  let lightSiteDescription = "";
  let showLightSite = false;
  let lightColor = "yellow";
  let lightTransparency = 40;
  let surfaceSelection = "";
  let carveRadius = 0;
  let surfaceStyle = "solid";
  let surfaceOpacity = 0.55;
  let sigmaMin = -5;
  let sigmaMax = 5;
  let surfaceSigmaLevel = 3;
  let colorScale = "signed";
  let contourIsSigma = true;
  let sigmaTimer = 0;
  let surfaceGeneration = 0;
  const SigmaColor = window.XFELSigma;
  let distanceComponent = null;
  let markerComponent = null;
  let distancePairs = [];
  let lightDistancePairs = [];
  let labeledAtoms = [];
  let lightLabeledAtoms = [];
  let labelRepresentation = null;
  let lightLabelRepresentation = null;
  let distanceLabelRepresentation = null;
  let lightDistanceLabelRepresentation = null;
  let lightExtraBonds = [];
  let lightRemovedBonds = [];
  let lightBondScene = null;
  let showDistanceLabels = false;
  let siteResidue = null;
  let neighbourResidues = [];
  let modelStick = null;
  let editMode = null;
  let pendingPick = null;
  let sceneNote = "";

  const LABEL_COLOR = "#ff7a18";
  const LABEL_EDGE = "#2a160c";

  const mapCenter = document.getElementById("map-center");
  const mapCenterXyz = document.getElementById("map-center-xyz");

  function syncMapCenter() {
    if (!mapCenter || !mapCenterXyz) return;
    const cartesian = mapCenter.value === "cartesian";
    mapCenterXyz.hidden = !cartesian;
    mapCenterXyz.querySelectorAll("input").forEach((input) => {
      input.required = cartesian;
      input.disabled = !cartesian;
    });
  }

  if (mapCenter) {
    mapCenter.addEventListener("change", syncMapCenter);
    syncMapCenter();
  }

  function showError(message) {
    errorBox.hidden = false;
    errorBox.textContent = message;
  }

  document.querySelectorAll('input[type="file"]').forEach((input) => {
    input.addEventListener("change", () => {
      const span = input.parentElement.querySelector("span");
      span.textContent = input.files[0] ? input.files[0].name : "Choose file";
    });
  });

  const modelNameInput = form.querySelector('[name="dark_id"], [name="pdb_id"]');
  const referenceNameInput = form.querySelector('[name="light_id"]');
  const outputDirInput = form.querySelector('[name="output_dir"]');
  const harmonicAdded = document.getElementById("harmonic-added");
  const harmonicSigma = document.getElementById("harmonic-sigma");

  function coordinateBasename(filename) {
    const base = String(filename || "").split(/[/\\]/).pop();
    if (!base) return "";
    return base.replace(/\.gz$/i, "").replace(/\.(pdb|ent|cif|mmcif)$/i, "");
  }

  function syncHarmonic() {
    if (!harmonicAdded || !harmonicSigma) return;
    harmonicSigma.disabled = !harmonicAdded.checked;
  }

  if (harmonicAdded) {
    harmonicAdded.addEventListener("change", syncHarmonic);
    syncHarmonic();
  }

  const initialOutput = outputDirInput ? outputDirInput.value.replace(/\/+$/, "") : "";
  const initialModelName = modelNameInput ? modelNameInput.value.trim() : "";
  const outputRoot = initialModelName && initialOutput.endsWith(`/${initialModelName}`)
    ? initialOutput.slice(0, -(initialModelName.length + 1))
    : initialOutput;
  let outputAuto = initialOutput;
  let outputTouched = false;

  function syncOutputDir() {
    if (!outputDirInput || !modelNameInput || outputTouched) return;
    const name = modelNameInput.value.trim().replace(/[\\/]/g, "");
    const next = name ? `${outputRoot}/${name}` : outputRoot;
    outputDirInput.value = next;
    outputAuto = next;
  }

  function useCoordinateName(filename, input) {
    const name = coordinateBasename(filename);
    if (!name || !input) return;
    input.value = name;
    if (input === modelNameInput) syncOutputDir();
  }

  function bindCoordinateName(fileKey, input) {
    if (!input) return;
    const file = form.querySelector(`[name="file_${fileKey}"]`);
    const path = form.querySelector(`[name="path_${fileKey}"]`);
    const update = () => {
      const fromFile = file && file.files[0] ? file.files[0].name : "";
      const fromPath = path && path.value.trim() ? path.value.trim() : "";
      if (fromFile || fromPath) useCoordinateName(fromFile || fromPath, input);
    };
    if (file) file.addEventListener("change", update);
    if (path) path.addEventListener("change", update);
  }

  bindCoordinateName(form.querySelector('[name="pdb_id"]') ? "pdb" : "dark_pdb", modelNameInput);
  bindCoordinateName("light_pdb", referenceNameInput);
  if (modelNameInput) modelNameInput.addEventListener("input", syncOutputDir);
  if (outputDirInput) {
    outputDirInput.addEventListener("input", () => {
      if (outputDirInput.value !== outputAuto) outputTouched = true;
    });
    syncOutputDir();
  }

  scanButton.addEventListener("click", async () => {
    const key = window.XFEL_SCAN_KEY;
    const data = new FormData();
    const file = form.querySelector(`[name="file_${key}"]`);
    const path = form.querySelector(`[name="path_${key}"]`);
    if (file && file.files[0]) data.append("file", file.files[0]);
    if (path) data.append("path", path.value);
    scanButton.disabled = true;
    try {
      const response = await fetch("/api/scan", { method: "POST", body: data });
      const payload = await response.json();
      if (!response.ok) {
        showError(payload.error || "Could not read that coordinate file.");
        return;
      }
      errorBox.hidden = true;
      scanSummary.textContent = `${payload.shown} residues listed from ${payload.file} (${payload.total} non-water hetero residues in the file; lipids and pigments are hidden).`;
      scanBody.replaceChildren();
      payload.rows.forEach((row) => {
        const tr = document.createElement("tr");
        [row.chain, row.resi, row.resname, row.natoms, row.occ, row.atoms].forEach((value) => {
          const td = document.createElement("td");
          td.textContent = value;
          tr.append(td);
        });
        tr.addEventListener("click", () => {
          form.querySelector('[name="site_chain"]').value = row.chain.trim();
          form.querySelector('[name="site_resi"]').value = row.resi;
          form.querySelector('[name="site_chain"]').focus();
        });
        scanBody.append(tr);
      });
      scanPanel.hidden = false;
      scanPanel.scrollIntoView({ behavior: "smooth", block: "start" });
    } catch (err) {
      showError("Could not reach the scan service.");
    } finally {
      scanButton.disabled = false;
    }
  });

  form.addEventListener("submit", async (event) => {
    event.preventDefault();
    errorBox.hidden = true;
    result.hidden = true;
    jobId = null;
    resultShown = false;
    shownView = null;
    setViewerTitle(null);
    if (continueButton) continueButton.hidden = true;
    disposeViewer();
    logBox.textContent = "";
    runPanel.hidden = false;
    runState.textContent = "Starting";
    runButton.disabled = true;
    runPanel.scrollIntoView({ behavior: "smooth", block: "start" });
    try {
      const response = await fetch(form.action, { method: "POST", body: new FormData(form) });
      const payload = await response.json();
      if (!response.ok) {
        runState.textContent = "Not started";
        showError(payload.error || "The form could not be submitted.");
        runButton.disabled = false;
        return;
      }
      jobId = payload.id;
      poll(payload.id, 0);
    } catch (err) {
      runState.textContent = "Not started";
      showError("Could not reach the calculation service.");
      runButton.disabled = false;
    }
  });

  if (continueButton) {
    continueButton.addEventListener("click", async () => {
      if (!jobId) return;
      continueButton.disabled = true;
      try {
        const response = await fetch(`/api/jobs/${jobId}/continue`, { method: "POST" });
        const payload = await response.json();
        if (!response.ok) {
          showError(payload.error || "The calculation could not be continued.");
          continueButton.disabled = false;
        }
      } catch (err) {
        showError("Could not reach the calculation service.");
        continueButton.disabled = false;
      }
    });
  }

  async function poll(id, offset) {
    let response;
    try {
      response = await fetch(`/api/jobs/${id}?offset=${offset}`);
    } catch (err) {
      runState.textContent = "Lost contact with the job";
      runButton.disabled = false;
      return;
    }
    if (!response.ok) {
      runState.textContent = "Lost contact with the job";
      runButton.disabled = false;
      return;
    }
    const payload = await response.json();
    if (payload.log) {
      const nearBottom = logBox.scrollHeight - logBox.scrollTop - logBox.clientHeight < 40;
      logBox.textContent += payload.log;
      if (nearBottom) logBox.scrollTop = logBox.scrollHeight;
    }
    runState.textContent = statusLabel(payload);
    if (continueButton) {
      const waiting = payload.status === "paused";
      continueButton.hidden = !waiting;
      continueButton.disabled = !waiting;
    }
    const viewId = payload.result && payload.result.view_id;
    if (payload.result && viewId && viewId !== shownView) {
      shownView = viewId;
      resultShown = true;
      showResult(payload.result);
      if (payload.result.checkpoint) result.scrollIntoView({ behavior: "smooth", block: "start" });
    }
    if (payload.status === "succeeded" || payload.status === "failed") {
      runButton.disabled = false;
      if (payload.status === "failed") showError(payload.error || "The calculation failed. The log above has the last CCP4 or PyMOL output.");
      if (payload.result && !resultShown) {
        resultShown = true;
        showResult(payload.result);
      }
      if (payload.log && payload.log.length >= 120000) {
        poll(id, payload.offset);
      }
      return;
    }
    window.setTimeout(() => poll(id, payload.offset), 1500);
  }

  function statusLabel(payload) {
    const step = Number(payload.step);
    const numbered = Number.isInteger(step) && step >= 1;
    if (payload.status === "running" && numbered) return `Running step ${step}`;
    if (payload.status === "paused" && numbered) return `Paused for inspection after step ${step}`;
    const labels = {
      queued: "Queued",
      running: "Running",
      paused: "Paused for inspection",
      succeeded: "Finished",
      failed: "Failed",
    };
    return labels[payload.status] || payload.status;
  }

  function appendFile(role, path) {
    const wrap = document.createElement("div");
    const dt = document.createElement("dt");
    const dd = document.createElement("dd");
    dt.textContent = role;
    dd.textContent = path;
    wrap.append(dt, dd);
    fileList.append(wrap);
  }

  function setViewerTitle(data) {
    if (!viewerTitle) return;
    const step = Number(data && data.iteration);
    viewerTitle.textContent = Number.isInteger(step) && step >= 1
      ? `Reconstruction Step: ${step}`
      : viewerTitleDefault;
  }

  function showReconstructionSteps(steps) {
    const box = document.getElementById("reconstruction-steps");
    const host = document.getElementById("reconstruction-steps-table");
    if (!box || !host) return;
    host.replaceChildren();
    if (!Array.isArray(steps) || !steps.length) {
      box.hidden = true;
      return;
    }
    const dash = "—";
    const fixed = (value, digits) => (
      typeof value === "number" && Number.isFinite(value) ? value.toFixed(digits) : dash
    );
    const coordinates = (value) => (
      Array.isArray(value) && value.length === 3
        ? value.map((item) => fixed(item, 3)).join(", ")
        : dash
    );
    const columns = [
      ["Step", ""],
      ["R-factor", "num"],
      ["R-free", "num"],
      ["Added atom", ""],
      ["Coordinates (Å)", "num"],
      ["Distance to reference (Å)", "num"],
      ["Voxel density", "num"],
      ["Note", "note"],
    ];
    const table = document.createElement("table");
    const header = document.createElement("tr");
    columns.forEach(([label, kind]) => {
      const cell = document.createElement("th");
      cell.textContent = label;
      if (kind) cell.className = kind;
      cell.scope = "col";
      header.append(cell);
    });
    const thead = document.createElement("thead");
    thead.append(header);
    const tbody = document.createElement("tbody");
    steps.forEach((step) => {
      const values = [
        [String(step.step), ""],
        [fixed(step.r_factor, 4), "num"],
        [fixed(step.r_free, 4), "num"],
        [step.atom ? String(step.atom) : dash, ""],
        [coordinates(step.coordinates), "num"],
        [fixed(step.distance, 2), "num"],
        [fixed(step.density, 3), "num"],
        [step.note ? String(step.note) : "", "note"],
      ];
      const row = document.createElement("tr");
      values.forEach(([value, kind]) => {
        const cell = document.createElement("td");
        cell.textContent = value;
        if (kind) cell.className = kind;
        row.append(cell);
      });
      tbody.append(row);
    });
    table.append(thead, tbody);
    host.append(table);
    box.hidden = false;
  }

  function showResult(data) {
    setViewerTitle(data);
    figure.src = `${data.figure_url}?t=${Date.now()}`;
    caption.textContent = data.caption;
    showReconstructionSteps(data.steps);
    loadViewer(data.viewer);
    fileList.replaceChildren();
    const rows = [
      ["Output directory", data.output_dir],
      ["Extended map", data.final_map],
      ["Site fragment", data.fragment],
      ["PyMOL script", data.pymol_script],
    ];
    rows.forEach(([role, path]) => appendFile(role, path));
    result.hidden = false;
  }

  function atomSelection(spec) {
    let sele = `.${spec.atom} and ${spec.resi}:${spec.chain}`;
    if (spec.alt) sele += ` and %${spec.alt}`;
    return sele;
  }

  function atomIndices(structure, spec) {
    const indices = structure.getAtomIndices(new NGL.Selection(atomSelection(spec)));
    return indices || [];
  }

  function applyScriptBonds(component, scene) {
    const bonds = scene.bonds || [];
    if (!bonds.length || !scene.site) return;
    const structure = component.structure;
    const store = structure.bondStore;
    const site = new Set(structure.getAtomIndices(new NGL.Selection(`${scene.site.resi}:${scene.site.chain}`)) || []);
    const kept = [];
    for (let index = 0; index < store.count; index += 1) {
      const left = store.atomIndex1[index];
      const right = store.atomIndex2[index];
      if (site.has(left) || site.has(right)) continue;
      kept.push([left, right, store.bondOrder[index] || 1]);
    }
    store.count = 0;
    const leftAtom = structure.getAtomProxy();
    const rightAtom = structure.getAtomProxy();
    kept.forEach(([left, right, order]) => {
      leftAtom.index = left;
      rightAtom.index = right;
      store.addBond(leftAtom, rightAtom, order);
    });
    bonds.forEach((bond) => {
      atomIndices(structure, bond.a).forEach((left) => {
        atomIndices(structure, bond.b).forEach((right) => {
          if (left === right) return;
          leftAtom.index = left;
          rightAtom.index = right;
          store.addBond(leftAtom, rightAtom, 1);
        });
      });
    });
    structure.finalizeBonds();
  }

  function pairKey(left, right) {
    return left < right ? [left, right] : [right, left];
  }

  function componentFor(source) {
    return source === "light" ? lightComponent : structureComponent;
  }

  function structureFor(source) {
    const component = componentFor(source);
    return component ? component.structure : null;
  }

  function describePick(pick) {
    const structure = structureFor(pick.source);
    const atom = structure.getAtomProxy();
    atom.index = pick.index;
    const alt = atom.altloc && atom.altloc.trim() ? ` alt ${atom.altloc.trim()}` : "";
    const chain = (atom.chainname || "").trim();
    const where = pick.source === "light" ? " light" : "";
    return `${atom.atomname.trim()} ${chain}/${atom.resno}${alt}${where}`;
  }

  function samePick(left, right) {
    return !!left && !!right && left.source === right.source && left.index === right.index;
  }

  function bondSlotIn(structure, left, right) {
    const [first, second] = pairKey(left, right);
    const store = structure.bondStore;
    for (let index = 0; index < store.count; index += 1) {
      if (store.atomIndex1[index] === first && store.atomIndex2[index] === second) return index;
    }
    return -1;
  }

  function bondSlot(left, right) {
    return bondSlotIn(structureComponent.structure, left, right);
  }

  function addBondIn(structure, left, right) {
    const first = structure.getAtomProxy();
    const second = structure.getAtomProxy();
    first.index = left;
    second.index = right;
    structure.bondStore.addBond(first, second, 1);
  }

  function removeBondIn(structure, left, right) {
    const slot = bondSlotIn(structure, left, right);
    if (slot < 0) return false;
    const store = structure.bondStore;
    const last = store.count - 1;
    if (slot !== last) {
      store.atomIndex1[slot] = store.atomIndex1[last];
      store.atomIndex2[slot] = store.atomIndex2[last];
      store.bondOrder[slot] = store.bondOrder[last];
    }
    store.count = last;
    return true;
  }

  function distanceSlot(left, right) {
    const [first, second] = pairKey(left, right);
    return distancePairs.findIndex((pair) => pair[0] === first && pair[1] === second);
  }

  function refreshComponentSticks(component) {
    if (!component) return;
    const structure = component.structure;
    structure.finalizeBonds();
    // The stick view keeps the bond list it had when it was built. A new bond
    // is stored on the structure, but it stays invisible until that list is rebuilt.
    const views = [];
    if (component.structureView) views.push(component.structureView);
    (component.reprList || []).forEach((entry) => {
      const repr = entry.repr || entry;
      if (repr.structureView) views.push(repr.structureView);
    });
    views.forEach((view) => view.refresh());
    component.rebuildRepresentations();
  }

  function refreshSticks() {
    refreshComponentSticks(structureComponent);
    syncLightBonds();
  }

  function atomIdentity(atom) {
    return {
      chain: (atom.chainname || "").trim(),
      resi: atom.resno,
      name: atom.atomname.trim(),
      alt: (atom.altloc || "").trim(),
    };
  }

  function selectedAtomSet(structure, selection) {
    return new Set(structure.getAtomIndices(new NGL.Selection(selection)) || []);
  }

  function atomsMatching(structure, selected, spec) {
    const atom = structure.getAtomProxy();
    const found = [];
    selected.forEach((index) => {
      atom.index = index;
      const identity = atomIdentity(atom);
      if (identity.chain !== spec.chain || identity.resi !== spec.resi || identity.name !== spec.name) return;
      if (spec.alt && identity.alt !== spec.alt) return;
      found.push(identity.alt ? { index, alt: identity.alt } : { index, alt: "" });
    });
    return found;
  }

  // Same-altloc copies stay paired. An atom with no altloc bonds to each
  // conformer of the other atom, which is how a single metal meets split ligands.
  function pairConformers(left, right, order, into) {
    left.forEach((first) => {
      right.forEach((second) => {
        if (first.alt && second.alt && first.alt !== second.alt) return;
        if (first.index === second.index) return;
        const low = Math.min(first.index, second.index);
        const high = Math.max(first.index, second.index);
        const key = `${low}:${high}`;
        if (into.seen.has(key)) return;
        into.seen.add(key);
        into.pairs.push([low, high, order]);
      });
    });
  }

  function selectionBondSpecs(structure, selected) {
    const store = structure.bondStore;
    const atom = structure.getAtomProxy();
    const specs = [];
    for (let index = 0; index < store.count; index += 1) {
      const left = store.atomIndex1[index];
      const right = store.atomIndex2[index];
      if (!selected.has(left) || !selected.has(right)) continue;
      atom.index = left;
      const a = atomIdentity(atom);
      atom.index = right;
      const b = atomIdentity(atom);
      specs.push({ a, b, order: store.bondOrder[index] || 1 });
    }
    return specs;
  }

  function replaceSelectionBonds(structure, selected, pairs) {
    const store = structure.bondStore;
    const kept = [];
    for (let index = 0; index < store.count; index += 1) {
      const left = store.atomIndex1[index];
      const right = store.atomIndex2[index];
      if (selected.has(left) && selected.has(right)) continue;
      kept.push([left, right, store.bondOrder[index] || 1]);
    }
    store.count = 0;
    const leftAtom = structure.getAtomProxy();
    const rightAtom = structure.getAtomProxy();
    kept.concat(pairs).forEach(([left, right, order]) => {
      leftAtom.index = left;
      rightAtom.index = right;
      store.addBond(leftAtom, rightAtom, order);
    });
    structure.finalizeBonds();
  }

  function copySelectionBonds(source, target, selection) {
    if (!selection) return;
    const sourceSelected = selectedAtomSet(source, selection);
    const targetSelected = selectedAtomSet(target, selection);
    const specs = selectionBondSpecs(source, sourceSelected);
    const into = { seen: new Set(), pairs: [] };
    specs.forEach((spec) => {
      const left = atomsMatching(target, targetSelected, spec.a);
      const right = atomsMatching(target, targetSelected, spec.b);
      pairConformers(left, right, spec.order, into);
    });
    replaceSelectionBonds(target, targetSelected, into.pairs);
  }

  function identityKey(identity) {
    return `${identity.chain}|${identity.resi}|${identity.name}|${identity.alt}`;
  }

  function pairIdentityKey(left, right) {
    const a = identityKey(left);
    const b = identityKey(right);
    return a < b ? `${a}~${b}` : `${b}~${a}`;
  }

  function identityAt(source, index) {
    const atom = structureFor(source).getAtomProxy();
    atom.index = index;
    return atomIdentity(atom);
  }

  function findAtom(source, identity) {
    const structure = structureFor(source);
    if (!structure || !surfaceSelection) return null;
    const selected = selectedAtomSet(structure, surfaceSelection);
    const atom = structure.getAtomProxy();
    let found = null;
    selected.forEach((index) => {
      if (found != null) return;
      atom.index = index;
      const id = atomIdentity(atom);
      if (id.chain === identity.chain && id.resi === identity.resi && id.name === identity.name && id.alt === identity.alt) {
        found = index;
      }
    });
    return found;
  }

  function omitIndexFor(pick) {
    if (pick.source === "model") return pick.index;
    return findAtom("model", identityAt("light", pick.index));
  }

  function applyLightExtraBonds() {
    if (!lightComponent) return;
    const structure = lightComponent.structure;
    lightExtraBonds.forEach((bond) => {
      const left = findAtom("light", bond.a);
      const right = findAtom("light", bond.b);
      if (left == null || right == null || left === right) return;
      if (bondSlotIn(structure, left, right) >= 0) return;
      addBondIn(structure, left, right);
    });
  }

  function forgetLightRemoval(left, right) {
    const key = pairIdentityKey(identityAt("light", left), identityAt("light", right));
    lightRemovedBonds = lightRemovedBonds.filter((bond) => pairIdentityKey(bond.a, bond.b) !== key);
  }

  function rememberLightRemoval(left, right) {
    forgetLightRemoval(left, right);
    lightRemovedBonds.push({ a: identityAt("light", left), b: identityAt("light", right) });
  }

  function applyLightRemovals() {
    lightRemovedBonds.forEach((bond) => {
      const left = findAtom("light", bond.a);
      const right = findAtom("light", bond.b);
      if (left == null || right == null || left === right) return;
      removeBondIn(lightComponent.structure, left, right);
    });
  }

  function syncLightBonds() {
    if (!lightComponent || !surfaceSelection) return;
    if (lightBondScene) {
      applyScriptBonds(lightComponent, lightBondScene);
      applyLightRemovals();
      applyLightExtraBonds();
      refreshComponentSticks(lightComponent);
      return;
    }
    if (!structureComponent) return;
    copySelectionBonds(structureComponent.structure, lightComponent.structure, surfaceSelection);
    applyLightExtraBonds();
    refreshComponentSticks(lightComponent);
  }

  function dropLightExtra(left, right) {
    const key = pairIdentityKey(identityAt("light", left), identityAt("light", right));
    const before = lightExtraBonds.length;
    lightExtraBonds = lightExtraBonds.filter((bond) => pairIdentityKey(bond.a, bond.b) !== key);
    return lightExtraBonds.length !== before;
  }

  function addBondBetween(left, right) {
    if (bondSlot(left, right) >= 0) return false;
    addBondIn(structureComponent.structure, left, right);
    refreshSticks();
    return true;
  }

  function removeBondBetween(left, right) {
    if (!removeBondIn(structureComponent.structure, left, right)) return false;
    refreshSticks();
    return true;
  }

  function addBondFromPicks(left, right) {
    if (lightBondScene && left.source === "light" && right.source === "light") {
      if (bondSlotIn(lightComponent.structure, left.index, right.index) >= 0) return false;
      forgetLightRemoval(left.index, right.index);
      lightExtraBonds.push({ a: identityAt("light", left.index), b: identityAt("light", right.index) });
      const omitLeft = omitIndexFor(left);
      const omitRight = omitIndexFor(right);
      if (omitLeft != null && omitRight != null) return addBondBetween(omitLeft, omitRight);
      syncLightBonds();
      return true;
    }
    const omitLeft = omitIndexFor(left);
    const omitRight = omitIndexFor(right);
    if (omitLeft != null && omitRight != null) return addBondBetween(omitLeft, omitRight);
    if (left.source !== "light" || right.source !== "light") return null;
    if (bondSlotIn(lightComponent.structure, left.index, right.index) >= 0) return false;
    lightExtraBonds.push({ a: identityAt("light", left.index), b: identityAt("light", right.index) });
    syncLightBonds();
    return true;
  }

  function removeBondFromPicks(left, right) {
    let removed = false;
    if (left.source === "light" && right.source === "light") {
      if (lightBondScene) rememberLightRemoval(left.index, right.index);
      if (dropLightExtra(left.index, right.index)) removed = true;
    }
    const omitLeft = omitIndexFor(left);
    const omitRight = omitIndexFor(right);
    if (omitLeft != null && omitRight != null && bondSlot(omitLeft, omitRight) >= 0) {
      return removeBondBetween(omitLeft, omitRight) || removed;
    }
    if (left.source === "light" && right.source === "light") {
      if (removeBondIn(lightComponent.structure, left.index, right.index)) removed = true;
      if (removed) syncLightBonds();
    }
    return removed;
  }

  function addDashed(shape, start, end, name) {
    const delta = [end[0] - start[0], end[1] - start[1], end[2] - start[2]];
    const length = Math.hypot(delta[0], delta[1], delta[2]);
    if (!length) return;
    const color = [1, 0.92, 0.1];
    const radius = 0.06;
    const dash = 0.18;
    const gap = 0.12;
    if (length <= dash + gap) {
      shape.addCylinder(start, end, color, radius, name);
      return;
    }
    const step = [delta[0] / length, delta[1] / length, delta[2] / length];
    for (let traveled = 0; traveled + dash <= length; traveled += dash + gap) {
      shape.addCylinder(
        [start[0] + step[0] * traveled, start[1] + step[1] * traveled, start[2] + step[2] * traveled],
        [start[0] + step[0] * (traveled + dash), start[1] + step[1] * (traveled + dash), start[2] + step[2] * (traveled + dash)],
        color,
        radius,
        name
      );
    }
  }

  function captureView() {
    if (!stage) return null;
    const viewer = stage.viewer;
    return {
      position: viewer.translationGroup.position.clone(),
      distance: viewer.cameraDistance,
      quaternion: viewer.rotationGroup.quaternion.clone(),
    };
  }

  function restoreView(saved) {
    if (!saved || !stage) return;
    const viewer = stage.viewer;
    viewer.translationGroup.position.copy(saved.position);
    viewer.rotationGroup.quaternion.copy(saved.quaternion);
    viewer.cameraDistance = saved.distance;
    if (viewer.updateZoom) viewer.updateZoom();
  }

  function appendDistances(shape, structure, pairs, source) {
    if (!structure || !pairs.length) return;
    const atom = structure.getAtomProxy();
    pairs.forEach(([left, right]) => {
      if (!pairResidueShown(source, [left, right])) return;
      atom.index = left;
      const start = [atom.x, atom.y, atom.z];
      atom.index = right;
      const end = [atom.x, atom.y, atom.z];
      addDashed(shape, start, end, `dist:${source}:${left}:${right}`);
    });
  }

  function redrawDistances() {
    const saved = captureView();
    if (distanceComponent) {
      stage.removeComponent(distanceComponent);
      distanceComponent = null;
    }
    const any = (distancePairs.length && structureComponent) || (lightDistancePairs.length && lightComponent);
    if (any) {
      const shape = new NGL.Shape("script-distances");
      appendDistances(shape, structureFor("model"), distancePairs, "model");
      appendDistances(shape, structureFor("light"), lightDistancePairs, "light");
      distanceComponent = stage.addComponentFromObject(shape);
      distanceComponent.addRepresentation("buffer", { disablePicking: true });
    }
    restoreView(saved);
    syncDistanceLabels();
  }

  function mountDistanceLabels(component, pairs) {
    const source = component === lightComponent ? "light" : "model";
    const visible = pairs.filter((pair) => pairResidueShown(source, pair));
    if (!component || !showDistanceLabels || !visible.length) return null;
    const representation = component.addRepresentation("distance", {
      atomPair: visible.map((pair) => pair.slice()),
      labelUnit: "",
      labelColor: LABEL_COLOR,
      labelBorder: true,
      labelBorderColor: LABEL_EDGE,
      labelBorderWidth: 0.18,
      labelSize: 1.35,
      labelFontWeight: "normal",
      color: "yellow",
      lineVisible: false,
      disablePicking: true,
    });
    const distanceLabels = representation.repr || representation;
    distanceLabels.lineVisible = false;
    if (distanceLabels.distanceBuffer) distanceLabels.distanceBuffer.setVisibility(false);
    return representation;
  }

  function syncDistanceLabels() {
    const saved = captureView();
    if (distanceLabelRepresentation && structureComponent) {
      structureComponent.removeRepresentation(distanceLabelRepresentation);
    }
    distanceLabelRepresentation = null;
    if (lightDistanceLabelRepresentation && lightComponent) {
      lightComponent.removeRepresentation(lightDistanceLabelRepresentation);
    }
    lightDistanceLabelRepresentation = null;
    distanceLabelRepresentation = mountDistanceLabels(structureComponent, distancePairs);
    lightDistanceLabelRepresentation = mountDistanceLabels(lightComponent, lightDistancePairs);
    restoreView(saved);
  }

  function atomName(source, index) {
    const atom = structureFor(source).getAtomProxy();
    atom.index = index;
    return atom.atomname.trim();
  }

  function labelParameters(indices, names) {
    return {
      sele: `@${indices.join(",")}`,
      labelType: "text",
      labelText: names,
      labelGrouping: "atom",
      color: LABEL_COLOR,
      showBorder: true,
      borderColor: LABEL_EDGE,
      borderWidth: 0.18,
      fontWeight: "normal",
      radiusType: "size",
      radiusSize: 1.125,
      xOffset: 0.15,
      yOffset: 0.3,
      zOffset: 0.7,
      attachment: "bottom-left",
    };
  }

  function redrawLabels(source) {
    const saved = captureView();
    const component = componentFor(source);
    const indices = (source === "light" ? lightLabeledAtoms : labeledAtoms)
      .filter((index) => atomResidueShown(source, index));
    const current = source === "light" ? lightLabelRepresentation : labelRepresentation;
    if (current && component) component.removeRepresentation(current);
    let next = null;
    if (indices.length && component) {
      const names = {};
      indices.forEach((index) => {
        names[index] = atomName(source, index);
      });
      next = component.addRepresentation("label", labelParameters(indices, names));
    }
    if (source === "light") lightLabelRepresentation = next;
    else labelRepresentation = next;
    restoreView(saved);
  }

  function redrawAtomLabels() {
    redrawLabels("model");
  }

  function toggleAtomLabel(pick) {
    const list = pick.source === "light" ? lightLabeledAtoms : labeledAtoms;
    const slot = list.indexOf(pick.index);
    const name = describePick(pick);
    if (slot >= 0) {
      list.splice(slot, 1);
      redrawLabels(pick.source);
      viewerStatus.textContent = `Removed the ${name} label. ${nextPrompt()}`;
      return;
    }
    list.push(pick.index);
    redrawLabels(pick.source);
    viewerStatus.textContent = `Labeled ${name}. ${nextPrompt()}`;
  }

  function addDistanceBetween(left, right) {
    if (distanceSlot(left, right) >= 0) return false;
    distancePairs.push(pairKey(left, right));
    redrawDistances();
    return true;
  }

  function removeDistanceBetween(left, right) {
    const slot = distanceSlot(left, right);
    if (slot < 0) return false;
    distancePairs.splice(slot, 1);
    redrawDistances();
    return true;
  }

  function clearMarker() {
    if (markerComponent && stage) stage.removeComponent(markerComponent);
    markerComponent = null;
  }

  function lightDistanceSlot(left, right) {
    const [first, second] = pairKey(left, right);
    return lightDistancePairs.findIndex((pair) => pair[0] === first && pair[1] === second);
  }

  function addDistanceFromPicks(left, right) {
    if (left.source !== right.source) return null;
    if (left.source === "light") {
      if (lightDistanceSlot(left.index, right.index) >= 0) return false;
      lightDistancePairs.push(pairKey(left.index, right.index));
      redrawDistances();
      return true;
    }
    return addDistanceBetween(left.index, right.index);
  }

  function removeDistanceFromPicks(left, right) {
    if (left.source !== right.source) return false;
    if (left.source === "light") {
      const slot = lightDistanceSlot(left.index, right.index);
      if (slot < 0) return false;
      lightDistancePairs.splice(slot, 1);
      redrawDistances();
      return true;
    }
    return removeDistanceBetween(left.index, right.index);
  }

  function showMarker(pick) {
    clearMarker();
    const atom = structureFor(pick.source).getAtomProxy();
    atom.index = pick.index;
    const shape = new NGL.Shape("pending-atom");
    const point = [atom.x, atom.y, atom.z];
    const arm = 1.15;
    const color = [0.2, 0.85, 1];
    const name = `pick:${pick.source}:${pick.index}`;
    shape.addCylinder([point[0] - arm, point[1], point[2]], [point[0] + arm, point[1], point[2]], color, 0.05, name);
    shape.addCylinder([point[0], point[1] - arm, point[2]], [point[0], point[1] + arm, point[2]], color, 0.05, name);
    shape.addCylinder([point[0], point[1], point[2] - arm], [point[0], point[1], point[2] + arm], color, 0.05, name);
    markerComponent = stage.addComponentFromObject(shape);
    markerComponent.addRepresentation("buffer");
  }

  function collectScriptDistances(component, scene) {
    const structure = component.structure;
    (scene.distances || []).forEach((distance) => {
      atomIndices(structure, distance.a).forEach((left) => {
        atomIndices(structure, distance.b).forEach((right) => {
          if (left === right || distanceSlot(left, right) >= 0) return;
          distancePairs.push(pairKey(left, right));
        });
      });
    });
    redrawDistances();
  }

  function nextPrompt() {
    if (editMode === "bond") return "Click the first atom of another bond.";
    if (editMode === "distance") return "Click the first atom of another distance.";
    if (editMode === "delete") return "Click a stick, or the first of two atoms.";
    if (editMode === "label") return "Click an atom to add or remove its name.";
    return "";
  }

  function setEditMode(mode) {
    editMode = editMode === mode ? null : mode;
    pendingPick = null;
    clearMarker();
    if (editBond) editBond.setAttribute("aria-pressed", editMode === "bond" ? "true" : "false");
    if (editDistance) editDistance.setAttribute("aria-pressed", editMode === "distance" ? "true" : "false");
    if (editDelete) editDelete.setAttribute("aria-pressed", editMode === "delete" ? "true" : "false");
    if (editLabel) editLabel.setAttribute("aria-pressed", editMode === "label" ? "true" : "false");
    if (!editMode) {
      viewerStatus.textContent = sceneNote;
      return;
    }
    if (editMode === "bond") viewerStatus.textContent = "Click the first atom of the bond.";
    else if (editMode === "distance") viewerStatus.textContent = "Click the first atom of the distance.";
    else if (editMode === "label") viewerStatus.textContent = "Click an atom to show its name. Click a labeled atom again to remove the label.";
    else viewerStatus.textContent = "Click a stick, or click the two atoms that define the bond or distance.";
  }

  function resetEditControls() {
    editMode = null;
    pendingPick = null;
    distancePairs = [];
    lightDistancePairs = [];
    labeledAtoms = [];
    lightLabeledAtoms = [];
    lightExtraBonds = [];
    lightRemovedBonds = [];
    lightBondScene = null;
    labelRepresentation = null;
    lightLabelRepresentation = null;
    distanceLabelRepresentation = null;
    lightDistanceLabelRepresentation = null;
    showDistanceLabels = false;
    distanceComponent = null;
    markerComponent = null;
    structureComponent = null;
    volumeComponent = null;
    [editBond, editDistance, editDelete, editLabel, editDistanceLabels, surfaceSolid, surfaceMesh].forEach((button) => {
      if (!button) return;
      button.disabled = true;
      button.setAttribute("aria-pressed", "false");
    });
    if (surfaceTransparency) surfaceTransparency.disabled = true;
  }

  function repickWithoutDistances(proxy) {
    if (!distanceComponent || !proxy || !proxy.canvasPosition) return null;
    distanceComponent.setVisibility(false);
    const again = stage.pickingControls.pick(proxy.canvasPosition.x, proxy.canvasPosition.y);
    distanceComponent.setVisibility(true);
    return again || null;
  }

  function pickFromProxy(proxy) {
    if (!proxy) return null;
    const named = proxy.cylinder || proxy.sphere;
    if (named && typeof named.name === "string" && named.name.startsWith("pick:")) {
      const parts = named.name.split(":");
      if (parts.length === 3 && (parts[1] === "model" || parts[1] === "light")) {
        const index = Number(parts[2]);
        if (!Number.isNaN(index)) return { source: parts[1], index };
      }
    }
    if (proxy.atom && proxy.component) {
      return { source: proxy.component === lightComponent ? "light" : "model", index: proxy.atom.index };
    }
    if ((editMode === "bond" || editMode === "distance" || editMode === "label") && proxy.bond && proxy.component) {
      const index = proxy.closestBondAtom ? proxy.closestBondAtom.index : proxy.bond.atomIndex1;
      return { source: proxy.component === lightComponent ? "light" : "model", index };
    }
    return null;
  }

  function nearestLabeled(proxy, component, indices) {
    if (!component || !indices.length || !proxy.canvasPosition) return null;
    const controls = stage.viewerControls;
    const atom = component.structure.getAtomProxy();
    let best = null;
    let bestDistance = 36;
    indices.forEach((index) => {
      atom.index = index;
      const point = atom.positionToVector3();
      point.applyMatrix4(component.matrix);
      const place = controls.getPositionOnCanvas(point);
      const distance = Math.hypot(place.x - proxy.canvasPosition.x, place.y - proxy.canvasPosition.y);
      if (distance < bestDistance) {
        bestDistance = distance;
        best = { index, distance };
      }
    });
    return best;
  }

  function labeledAtomNear(proxy) {
    if (!proxy || !proxy.canvasPosition) return null;
    const model = nearestLabeled(proxy, structureComponent, labeledAtoms);
    const light = showLightSite ? nearestLabeled(proxy, lightComponent, lightLabeledAtoms) : null;
    if (model && light) {
      return model.distance <= light.distance
        ? { source: "model", index: model.index }
        : { source: "light", index: light.index };
    }
    if (model) return { source: "model", index: model.index };
    if (light) return { source: "light", index: light.index };
    return null;
  }

  function differentModelsMessage() {
    return "Choose both atoms from the light site, or both from the omit model.";
  }

  function finishPair(left, right) {
    pendingPick = null;
    clearMarker();
    const first = describePick(left);
    const second = describePick(right);
    if (editMode === "bond") {
      const added = addBondFromPicks(left, right);
      if (added === null) viewerStatus.textContent = `${differentModelsMessage()} ${nextPrompt()}`;
      else if (added) viewerStatus.textContent = `Added a bond between ${first} and ${second}. ${nextPrompt()}`;
      else viewerStatus.textContent = `${first} and ${second} already have a bond. ${nextPrompt()}`;
      return;
    }
    if (editMode === "distance") {
      const added = addDistanceFromPicks(left, right);
      if (added === null) viewerStatus.textContent = `${differentModelsMessage()} ${nextPrompt()}`;
      else if (added) viewerStatus.textContent = `Added a distance between ${first} and ${second}. ${nextPrompt()}`;
      else viewerStatus.textContent = `Those atoms already have a distance. ${nextPrompt()}`;
      return;
    }
    const removedBond = removeBondFromPicks(left, right);
    const removedDistance = removeDistanceFromPicks(left, right);
    if (removedBond && removedDistance) viewerStatus.textContent = `Removed the bond and the distance between ${first} and ${second}. ${nextPrompt()}`;
    else if (removedBond) viewerStatus.textContent = `Removed the bond between ${first} and ${second}. ${nextPrompt()}`;
    else if (removedDistance) viewerStatus.textContent = `Removed the distance between ${first} and ${second}. ${nextPrompt()}`;
    else viewerStatus.textContent = `Those atoms have no bond or distance to remove. ${nextPrompt()}`;
  }

  function onStageClick(proxy) {
    if (!editMode || !structureComponent || !proxy) return;
    if (editMode === "delete" && proxy.cylinder && typeof proxy.cylinder.name === "string") {
      const match = /^dist:(model|light):(\d+):(\d+)$/.exec(proxy.cylinder.name);
      if (match) {
        const left = { source: match[1], index: Number(match[2]) };
        const right = { source: match[1], index: Number(match[3]) };
        pendingPick = null;
        clearMarker();
        removeDistanceFromPicks(left, right);
        viewerStatus.textContent = `Removed the distance between ${describePick(left)} and ${describePick(right)}. ${nextPrompt()}`;
        return;
      }
    }
    if (editMode === "delete" && proxy.bond && proxy.component) {
      const source = proxy.component === lightComponent ? "light" : "model";
      const left = { source, index: proxy.bond.atomIndex1 };
      const right = { source, index: proxy.bond.atomIndex2 };
      pendingPick = null;
      clearMarker();
      removeBondFromPicks(left, right);
      viewerStatus.textContent = `Removed the bond between ${describePick(left)} and ${describePick(right)}. ${nextPrompt()}`;
      return;
    }
    let picked = proxy;
    if ((editMode === "bond" || editMode === "distance" || editMode === "label") && proxy.cylinder && /^dist:/.test(proxy.cylinder.name || "")) {
      picked = repickWithoutDistances(proxy) || proxy;
    }
    const pick = pickFromProxy(picked);
    if (editMode === "label") {
      const labeled = pick || labeledAtomNear(proxy);
      if (!labeled || Number.isNaN(labeled.index)) return;
      toggleAtomLabel(labeled);
      return;
    }
    if (!pick || Number.isNaN(pick.index)) return;
    if (pendingPick == null) {
      pendingPick = pick;
      showMarker(pick);
      viewerStatus.textContent = `${describePick(pick)} selected. Click the second atom.`;
      return;
    }
    if (samePick(pendingPick, pick)) {
      viewerStatus.textContent = `${describePick(pick)} is already selected. Click a different atom.`;
      return;
    }
    finishPair(pendingPick, pick);
  }

  function residueSele(residue) {
    return `${residue.resi}:${residue.chain}`;
  }

  function structureSelection(includeHidden) {
    const residues = [];
    if (siteResidue) residues.push(residueSele(siteResidue));
    neighbourResidues.forEach((residue) => {
      if (includeHidden || residue.shown) residues.push(residueSele(residue));
    });
    return residues.join(" or ");
  }

  function atomResidueShown(source, index) {
    const structure = structureFor(source);
    if (!structure || !neighbourResidues.length) return true;
    const atom = structure.getAtomProxy();
    atom.index = index;
    const chain = String(atom.chainname || "").trim();
    const resi = Number(atom.resno);
    const match = neighbourResidues.find((residue) => residue.chain === chain && residue.resi === resi);
    return !match || match.shown;
  }

  function pairResidueShown(source, pair) {
    return atomResidueShown(source, pair[0]) && atomResidueShown(source, pair[1]);
  }

  function neighbourButtonLabel(residue) {
    const same = neighbourResidues.filter((item) => item.label === residue.label).length;
    return same > 1 ? `${residue.label} ${residue.chain}${residue.resi}` : residue.label;
  }

  function syncNeighbourControls() {
    const box = document.getElementById("viewer-neighbours");
    const host = document.getElementById("neighbour-toggles");
    if (!box || !host) return;
    host.replaceChildren();
    if (!neighbourResidues.length) {
      box.hidden = true;
      return;
    }
    box.hidden = false;
    neighbourResidues.forEach((residue) => {
      const button = document.createElement("button");
      button.type = "button";
      button.className = "pill";
      const name = neighbourButtonLabel(residue);
      button.textContent = name;
      const describe = () => `${residue.shown ? "Hide" : "Show"} ${name}, chain ${residue.chain} residue ${residue.resi}`;
      button.setAttribute("aria-pressed", residue.shown ? "true" : "false");
      button.setAttribute("aria-label", describe());
      button.addEventListener("click", () => {
        if (!structureComponent) return;
        residue.shown = !residue.shown;
        button.setAttribute("aria-pressed", residue.shown ? "true" : "false");
        button.setAttribute("aria-label", describe());
        applyNeighbourVisibility();
        viewerStatus.textContent = residue.shown ? `${name} is shown.` : `${name} is hidden.`;
      });
      host.append(button);
    });
  }

  function applyNeighbourVisibility() {
    const sele = structureSelection(false);
    const saved = captureView();
    if (modelStick) modelStick.setSelection(sele);
    if (lightRepresentation) lightRepresentation.setSelection(sele);
    if (pendingPick && !atomResidueShown(pendingPick.source, pendingPick.index)) {
      pendingPick = null;
      clearMarker();
    }
    redrawDistances();
    redrawLabels("model");
    redrawLabels("light");
    whenStageIdle(viewerToken, () => restoreView(saved));
  }

  function lightStatusText() {
    const color = lightColor === "element" ? "by element name" : lightColor;
    if (lightSiteName) {
      if (!lightAtomCount) return `${lightSiteName} has no atoms in ${lightSiteDescription}.`;
      if (showLightSite) return `${lightSiteName} is ball and stick, colored ${color}.`;
      return `${lightSiteName} is hidden. Its color is ${color}.`;
    }
    if (!lightAtomCount) return `The light coordinates have no atoms in ${lightSiteDescription}.`;
    if (showLightSite) return `Light site is ball and stick, colored ${color}.`;
    return `Light site is hidden. Its color is ${color}.`;
  }

  function lightOpacityValue() {
    return 1 - lightTransparency / 100;
  }

  function syncLightControls(available) {
    if (!viewerLight) return;
    viewerLight.hidden = !available;
    if (lightSiteButton) {
      lightSiteButton.disabled = !available;
      lightSiteButton.setAttribute("aria-pressed", showLightSite ? "true" : "false");
    }
    lightColorButtons.forEach((button) => {
      button.disabled = !available;
      button.setAttribute("aria-pressed", button.dataset.lightColor === lightColor ? "true" : "false");
    });
    if (lightTransparencyInput) {
      lightTransparencyInput.disabled = !available;
      lightTransparencyInput.value = String(lightTransparency);
    }
    if (lightTransparencyValue) lightTransparencyValue.textContent = `${lightTransparency}%`;
  }

  function applyLightAppearance() {
    if (!lightRepresentation) return;
    const show = showLightSite && lightAtomCount > 0;
    const opacity = lightOpacityValue();
    // Hide the whole component while it is off. A deposited model is much
    // larger than the site, and a visible component would pull the camera out.
    if (lightComponent) lightComponent.setVisibility(show);
    lightRepresentation.setColor(lightColor);
    lightRepresentation.setParameters({ opacity, depthWrite: opacity >= 0.999 });
    lightRepresentation.setVisibility(show);
    if (stage) stage.viewer.requestRender();
  }

  function modelPoints(structureComponent, selection) {
    const points = [];
    structureComponent.structure.eachAtom((atom) => {
      points.push(atom.positionToVector3().applyMatrix4(structureComponent.matrix));
    }, new NGL.Selection(selection));
    return points;
  }

  function pointsInVolume(points, volumeComponent) {
    const inverse = volumeComponent.matrix.clone().invert();
    return points.map((point) => point.clone().applyMatrix4(inverse));
  }

  // PyMOL `carve` keeps isosurface vertices within this radius of the model.
  // The index is groups of 3 for a solid surface and groups of 2 for a mesh.
  function carvedIndex(position, indexArray, atoms, radius, stride, count) {
    const radius2 = radius * radius;
    const near = new Uint8Array(position.length / 3);
    for (let vertex = 0; vertex < near.length; vertex += 1) {
      const x = position[vertex * 3];
      const y = position[vertex * 3 + 1];
      const z = position[vertex * 3 + 2];
      for (let a = 0; a < atoms.length; a += 1) {
        const atom = atoms[a];
        const dx = x - atom.x;
        const dy = y - atom.y;
        const dz = z - atom.z;
        if (dx * dx + dy * dy + dz * dz <= radius2) {
          near[vertex] = 1;
          break;
        }
      }
    }
    const limit = count == null ? indexArray.length : count;
    const kept = [];
    for (let i = 0; i + stride <= limit; i += stride) {
      let inside = true;
      for (let k = 0; k < stride; k += 1) {
        if (!near[indexArray[i + k]]) {
          inside = false;
          break;
        }
      }
      if (!inside) continue;
      for (let k = 0; k < stride; k += 1) kept.push(indexArray[i + k]);
    }
    return kept;
  }

  function carveBuffer(buffer, atoms, radius, stride) {
    const geometry = buffer && buffer.geometry;
    if (!geometry || typeof geometry.getAttribute !== "function" || typeof geometry.getIndex !== "function") {
      return null;
    }
    const position = geometry.getAttribute("position");
    const index = geometry.getIndex();
    if (!position || !index || !index.array.length) return null;
    const drawCount = geometry.drawRange && Number.isFinite(geometry.drawRange.count)
      ? geometry.drawRange.count
      : index.array.length;
    const visible = Math.min(index.array.length, Math.max(0, drawCount));
    const kept = carvedIndex(position.array, index.array, atoms, radius, stride, visible);
    // Replace the index outright. Copying into the old buffer leaves the previous
    // contour in the unused tail, and the next carve draws those triangles again.
    const install = (values) => {
      geometry.setIndex(new index.constructor(values, 1));
      if (typeof geometry.clearGroups === "function") geometry.clearGroups();
      geometry.setDrawRange(0, values.length);
    };
    if (!kept.length) {
      install(new index.array.constructor(0));
      return { kept: 0, total: index.array.length };
    }
    if (kept.length === visible && visible === index.array.length) return { kept: kept.length, total: index.array.length };
    install(new index.array.constructor(kept));
    return { kept: kept.length, total: index.array.length };
  }

  function limitDensityToModel(volumeComponent, structureComponent, selection, radius) {
    const world = modelPoints(structureComponent, selection);
    if (!world.length) return false;
    const local = pointsInVolume(world, volumeComponent);
    const seen = new Set();
    let carved = false;
    volumeComponent.reprList.forEach((entry) => {
      const repr = entry.repr || entry;
      const stride = repr.contour ? 2 : 3;
      (repr.bufferList || []).forEach((buffer) => {
        const geometry = buffer.geometry;
        if (!geometry || seen.has(geometry)) return;
        seen.add(geometry);
        let result = carveBuffer(buffer, local, radius, stride);
        if (result && result.kept === 0) result = carveBuffer(buffer, world, radius, stride);
        if (result && result.kept > 0) carved = true;
      });
    });
    return carved;
  }

  function surfaceRepresentations() {
    if (!volumeComponent) return [];
    return volumeComponent.reprList.filter((entry) => (entry.repr || entry).type === "surface");
  }

  function sigmaSetup(scene) {
    let min = Number(scene && scene.sigma_min);
    let max = Number(scene && scene.sigma_max);
    let sigma = Number(scene && scene.sigma);
    const levels = []
      .concat(scene && scene.positive_levels || [])
      .map(Number)
      .filter(Number.isFinite);
    if (!(Number.isFinite(min) && Number.isFinite(max) && min < max)) {
      const peak = levels.length ? Math.max(...levels.map((level) => Math.abs(level))) : 5;
      const hasNegative = ((scene && scene.negative_levels) || []).length > 0 || levels.some((level) => level < 0);
      min = hasNegative ? -peak : Math.min(...levels.filter((level) => level > 0).concat([peak]));
      max = Math.max(peak, ...levels);
      if (!(min < max)) {
        min = -Math.abs(peak || 5);
        max = Math.abs(peak || 5) || 5;
      }
    }
    if (!Number.isFinite(sigma) || sigma < min || sigma > max) {
      const inside = levels.filter((level) => level >= min && level <= max);
      sigma = inside.length
        ? inside.slice().sort((a, b) => Math.abs(a) - Math.abs(b) || Number(a < 0) - Number(b < 0) || a - b)[0]
        : (min + max) / 2;
    }
    return { min, max, sigma: Math.min(max, Math.max(min, sigma)) };
  }

  function colorForSigma(sigma) {
    if (!SigmaColor) return sigma < 0 ? "#ff2a2a" : "#0033ff";
    return SigmaColor.colorForSigma(sigma, sigmaMin, sigmaMax, colorScale);
  }

  function formatSigma(value) {
    return SigmaColor ? SigmaColor.formatSigma(value) : String(value);
  }

  function sigmaScaleGradient() {
    const steps = 16;
    const stops = [];
    for (let i = 0; i <= steps; i += 1) {
      const t = i / steps;
      const sigma = sigmaMin + (sigmaMax - sigmaMin) * t;
      stops.push(`${colorForSigma(sigma)} ${(t * 100).toFixed(2)}%`);
    }
    return `linear-gradient(to top, ${stops.join(", ")})`;
  }

  function sigmaNoteText() {
    if (colorScale === "temperature") {
      return "Temperature scale: cold blue at the low end of the range, hot red at the high end.";
    }
    return "Blue is positive density and red is negative density. Darker shades are farther from zero.";
  }

  function parseSigmaText(text) {
    const trimmed = String(text).trim();
    if (!/^[+-]?(?:\d+(?:\.\d+)?|\.\d+)$/.test(trimmed)) return null;
    const value = Number(trimmed);
    return Number.isFinite(value) ? value : null;
  }

  function sigmaInRange(value) {
    return value >= sigmaMin - 1e-6 && value <= sigmaMax + 1e-6;
  }

  function paintSigmaReadout(source) {
    const color = colorForSigma(surfaceSigmaLevel);
    if (surfaceSigma && source !== "slider") surfaceSigma.value = String(surfaceSigmaLevel);
    if (surfaceSigmaValue && source !== "number" && document.activeElement !== surfaceSigmaValue) {
      surfaceSigmaValue.value = formatSigma(surfaceSigmaLevel);
    }
    if (surfaceSigmaMin) surfaceSigmaMin.textContent = formatSigma(sigmaMin);
    if (surfaceSigmaMax) surfaceSigmaMax.textContent = formatSigma(sigmaMax);
    if (surfaceSigma) {
      surfaceSigma.min = String(sigmaMin);
      surfaceSigma.max = String(sigmaMax);
      surfaceSigma.step = "0.1";
      surfaceSigma.style.accentColor = color;
      const unit = contourIsSigma ? "sigma" : "map units";
      surfaceSigma.setAttribute("aria-valuetext", `${formatSigma(surfaceSigmaLevel)} ${unit}`);
    }
    if (surfaceSigmaValue) {
      surfaceSigmaValue.min = String(sigmaMin);
      surfaceSigmaValue.max = String(sigmaMax);
      surfaceSigmaValue.step = "0.1";
      surfaceSigmaValue.setAttribute("aria-label", contourIsSigma ? "Exact sigma" : "Exact contour level");
    }
    if (surfaceSigmaLabel) surfaceSigmaLabel.textContent = contourIsSigma ? "Sigma" : "Level";
    if (surfaceSigmaUnit) surfaceSigmaUnit.textContent = contourIsSigma ? "σ" : "";
    if (surfaceSigmaSwatch) {
      surfaceSigmaSwatch.style.background = color;
      surfaceSigmaSwatch.title = contourIsSigma
        ? `${formatSigma(surfaceSigmaLevel)} σ`
        : formatSigma(surfaceSigmaLevel);
    }
    if (surfaceSigmaScale) surfaceSigmaScale.style.background = sigmaScaleGradient();
    if (surfaceSigmaNote) surfaceSigmaNote.textContent = sigmaNoteText();
  }

  function syncSigmaControls(enabled) {
    [surfaceSigma, surfaceSigmaValue, colorSigned, colorTemperature].forEach((control) => {
      if (control) control.disabled = !enabled;
    });
    if (colorSigned) colorSigned.setAttribute("aria-pressed", colorScale === "signed" ? "true" : "false");
    if (colorTemperature) colorTemperature.setAttribute("aria-pressed", colorScale === "temperature" ? "true" : "false");
    paintSigmaReadout();
  }

  function mountSigmaSurface() {
    if (!volumeComponent) return null;
    const generation = ++surfaceGeneration;
    const token = viewerToken;
    // Drop the previous mesh completely. Updating a carved surface in place
    // keeps the old contour in the index buffer, so both stay on screen.
    surfaceRepresentations().slice().forEach((entry) => {
      const tasks = entry.repr && entry.repr.tasks;
      if (tasks && tasks.count) tasks.change(-tasks.count);
      volumeComponent.removeRepresentation(entry);
    });
    const created = volumeComponent.addRepresentation("surface", {
      isolevelType: contourIsSigma ? "sigma" : "value",
      contour: surfaceStyle === "mesh",
      opacity: surfaceOpacity,
      useWorker: true,
      isolevel: surfaceSigmaLevel,
      negateIsolevel: false,
      color: colorForSigma(surfaceSigmaLevel),
    });
    return { created, generation, token };
  }

  function whenSurfaceReady(mounted, callback) {
    if (!mounted) return;
    const tasks = mounted.created && mounted.created.repr && mounted.created.repr.tasks;
    const run = () => {
      if (mounted.token !== viewerToken || mounted.generation !== surfaceGeneration) return;
      callback();
    };
    if (!tasks || !tasks.count) {
      run();
      return;
    }
    const onCount = (_delta, count) => {
      if (mounted.token !== viewerToken || mounted.generation !== surfaceGeneration) {
        tasks.signals.countChanged.remove(onCount);
        return;
      }
      if (count === 0) {
        tasks.signals.countChanged.remove(onCount);
        run();
      }
    };
    tasks.signals.countChanged.add(onCount);
  }

  function finishSigmaSurface(mounted) {
    whenSurfaceReady(mounted, () => {
      if (carveRadius > 0 && volumeComponent && structureComponent) {
        limitDensityToModel(volumeComponent, structureComponent, surfaceSelection, carveRadius);
      }
      if (stage) stage.viewer.requestRender();
    });
  }

  function rebuildSigmaSurface() {
    finishSigmaSurface(mountSigmaSurface());
  }

  function quantizeSigma(value) {
    const stepped = Math.round(value * 10) / 10;
    return Math.min(sigmaMax, Math.max(sigmaMin, stepped));
  }

  function previewSigma(value, source) {
    surfaceSigmaLevel = source === "slider"
      ? quantizeSigma(value)
      : Math.min(sigmaMax, Math.max(sigmaMin, value));
    paintSigmaReadout(source);
    window.clearTimeout(sigmaTimer);
    sigmaTimer = window.setTimeout(() => rebuildSigmaSurface(), 150);
  }

  function commitSigma(value, source) {
    surfaceSigmaLevel = source === "slider"
      ? quantizeSigma(value)
      : Math.min(sigmaMax, Math.max(sigmaMin, value));
    paintSigmaReadout(source);
    window.clearTimeout(sigmaTimer);
    rebuildSigmaSurface();
    if (viewerStatus.textContent && viewerStatus.textContent.includes("must be between")) {
      viewerStatus.textContent = sceneNote;
    }
  }

  function rejectSigmaText() {
    if (surfaceSigmaValue) surfaceSigmaValue.value = formatSigma(surfaceSigmaLevel);
    const name = contourIsSigma ? "Sigma" : "Level";
    const unit = contourIsSigma ? "σ" : "map units";
    viewerStatus.textContent = `${name} must be between ${formatSigma(sigmaMin)} and ${formatSigma(sigmaMax)} ${unit}.`;
  }

  function commitSigmaText() {
    const parsed = surfaceSigmaValue ? parseSigmaText(surfaceSigmaValue.value) : null;
    if (parsed == null || !sigmaInRange(parsed)) {
      rejectSigmaText();
      return;
    }
    commitSigma(parsed, "number");
    if (surfaceSigmaValue) surfaceSigmaValue.value = formatSigma(surfaceSigmaLevel);
  }

  function setColorScale(scale) {
    if (colorScale === scale) return;
    colorScale = scale;
    if (colorSigned) colorSigned.setAttribute("aria-pressed", scale === "signed" ? "true" : "false");
    if (colorTemperature) colorTemperature.setAttribute("aria-pressed", scale === "temperature" ? "true" : "false");
    paintSigmaReadout();
    if (!volumeComponent) return;
    surfaceRepresentations().forEach((entry) => entry.setColor(colorForSigma(surfaceSigmaLevel)));
    if (stage) stage.viewer.requestRender();
  }

  function syncSurfaceControls(enabled) {
    if (surfaceSolid) {
      surfaceSolid.disabled = !enabled;
      surfaceSolid.setAttribute("aria-pressed", surfaceStyle === "solid" ? "true" : "false");
    }
    if (surfaceMesh) {
      surfaceMesh.disabled = !enabled;
      surfaceMesh.setAttribute("aria-pressed", surfaceStyle === "mesh" ? "true" : "false");
    }
    if (surfaceTransparency) {
      surfaceTransparency.disabled = !enabled;
      const percent = Math.round((1 - surfaceOpacity) * 100);
      surfaceTransparency.value = String(percent);
      if (surfaceTransparencyValue) surfaceTransparencyValue.textContent = `${percent}%`;
    }
  }

  // Switching contour rebuilds the isosurface, which drops the carve until the
  // new mesh is ready. Opacity updates the existing buffers in place.
  function applySurfaceStyle() {
    finishSigmaSurface(mountSigmaSurface());
  }

  function applySurfaceOpacity(opacity) {
    surfaceOpacity = opacity;
    surfaceRepresentations().forEach((entry) => entry.setParameters({ opacity }));
    if (stage) stage.viewer.requestRender();
  }

  function whenStageIdle(token, callback) {
    const current = stage;
    const run = () => {
      if (token !== viewerToken || stage !== current) return;
      callback();
    };
    if (!current.tasks.count) {
      run();
      return;
    }
    const onCount = (_delta, count) => {
      if (token !== viewerToken || stage !== current) {
        current.tasks.signals.countChanged.remove(onCount);
        return;
      }
      if (count === 0) {
        current.tasks.signals.countChanged.remove(onCount);
        run();
      }
    };
    current.tasks.signals.countChanged.add(onCount);
  }

  function disposeViewer() {
    viewerToken += 1;
    window.clearTimeout(sigmaTimer);
    surfaceGeneration += 1;
    resetEditControls();
    if (stage) stage.dispose();
    stage = null;
    structureComponent = null;
    volumeComponent = null;
    lightComponent = null;
    lightRepresentation = null;
    lightAtomCount = 0;
    lightSiteDescription = "";
    lightSiteName = "";
    showLightSite = false;
    lightColor = "yellow";
    lightTransparency = 40;
    siteResidue = null;
    neighbourResidues = [];
    modelStick = null;
    syncNeighbourControls();
    viewer.replaceChildren();
    if (viewerSave) viewerSave.disabled = true;
    if (viewerPymol) viewerPymol.disabled = true;
    if (viewerOrientation) viewerOrientation.disabled = true;
    syncSurfaceControls(false);
    syncSigmaControls(false);
    syncLightControls(false);
  }

  async function loadViewer(scene) {
    disposeViewer();
    lightSiteName = (scene && scene.light_label) || "";
    if (lightSiteButton) lightSiteButton.textContent = lightSiteName || "Light site";
    if (lightHint) {
      lightHint.textContent = lightSiteName
        ? "Ball and stick of this site, including neighbouring residues, from the reference model. Add bond, Add distance, and Add label work on these atoms. Hide this site to click an atom of the model being reconstructed."
        : lightHintDefault;
    }
    const token = viewerToken;
    viewerStatus.textContent = "Loading 3D coordinates and density map…";
    if (!scene || !window.NGL) {
      viewerStatus.textContent = "The interactive viewer could not be loaded. The saved PNG is available below.";
      return;
    }
    try {
      stage = new NGL.Stage(viewer, { backgroundColor: "white" });
      // A left click otherwise recenters the camera on the pick after this
      // handler runs. Deleting a stick or dash invalidates that pick, and the
      // camera leaves the site.
      stage.mouseControls.remove("clickPick-*");
      siteResidue = scene.site
        ? { chain: String(scene.site.chain || "").trim(), resi: Number(scene.site.resi) }
        : null;
      neighbourResidues = (scene.extra_residues || []).map((residue) => ({
        chain: String(residue.chain || "").trim(),
        resi: Number(residue.resi),
        label: residue.label || `${residue.chain}${residue.resi}`,
        shown: true,
      }));
      const shownSelection = structureSelection(true);
      syncNeighbourControls();
      const [structure, volume] = await Promise.all([
        stage.loadFile(scene.model_url, { ext: scene.model_ext || "pdb", defaultRepresentation: false }),
        stage.loadFile(scene.map_url, { ext: "ccp4", defaultRepresentation: false }),
      ]);
      if (token !== viewerToken) return;
      structureComponent = structure;
      volumeComponent = volume;
      surfaceSelection = shownSelection;
      carveRadius = Number(scene.carve) > 0 ? Number(scene.carve) : 0;
      surfaceStyle = scene.surface === "mesh" ? "mesh" : "solid";
      surfaceOpacity = surfaceStyle === "solid" ? 0.55 : 1;
      const sigma = sigmaSetup(scene);
      sigmaMin = sigma.min;
      sigmaMax = sigma.max;
      surfaceSigmaLevel = sigma.sigma;
      colorScale = "signed";
      contourIsSigma = scene.contour_mode !== "absolute";
      applyScriptBonds(structure, scene);
      modelStick = structure.addRepresentation("ball+stick", {
        sele: shownSelection,
        colorScheme: "element",
        multipleBond: false,
        aspectRatio: 2,
        scale: 2,
      });
      if (scene.light_model_url) {
        try {
          lightComponent = await stage.loadFile(scene.light_model_url, {
            ext: scene.light_model_ext || "pdb",
            defaultRepresentation: false,
          });
        } catch (err) {
          lightComponent = null;
        }
        if (token !== viewerToken) return;
        if (lightComponent) {
          const lightAtoms = lightComponent.structure.getAtomIndices(new NGL.Selection(shownSelection));
          lightAtomCount = lightAtoms ? lightAtoms.length : 0;
          const neighbours = (scene.extra_residues || []).map((residue) => `${residue.chain}/${residue.resi}`);
          lightSiteDescription = `chain ${scene.site.chain} residue ${scene.site.resi}`;
          if (neighbours.length) lightSiteDescription += ` and ${neighbours.join(", ")}`;
          lightBondScene = Array.isArray(scene.light_bonds)
            ? { bonds: scene.light_bonds, site: scene.site }
            : null;
          if (lightBondScene) applyScriptBonds(lightComponent, lightBondScene);
          else copySelectionBonds(structure.structure, lightComponent.structure, shownSelection);
          lightComponent.setVisibility(false);
          lightRepresentation = lightComponent.addRepresentation("ball+stick", {
            sele: shownSelection,
            color: lightColor,
            opacity: lightOpacityValue(),
            depthWrite: lightOpacityValue() >= 0.999,
            multipleBond: "off",
            aspectRatio: 2,
            scale: 2,
            visible: false,
          });
          syncLightControls(true);
        }
      }
      collectScriptDistances(structure, scene);
      stage.signals.clicked.add(onStageClick);
      const mounted = mountSigmaSurface();
      if (token !== viewerToken) return;
      stage.autoView(0);
      stage.handleResize();
      viewerSave.disabled = false;
      if (viewerPymol) viewerPymol.disabled = false;
      if (viewerOrientation) viewerOrientation.disabled = false;
      [editBond, editDistance, editDelete, editLabel, editDistanceLabels].forEach((button) => {
        if (button) button.disabled = false;
      });
      syncSurfaceControls(true);
      syncSigmaControls(true);
      const notes = [];
      if (scene.map_file) notes.push(`Density is the Coot site fragment ${scene.map_file}.`);
      if ((scene.bonds || []).length || (scene.distances || []).length) {
        notes.push("Sticks and dashed distances follow the PyMOL script.");
      }
      if (scene.light_model_url && !lightComponent) {
        notes.push("The light coordinates could not be loaded.");
      } else if (lightComponent) {
        notes.push("Light site draws that file in ball and stick.");
      }
      notes.push("Drag to rotate, scroll to zoom, and right-drag to pan. Save PNG keeps this orientation. PyMOL script writes a .pml file. When Current orientation is on, that script uses this view.");
      sceneNote = notes.join(" ");
      viewerStatus.textContent = sceneNote;
      if (carveRadius > 0) {
        whenSurfaceReady(mounted, () => {
          const carved = limitDensityToModel(volume, structure, shownSelection, carveRadius);
          if (!carved || token !== viewerToken) return;
          const radiusText = Number.isInteger(carveRadius) ? String(carveRadius) : String(carveRadius);
          sceneNote = `${sceneNote} Density is kept within ${radiusText} Å of the model atoms, as in the PyMOL figure.`;
          viewerStatus.textContent = sceneNote;
          stage.viewer.requestRender();
        });
      }
    } catch (err) {
      if (token !== viewerToken) return;
      disposeViewer();
      viewerStatus.textContent = "The interactive scene could not be drawn. The saved PNG is available below.";
    }
  }

  function downloadBlob(blob, filename) {
    const link = document.createElement("a");
    const url = URL.createObjectURL(blob);
    link.href = url;
    link.download = filename;
    document.body.append(link);
    link.click();
    link.remove();
    window.setTimeout(() => URL.revokeObjectURL(url), 2000);
  }

  function orientationIsOn() {
    return viewerOrientation && viewerOrientation.getAttribute("aria-pressed") === "true";
  }

  function captureViewerOrientation() {
    const camera = stage.viewer;
    camera.rotationGroup.updateMatrix();
    const elements = camera.rotationGroup.matrix.elements;
    const pan = camera.translationGroup.position;
    return {
      use_orientation: true,
      rotation: [
        elements[0], elements[4], elements[8],
        elements[1], elements[5], elements[9],
        elements[2], elements[6], elements[10],
      ],
      translation: [pan.x, pan.y, pan.z],
      distance: camera.cameraDistance,
      fov: camera.perspectiveCamera.fov,
      near: camera.camera.near,
      far: camera.camera.far,
      width: Math.round(camera.width),
      height: Math.round(camera.height),
    };
  }

  if (viewerOrientation) {
    viewerOrientation.addEventListener("click", () => {
      const turnOn = !orientationIsOn();
      viewerOrientation.setAttribute("aria-pressed", turnOn ? "true" : "false");
      viewerStatus.textContent = turnOn
        ? "Current orientation is on. PyMOL script will use the rotation, zoom, and pan of the 3D viewer."
        : "Current orientation is off. PyMOL script will keep the framing from the calculation.";
    });
  }

  if (viewerPymol) {
    viewerPymol.addEventListener("click", async () => {
      if (!stage || !jobId) {
        viewerStatus.textContent = "The scene is not ready to save.";
        return;
      }
      const useOrientation = orientationIsOn();
      viewerPymol.disabled = true;
      viewerStatus.textContent = useOrientation
        ? "Writing a PyMOL script from the current view…"
        : "Writing a PyMOL script with the calculated framing…";
      try {
        const payload = useOrientation ? captureViewerOrientation() : { use_orientation: false };
        const response = await fetch(`/api/jobs/${jobId}/pymol-script`, {
          method: "POST",
          headers: { "Content-Type": "application/json" },
          body: JSON.stringify(payload),
        });
        const data = await response.json().catch(() => ({}));
        if (!response.ok) {
          viewerStatus.textContent = data.error || "The PyMOL script could not be written.";
          return;
        }
        downloadBlob(new Blob([data.script || ""], { type: "text/plain" }), data.name || "scene.pml");
        appendFile(data.role || "PyMOL script", data.path);
        viewerStatus.textContent = useOrientation
          ? `Saved ${data.name} using the orientation in the 3D viewer. A copy was also downloaded.`
          : `Saved ${data.name} with the framing from the calculation. A copy was also downloaded.`;
      } catch (err) {
        viewerStatus.textContent = "The PyMOL script could not be written.";
      } finally {
        viewerPymol.disabled = !stage;
      }
    });
  }

  viewerSave.addEventListener("click", async () => {
    if (!stage || !jobId) {
      viewerStatus.textContent = "The scene is not ready to save.";
      return;
    }
    viewerSave.disabled = true;
    viewerStatus.textContent = "Saving the current view…";
    try {
      const blob = await stage.makeImage({ factor: 2, antialias: true, transparent: false });
      const body = new FormData();
      body.append("image", blob, "viewer.png");
      const response = await fetch(`/api/jobs/${jobId}/snapshot`, { method: "POST", body });
      const payload = await response.json().catch(() => ({}));
      if (!response.ok) {
        viewerStatus.textContent = payload.error || "The view could not be saved.";
        return;
      }
      downloadBlob(blob, payload.name || "viewer.png");
      appendFile("Viewer snapshot", payload.path);
      viewerStatus.textContent = `Saved ${payload.name} in the output directory. A copy was also downloaded.`;
    } catch (err) {
      viewerStatus.textContent = "The view could not be saved.";
    } finally {
      viewerSave.disabled = !stage;
    }
  });

  viewerReset.addEventListener("click", () => {
    if (stage) stage.autoView(250);
  });

  function bindEditButton(button, mode) {
    if (!button) return;
    button.addEventListener("click", () => {
      if (!structureComponent) return;
      setEditMode(mode);
    });
  }
  bindEditButton(editBond, "bond");
  bindEditButton(editDistance, "distance");
  bindEditButton(editDelete, "delete");
  bindEditButton(editLabel, "label");
  if (editDistanceLabels) {
    editDistanceLabels.addEventListener("click", () => {
      if (!structureComponent) return;
      showDistanceLabels = !showDistanceLabels;
      editDistanceLabels.setAttribute("aria-pressed", showDistanceLabels ? "true" : "false");
      syncDistanceLabels();
      viewerStatus.textContent = showDistanceLabels
        ? "Distance labels show the length of each dashed line in angstroms."
        : "Distance labels are hidden.";
    });
  }
  if (surfaceSolid) {
    surfaceSolid.addEventListener("click", () => {
      if (!volumeComponent || surfaceStyle === "solid") return;
      surfaceStyle = "solid";
      syncSurfaceControls(true);
      applySurfaceStyle();
      viewerStatus.textContent = "The isosurface is solid.";
    });
  }
  if (surfaceMesh) {
    surfaceMesh.addEventListener("click", () => {
      if (!volumeComponent || surfaceStyle === "mesh") return;
      surfaceStyle = "mesh";
      syncSurfaceControls(true);
      applySurfaceStyle();
      viewerStatus.textContent = "The isosurface is mesh.";
    });
  }
  if (surfaceTransparency) {
    surfaceTransparency.addEventListener("input", () => {
      if (!volumeComponent) return;
      const percent = Number(surfaceTransparency.value);
      if (surfaceTransparencyValue) surfaceTransparencyValue.textContent = `${percent}%`;
      applySurfaceOpacity(1 - percent / 100);
    });
  }
  if (surfaceSigma) {
    surfaceSigma.addEventListener("input", () => {
      previewSigma(Number(surfaceSigma.value), "slider");
    });
    surfaceSigma.addEventListener("change", () => {
      commitSigma(Number(surfaceSigma.value), "slider");
    });
  }
  if (surfaceSigmaValue) {
    surfaceSigmaValue.addEventListener("input", () => {
      const parsed = parseSigmaText(surfaceSigmaValue.value);
      if (parsed == null || !sigmaInRange(parsed)) return;
      previewSigma(parsed, "number");
    });
    surfaceSigmaValue.addEventListener("change", () => {
      commitSigmaText();
    });
  }
  if (colorSigned) {
    colorSigned.addEventListener("click", () => setColorScale("signed"));
  }
  if (colorTemperature) {
    colorTemperature.addEventListener("click", () => setColorScale("temperature"));
  }
  if (lightSiteButton) {
    lightSiteButton.addEventListener("click", () => {
      if (!lightRepresentation) return;
      showLightSite = !showLightSite;
      syncLightControls(true);
      applyLightAppearance();
      viewerStatus.textContent = lightStatusText();
    });
  }
  lightColorButtons.forEach((button) => {
    button.addEventListener("click", () => {
      if (!lightRepresentation) return;
      lightColor = button.dataset.lightColor;
      syncLightControls(true);
      applyLightAppearance();
      viewerStatus.textContent = lightStatusText();
    });
  });
  if (lightTransparencyInput) {
    lightTransparencyInput.addEventListener("input", () => {
      if (!lightRepresentation) return;
      lightTransparency = Number(lightTransparencyInput.value);
      if (lightTransparencyValue) lightTransparencyValue.textContent = `${lightTransparency}%`;
      applyLightAppearance();
    });
  }

  window.addEventListener("resize", () => {
    if (stage) stage.handleResize();
  });
})();
