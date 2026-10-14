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
  const editBond = document.getElementById("edit-bond");
  const editDistance = document.getElementById("edit-distance");
  const editDelete = document.getElementById("edit-delete");
  const editLabel = document.getElementById("edit-label");
  const editDistanceLabels = document.getElementById("edit-distance-labels");
  const surfaceSolid = document.getElementById("surface-solid");
  const surfaceMesh = document.getElementById("surface-mesh");
  const surfaceTransparency = document.getElementById("surface-transparency");
  const surfaceTransparencyValue = document.getElementById("surface-transparency-value");
  const caption = document.getElementById("caption");
  const fileList = document.getElementById("file-list");
  const runButton = document.getElementById("run-btn");
  const scanButton = document.getElementById("scan-btn");
  const scanPanel = document.getElementById("scan-panel");
  const scanBody = document.getElementById("scan-body");
  const scanSummary = document.getElementById("scan-summary");
  let stage = null;
  let jobId = null;
  let viewerToken = 0;
  let resultShown = false;
  let structureComponent = null;
  let volumeComponent = null;
  let surfaceSelection = "";
  let carveRadius = 0;
  let surfaceStyle = "solid";
  let surfaceOpacity = 0.55;
  let distanceComponent = null;
  let markerComponent = null;
  let distancePairs = [];
  let labeledAtoms = [];
  let labelRepresentation = null;
  let distanceLabelRepresentation = null;
  let showDistanceLabels = false;
  let editMode = null;
  let pendingIndex = null;
  let sceneNote = "";

  const NGL_COLORS = {
    forest: "#228b22",
    cyan: "#00e5ff",
    blue: "#0033ff",
    density: "#0033ff",
    red: "#ff2a2a",
    salmon: "#fa8072",
    magenta: "#ff00ff",
    hotpink: "#ff69b4",
    marine: "#1a6ca8",
    deepblue: "#00008b",
  };
  const LABEL_COLOR = "#ff7a18";
  const LABEL_EDGE = "#2a160c";

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
    const labels = { queued: "Queued", running: "Running", succeeded: "Finished", failed: "Failed" };
    runState.textContent = labels[payload.status] || payload.status;
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

  function appendFile(role, path) {
    const wrap = document.createElement("div");
    const dt = document.createElement("dt");
    const dd = document.createElement("dd");
    dt.textContent = role;
    dd.textContent = path;
    wrap.append(dt, dd);
    fileList.append(wrap);
  }

  function showResult(data) {
    figure.src = `${data.figure_url}?t=${Date.now()}`;
    caption.textContent = data.caption;
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

  function describeAtom(index) {
    const atom = structureComponent.structure.getAtomProxy();
    atom.index = index;
    const alt = atom.altloc && atom.altloc.trim() ? ` alt ${atom.altloc.trim()}` : "";
    const chain = (atom.chainname || "").trim();
    return `${atom.atomname.trim()} ${chain}/${atom.resno}${alt}`;
  }

  function bondSlot(left, right) {
    const [first, second] = pairKey(left, right);
    const store = structureComponent.structure.bondStore;
    for (let index = 0; index < store.count; index += 1) {
      if (store.atomIndex1[index] === first && store.atomIndex2[index] === second) return index;
    }
    return -1;
  }

  function distanceSlot(left, right) {
    const [first, second] = pairKey(left, right);
    return distancePairs.findIndex((pair) => pair[0] === first && pair[1] === second);
  }

  function refreshSticks() {
    const structure = structureComponent.structure;
    structure.finalizeBonds();
    // The stick view keeps the bond list it had when it was built. A new bond
    // is stored on the structure, but it stays invisible until that list is rebuilt.
    const views = [];
    if (structureComponent.structureView) views.push(structureComponent.structureView);
    (structureComponent.reprList || []).forEach((entry) => {
      const repr = entry.repr || entry;
      if (repr.structureView) views.push(repr.structureView);
    });
    views.forEach((view) => view.refresh());
    structureComponent.rebuildRepresentations();
  }

  function addBondBetween(left, right) {
    if (bondSlot(left, right) >= 0) return false;
    const structure = structureComponent.structure;
    const first = structure.getAtomProxy();
    const second = structure.getAtomProxy();
    first.index = left;
    second.index = right;
    structure.bondStore.addBond(first, second, 1);
    refreshSticks();
    return true;
  }

  function removeBondBetween(left, right) {
    const slot = bondSlot(left, right);
    if (slot < 0) return false;
    const store = structureComponent.structure.bondStore;
    const last = store.count - 1;
    if (slot !== last) {
      store.atomIndex1[slot] = store.atomIndex1[last];
      store.atomIndex2[slot] = store.atomIndex2[last];
      store.bondOrder[slot] = store.bondOrder[last];
    }
    store.count = last;
    refreshSticks();
    return true;
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

  function redrawDistances() {
    const saved = captureView();
    if (distanceComponent) {
      stage.removeComponent(distanceComponent);
      distanceComponent = null;
    }
    if (distancePairs.length && structureComponent) {
      const structure = structureComponent.structure;
      const atom = structure.getAtomProxy();
      const shape = new NGL.Shape("script-distances");
      distancePairs.forEach(([left, right]) => {
        atom.index = left;
        const start = [atom.x, atom.y, atom.z];
        atom.index = right;
        const end = [atom.x, atom.y, atom.z];
        addDashed(shape, start, end, `dist:${left}:${right}`);
      });
      distanceComponent = stage.addComponentFromObject(shape);
      distanceComponent.addRepresentation("buffer", { disablePicking: true });
    }
    restoreView(saved);
    syncDistanceLabels();
  }

  function syncDistanceLabels() {
    const saved = captureView();
    if (distanceLabelRepresentation && structureComponent) {
      structureComponent.removeRepresentation(distanceLabelRepresentation);
      distanceLabelRepresentation = null;
    }
    if (showDistanceLabels && distancePairs.length && structureComponent) {
      distanceLabelRepresentation = structureComponent.addRepresentation("distance", {
        atomPair: distancePairs.map((pair) => pair.slice()),
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
      const distanceLabels = distanceLabelRepresentation.repr || distanceLabelRepresentation;
      distanceLabels.lineVisible = false;
      if (distanceLabels.distanceBuffer) distanceLabels.distanceBuffer.setVisibility(false);
    }
    restoreView(saved);
  }

  function atomName(index) {
    const atom = structureComponent.structure.getAtomProxy();
    atom.index = index;
    return atom.atomname.trim();
  }

  function redrawAtomLabels() {
    const saved = captureView();
    if (labelRepresentation && structureComponent) {
      structureComponent.removeRepresentation(labelRepresentation);
      labelRepresentation = null;
    }
    if (labeledAtoms.length && structureComponent) {
      const names = {};
      labeledAtoms.forEach((index) => {
        names[index] = atomName(index);
      });
      labelRepresentation = structureComponent.addRepresentation("label", {
        sele: `@${labeledAtoms.join(",")}`,
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
      });
    }
    restoreView(saved);
  }

  function toggleAtomLabel(index) {
    const slot = labeledAtoms.indexOf(index);
    const name = atomName(index);
    if (slot >= 0) {
      labeledAtoms.splice(slot, 1);
      redrawAtomLabels();
      viewerStatus.textContent = `Removed the ${name} label. ${nextPrompt()}`;
      return;
    }
    labeledAtoms.push(index);
    redrawAtomLabels();
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

  function showMarker(index) {
    clearMarker();
    const atom = structureComponent.structure.getAtomProxy();
    atom.index = index;
    const shape = new NGL.Shape("pending-atom");
    const point = [atom.x, atom.y, atom.z];
    const arm = 1.15;
    const color = [0.2, 0.85, 1];
    const name = `pick:${index}`;
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
    pendingIndex = null;
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
    pendingIndex = null;
    distancePairs = [];
    labeledAtoms = [];
    labelRepresentation = null;
    distanceLabelRepresentation = null;
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

  function atomIndexFromProxy(proxy) {
    if (!proxy) return null;
    if (proxy.atom) return proxy.atom.index;
    const named = proxy.cylinder || proxy.sphere;
    if (named && typeof named.name === "string" && named.name.startsWith("pick:")) {
      return Number(named.name.slice(5));
    }
    if ((editMode === "bond" || editMode === "distance" || editMode === "label") && proxy.bond) {
      if (proxy.closestBondAtom) return proxy.closestBondAtom.index;
      return proxy.bond.atomIndex1;
    }
    return null;
  }

  function labeledAtomNear(proxy) {
    if (!proxy || !proxy.canvasPosition || !labeledAtoms.length) return null;
    const controls = stage.viewerControls;
    const atom = structureComponent.structure.getAtomProxy();
    let best = null;
    let bestDistance = 36;
    labeledAtoms.forEach((index) => {
      atom.index = index;
      const point = atom.positionToVector3();
      point.applyMatrix4(structureComponent.matrix);
      const place = controls.getPositionOnCanvas(point);
      const distance = Math.hypot(place.x - proxy.canvasPosition.x, place.y - proxy.canvasPosition.y);
      if (distance < bestDistance) {
        bestDistance = distance;
        best = index;
      }
    });
    return best;
  }

  function finishPair(left, right) {
    pendingIndex = null;
    clearMarker();
    const first = describeAtom(left);
    const second = describeAtom(right);
    if (editMode === "bond") {
      if (addBondBetween(left, right)) viewerStatus.textContent = `Added a bond between ${first} and ${second}. ${nextPrompt()}`;
      else viewerStatus.textContent = `${first} and ${second} already have a bond. ${nextPrompt()}`;
      return;
    }
    if (editMode === "distance") {
      if (addDistanceBetween(left, right)) viewerStatus.textContent = `Added a distance between ${first} and ${second}. ${nextPrompt()}`;
      else viewerStatus.textContent = `Those atoms already have a distance. ${nextPrompt()}`;
      return;
    }
    const removedBond = removeBondBetween(left, right);
    const removedDistance = removeDistanceBetween(left, right);
    if (removedBond && removedDistance) viewerStatus.textContent = `Removed the bond and the distance between ${first} and ${second}. ${nextPrompt()}`;
    else if (removedBond) viewerStatus.textContent = `Removed the bond between ${first} and ${second}. ${nextPrompt()}`;
    else if (removedDistance) viewerStatus.textContent = `Removed the distance between ${first} and ${second}. ${nextPrompt()}`;
    else viewerStatus.textContent = `Those atoms have no bond or distance to remove. ${nextPrompt()}`;
  }

  function onStageClick(proxy) {
    if (!editMode || !structureComponent || !proxy) return;
    if (editMode === "delete" && proxy.cylinder && typeof proxy.cylinder.name === "string") {
      const match = /^dist:(\d+):(\d+)$/.exec(proxy.cylinder.name);
      if (match) {
        const left = Number(match[1]);
        const right = Number(match[2]);
        pendingIndex = null;
        clearMarker();
        removeDistanceBetween(left, right);
        viewerStatus.textContent = `Removed the distance between ${describeAtom(left)} and ${describeAtom(right)}. ${nextPrompt()}`;
        return;
      }
    }
    if (editMode === "delete" && proxy.bond) {
      const left = proxy.bond.atomIndex1;
      const right = proxy.bond.atomIndex2;
      pendingIndex = null;
      clearMarker();
      removeBondBetween(left, right);
      viewerStatus.textContent = `Removed the bond between ${describeAtom(left)} and ${describeAtom(right)}. ${nextPrompt()}`;
      return;
    }
    let picked = proxy;
    if ((editMode === "bond" || editMode === "distance" || editMode === "label") && proxy.cylinder && /^dist:/.test(proxy.cylinder.name || "")) {
      picked = repickWithoutDistances(proxy) || proxy;
    }
    const index = atomIndexFromProxy(picked);
    if (editMode === "label") {
      const labeled = index == null || Number.isNaN(index) ? labeledAtomNear(proxy) : index;
      if (labeled == null || Number.isNaN(labeled)) return;
      toggleAtomLabel(labeled);
      return;
    }
    if (index == null || Number.isNaN(index)) return;
    if (pendingIndex == null) {
      pendingIndex = index;
      showMarker(index);
      viewerStatus.textContent = `${describeAtom(index)} selected. Click the second atom.`;
      return;
    }
    if (pendingIndex === index) {
      viewerStatus.textContent = `${describeAtom(index)} is already selected. Click a different atom.`;
      return;
    }
    finishPair(pendingIndex, index);
  }

  function selectionFor(scene) {
    const residues = [`${scene.site.resi}:${scene.site.chain}`];
    (scene.extra_residues || []).forEach((residue) => residues.push(`${residue.resi}:${residue.chain}`));
    return residues.join(" or ");
  }

  function colorFor(colors, index) {
    const color = colors[Math.min(index, colors.length - 1)];
    return NGL_COLORS[color] || color || "#0033ff";
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
  function carvedIndex(position, indexArray, atoms, radius, stride) {
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
    const kept = [];
    for (let i = 0; i + stride <= indexArray.length; i += stride) {
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
    const kept = carvedIndex(position.array, index.array, atoms, radius, stride);
    if (!kept.length) return { kept: 0, total: index.array.length };
    if (kept.length === index.array.length) return { kept: kept.length, total: index.array.length };
    const filtered = new index.array.constructor(kept);
    index.set(filtered);
    index.updateRange.offset = 0;
    index.updateRange.count = filtered.length;
    index.needsUpdate = true;
    geometry.setDrawRange(0, filtered.length);
    return { kept: filtered.length, total: index.array.length };
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
    const contour = surfaceStyle === "mesh";
    let rebuild = false;
    surfaceRepresentations().forEach((entry) => {
      const current = entry.getParameters ? entry.getParameters() : {};
      if (!!current.contour !== contour) rebuild = true;
      entry.setParameters({ contour, opacity: surfaceOpacity });
    });
    const token = viewerToken;
    if (!rebuild) {
      if (stage) stage.viewer.requestRender();
      return;
    }
    whenStageIdle(token, () => {
      if (carveRadius > 0 && volumeComponent && structureComponent) {
        limitDensityToModel(volumeComponent, structureComponent, surfaceSelection, carveRadius);
      }
      if (stage) stage.viewer.requestRender();
    });
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
    resetEditControls();
    if (stage) stage.dispose();
    stage = null;
    viewer.replaceChildren();
    if (viewerSave) viewerSave.disabled = true;
  }

  async function loadViewer(scene) {
    disposeViewer();
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
      const structureSelection = selectionFor(scene);
      const [structure, volume] = await Promise.all([
        stage.loadFile(scene.model_url, { ext: scene.model_ext || "pdb", defaultRepresentation: false }),
        stage.loadFile(scene.map_url, { ext: "ccp4", defaultRepresentation: false }),
      ]);
      if (token !== viewerToken) return;
      structureComponent = structure;
      volumeComponent = volume;
      surfaceSelection = structureSelection;
      carveRadius = Number(scene.carve) > 0 ? Number(scene.carve) : 0;
      surfaceStyle = scene.surface === "mesh" ? "mesh" : "solid";
      surfaceOpacity = surfaceStyle === "solid" ? 0.55 : 1;
      applyScriptBonds(structure, scene);
      structure.addRepresentation("ball+stick", {
        sele: structureSelection,
        colorScheme: "element",
        multipleBond: false,
        aspectRatio: 2,
        scale: 2,
      });
      collectScriptDistances(structure, scene);
      stage.signals.clicked.add(onStageClick);
      const surfaceOptions = {
        isolevelType: scene.contour_mode === "absolute" ? "value" : "sigma",
        contour: surfaceStyle === "mesh",
        opacity: surfaceOpacity,
        useWorker: true,
      };
      (scene.positive_levels || []).forEach((level, index) => {
        volume.addRepresentation("surface", {
          ...surfaceOptions,
          isolevel: level,
          color: colorFor(scene.positive_colors || [], index),
        });
      });
      (scene.negative_levels || []).forEach((level, index) => {
        volume.addRepresentation("surface", {
          ...surfaceOptions,
          isolevel: level,
          negateIsolevel: true,
          color: colorFor(scene.negative_colors || [], index),
        });
      });
      if (token !== viewerToken) return;
      stage.autoView(0);
      stage.handleResize();
      viewerSave.disabled = false;
      [editBond, editDistance, editDelete, editLabel, editDistanceLabels].forEach((button) => {
        if (button) button.disabled = false;
      });
      syncSurfaceControls(true);
      const notes = [];
      if (scene.map_file) notes.push(`Density is the Coot site fragment ${scene.map_file}.`);
      if ((scene.bonds || []).length || (scene.distances || []).length) {
        notes.push("Sticks and dashed distances follow the PyMOL script.");
      }
      notes.push("Drag to rotate, scroll to zoom, and right-drag to pan. Save PNG keeps this orientation.");
      sceneNote = notes.join(" ");
      viewerStatus.textContent = sceneNote;
      if (carveRadius > 0) {
        whenStageIdle(token, () => {
          const carved = limitDensityToModel(volume, structure, structureSelection, carveRadius);
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
      viewerStatus.textContent = "Isosurfaces are solid.";
    });
  }
  if (surfaceMesh) {
    surfaceMesh.addEventListener("click", () => {
      if (!volumeComponent || surfaceStyle === "mesh") return;
      surfaceStyle = "mesh";
      syncSurfaceControls(true);
      applySurfaceStyle();
      viewerStatus.textContent = "Isosurfaces are mesh.";
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

  window.addEventListener("resize", () => {
    if (stage) stage.handleResize();
  });
})();
