const canvas = document.getElementById("structureCanvas");
const ctx = canvas.getContext("2d");
const introScreen = document.getElementById("introScreen");
const introCanvas = document.getElementById("introCanvas");
const introCtx = introCanvas.getContext("2d");
const enterLab = document.getElementById("enterLab");

const runStatus = document.getElementById("runStatus");
const proteinName = document.getElementById("proteinName");
const sequenceInput = document.getElementById("sequenceInput");
const sequenceMetrics = document.getElementById("sequenceMetrics");
const sequenceMap = document.getElementById("sequenceMap");
const analyzeSequence = document.getElementById("analyzeSequence");
const cleanSequence = document.getElementById("cleanSequence");
const loadExample = document.getElementById("loadExample");
const moleculeType = document.getElementById("moleculeType");
const previewMode = document.getElementById("previewMode");
const presetButtons = document.querySelectorAll(".preset-button");
const ligandCode = document.getElementById("ligandCode");
const chainCopies = document.getElementById("chainCopies");
const phValue = document.getElementById("phValue");
const temperaturePreset = document.getElementById("temperaturePreset");
const mutationPosition = document.getElementById("mutationPosition");
const mutationResidue = document.getElementById("mutationResidue");
const addMutation = document.getElementById("addMutation");
const mutationList = document.getElementById("mutationList");
const runPreview = document.getElementById("runPreview");
const downloadJson = document.getElementById("downloadJson");
const predictionCard = document.getElementById("predictionCard");
const resetView = document.getElementById("resetView");
const toggleSpin = document.getElementById("toggleSpin");
const structureFile = document.getElementById("structureFile");
const viewerReadout = document.getElementById("viewerReadout");
const eventLog = document.getElementById("eventLog");
const clearLog = document.getElementById("clearLog");

const validResidues = new Set("ACDEFGHIKLMNPQRSTVWY".split(""));
const hydrophobic = new Set("AILMFWVY".split(""));
const charged = new Set("DEKRH".split(""));
const polar = new Set("CNQST".split(""));
const special = new Set("GP".split(""));

const residueNames = {
  A: "ALA", C: "CYS", D: "ASP", E: "GLU", F: "PHE",
  G: "GLY", H: "HIS", I: "ILE", K: "LYS", L: "LEU",
  M: "MET", N: "ASN", P: "PRO", Q: "GLN", R: "ARG",
  S: "SER", T: "THR", V: "VAL", W: "TRP", Y: "TYR"
};

const examples = [
  {
    name: "zinc_finger_domain",
    sequence: "MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAG"
  },
  {
    name: "small_helical_binder",
    sequence: "MKTAYIAKQRQISFVKSHFSRQDILDLWQNPFQRLVAFWQFGLK"
  },
  {
    name: "enzyme_like_test_chain",
    sequence: "MGSSHHHHHHSSGLVPRGSHMSEQNNTEMTFQIQRIYTKDISFEAPNAPHVFQKDWQPEVK"
  }
];

const state = {
  residues: [],
  imported: false,
  mutations: [],
  experiment: "single_fold",
  rotX: -0.55,
  rotY: 0.74,
  zoom: 42,
  spin: true,
  dragging: false,
  lastX: 0,
  lastY: 0,
  hovered: null,
  projected: []
};

const experimentDefaults = {
  single_fold: { moleculeType: "protein", previewMode: "fold", ligand: "", copies: 1, ph: 7.4, temp: "room" },
  mutagenesis: { moleculeType: "protein", previewMode: "mutation", ligand: "", copies: 1, ph: 7.4, temp: "room" },
  ligand: { moleculeType: "protein_ligand", previewMode: "binding", ligand: "ATP", copies: 1, ph: 7.4, temp: "physio" },
  complex: { moleculeType: "complex", previewMode: "fold", ligand: "", copies: 2, ph: 7.4, temp: "physio" },
  stability: { moleculeType: "protein", previewMode: "mutation", ligand: "", copies: 1, ph: 6.8, temp: "heat" },
  disorder: { moleculeType: "protein", previewMode: "fold", ligand: "", copies: 1, ph: 7.4, temp: "room" }
};

function clean(seq) {
  return seq.toUpperCase().replace(/[^A-Z]/g, "");
}

function classify(residue) {
  if (hydrophobic.has(residue)) return "hydrophobic";
  if (charged.has(residue)) return "charged";
  if (polar.has(residue)) return "polar";
  return "special";
}

function formatPercent(value) {
  return `${Math.round(value)}%`;
}

function molecularWeight(sequence) {
  const weights = {
    A: 89.09, C: 121.16, D: 133.1, E: 147.13, F: 165.19,
    G: 75.07, H: 155.16, I: 131.17, K: 146.19, L: 131.17,
    M: 149.21, N: 132.12, P: 115.13, Q: 146.15, R: 174.2,
    S: 105.09, T: 119.12, V: 117.15, W: 204.23, Y: 181.19
  };
  if (!sequence.length) return 0;
  const total = sequence.split("").reduce((sum, residue) => sum + (weights[residue] || 0), 18.02);
  return Math.max(0, total - (sequence.length - 1) * 18.02);
}

function estimateIsoelectricPoint(sequence) {
  if (!sequence.length) return 0;
  const positives = (sequence.match(/[KRH]/g) || []).length;
  const negatives = (sequence.match(/[DE]/g) || []).length;
  return Math.max(3.8, Math.min(11.2, 7 + (positives - negatives) * 0.18));
}

function analyze(seq) {
  const invalid = seq.split("").filter((residue) => !validResidues.has(residue));
  const usable = seq.split("").filter((residue) => validResidues.has(residue)).join("");
  const counts = { hydrophobic: 0, charged: 0, polar: 0, special: 0 };
  usable.split("").forEach((residue) => {
    counts[classify(residue)] += 1;
  });
  const disorderRisk = usable.length
    ? Math.min(92, Math.max(8, counts.special * 2.4 + counts.charged * 0.9 - counts.hydrophobic * 0.35))
    : 0;
  const pocketSignal = usable.length
    ? Math.min(96, Math.max(12, counts.hydrophobic * 1.5 + (usable.match(/[FWY]/g) || []).length * 2.5))
    : 0;

  return {
    sequence: usable,
    invalid,
    counts,
    length: usable.length,
    mass: molecularWeight(usable),
    pI: estimateIsoelectricPoint(usable),
    hydrophobicity: usable.length ? counts.hydrophobic / usable.length * 100 : 0,
    chargeDensity: usable.length ? counts.charged / usable.length * 100 : 0,
    disorderRisk,
    pocketSignal
  };
}

function buildSyntheticStructure(sequence) {
  const residues = [];
  const len = Math.max(sequence.length, 1);
  for (let i = 0; i < sequence.length; i += 1) {
    const t = i / Math.max(1, len - 1);
    const turn = i * 0.72;
    const helixBias = Math.sin(i / 8) * 1.2;
    const radius = 2.2 + (classify(sequence[i]) === "hydrophobic" ? 0.32 : 0);
    residues.push({
      id: i + 1,
      code: sequence[i],
      name: residueNames[sequence[i]],
      type: classify(sequence[i]),
      x: Math.cos(turn) * radius + (t - 0.5) * 10.5,
      y: Math.sin(turn) * radius + helixBias,
      z: Math.sin(i / 5) * 2.2 + Math.cos(turn) * 0.7
    });
  }
  return residues;
}

function normalizeImported(points) {
  if (!points.length) return [];
  const center = points.reduce((acc, point) => {
    acc.x += point.x;
    acc.y += point.y;
    acc.z += point.z;
    return acc;
  }, { x: 0, y: 0, z: 0 });
  center.x /= points.length;
  center.y /= points.length;
  center.z /= points.length;
  const maxDistance = Math.max(...points.map((point) => Math.hypot(point.x - center.x, point.y - center.y, point.z - center.z)), 1);
  return points.map((point, index) => ({
    ...point,
    id: point.id || index + 1,
    x: (point.x - center.x) / maxDistance * 7,
    y: (point.y - center.y) / maxDistance * 7,
    z: (point.z - center.z) / maxDistance * 7,
    type: classify(point.code || "A")
  }));
}

function parsePdb(text) {
  const residues = [];
  const seen = new Set();
  text.split(/\r?\n/).forEach((line) => {
    if (!line.startsWith("ATOM") && !line.startsWith("HETATM")) return;
    const atom = line.slice(12, 16).trim();
    if (atom !== "CA" && atom !== "P") return;
    const chain = line.slice(21, 22).trim() || "A";
    const seqId = Number.parseInt(line.slice(22, 26), 10);
    const key = `${chain}:${seqId}`;
    if (seen.has(key)) return;
    seen.add(key);
    const name = line.slice(17, 20).trim() || "UNK";
    residues.push({
      id: residues.length + 1,
      code: Object.keys(residueNames).find((code) => residueNames[code] === name) || "A",
      name,
      chain,
      x: Number.parseFloat(line.slice(30, 38)),
      y: Number.parseFloat(line.slice(38, 46)),
      z: Number.parseFloat(line.slice(46, 54))
    });
  });
  return normalizeImported(residues.filter((point) => Number.isFinite(point.x + point.y + point.z)));
}

function parseCif(text) {
  const residues = [];
  text.split(/\r?\n/).forEach((line) => {
    const trimmed = line.trim();
    if (!trimmed.startsWith("ATOM ") && !trimmed.startsWith("HETATM ")) return;
    const parts = trimmed.split(/\s+/);
    const atomIndex = parts.findIndex((part) => part === "CA" || part === "P");
    if (atomIndex === -1) return;
    const floats = parts.map((part) => Number.parseFloat(part));
    const numeric = floats.map((value, index) => ({ value, index })).filter((item) => Number.isFinite(item.value));
    if (numeric.length < 3) return;
    const coords = numeric.slice(-5, -2).length === 3 ? numeric.slice(-5, -2) : numeric.slice(-3);
    const name = parts.find((part) => Object.values(residueNames).includes(part)) || "UNK";
    residues.push({
      id: residues.length + 1,
      code: Object.keys(residueNames).find((code) => residueNames[code] === name) || "A",
      name,
      x: coords[0].value,
      y: coords[1].value,
      z: coords[2].value
    });
  });
  return normalizeImported(residues);
}

function getSequence() {
  return analyze(clean(sequenceInput.value)).sequence;
}

function renderMetrics(stats) {
  const status = stats.invalid.length
    ? `${stats.invalid.length} unsupported symbols removed`
    : "Sequence valid";
  sequenceMetrics.innerHTML = `
    <div class="metric"><span>Length</span><strong>${stats.length}</strong></div>
    <div class="metric"><span>Mass</span><strong>${(stats.mass / 1000).toFixed(1)} kDa</strong></div>
    <div class="metric"><span>Hydrophobic</span><strong>${formatPercent(stats.hydrophobicity)}</strong></div>
    <div class="metric"><span>pI estimate</span><strong>${stats.pI.toFixed(1)}</strong></div>
    <div class="metric"><span>Disorder risk</span><strong>${formatPercent(stats.disorderRisk)}</strong></div>
    <div class="metric"><span>Input status</span><strong>${status}</strong></div>
  `;
}

function renderSequenceMap(sequence) {
  sequenceMap.innerHTML = "";
  if (!sequence.length) {
    sequenceMap.textContent = "Enter a protein sequence to inspect residue classes.";
    return;
  }
  sequence.split("").forEach((residue, index) => {
    const node = document.createElement("button");
    node.type = "button";
    node.className = `residue-cell ${classify(residue)}`;
    node.textContent = residue;
    node.title = `${index + 1}: ${residueNames[residue]} (${classify(residue)})`;
    node.addEventListener("click", () => {
      mutationPosition.value = index + 1;
      state.hovered = state.residues[index] || null;
      viewerReadout.textContent = `${index + 1} ${residueNames[residue]} selected for mutation planning.`;
    });
    sequenceMap.appendChild(node);
  });
}

function renderMutations() {
  mutationList.innerHTML = "";
  if (!state.mutations.length) {
    const item = document.createElement("li");
    item.textContent = "No variants queued.";
    mutationList.appendChild(item);
    return;
  }
  state.mutations.forEach((mutation, index) => {
    const item = document.createElement("li");
    item.innerHTML = `<strong>${mutation.from}${mutation.position}${mutation.to}</strong> ${mutation.effect}`;
    const remove = document.createElement("button");
    remove.type = "button";
    remove.textContent = "Remove";
    remove.addEventListener("click", () => {
      state.mutations.splice(index, 1);
      renderMutations();
      runLocalPreview();
    });
    item.append(" ");
    item.appendChild(remove);
    mutationList.appendChild(item);
  });
}

function mutationEffect(from, to) {
  if (from === to) return "same residue; neutral control";
  if (special.has(to)) return "may increase flexibility or disorder";
  if (hydrophobic.has(from) && charged.has(to)) return "may expose a charge in the core";
  if (charged.has(from) && hydrophobic.has(to)) return "may stabilize a buried pocket";
  if (polar.has(to)) return "adds hydrogen-bonding potential";
  return "changes local packing";
}

function addMutationToQueue() {
  const sequence = getSequence();
  const position = Number.parseInt(mutationPosition.value, 10);
  if (!sequence.length || position < 1 || position > sequence.length) {
    logEvent("Mutation rejected.", "Choose a residue position inside the current sequence.");
    return;
  }
  const from = sequence[position - 1];
  const to = mutationResidue.value;
  state.mutations.push({ position, from, to, effect: mutationEffect(from, to) });
  renderMutations();
  logEvent("Variant queued.", `${from}${position}${to} added to the experiment set.`);
}

function selectExperiment(experiment) {
  state.experiment = experiment;
  const defaults = experimentDefaults[experiment];
  moleculeType.value = defaults.moleculeType;
  previewMode.value = defaults.previewMode;
  ligandCode.value = defaults.ligand;
  chainCopies.value = defaults.copies;
  phValue.value = defaults.ph;
  temperaturePreset.value = defaults.temp;
  presetButtons.forEach((button) => {
    button.classList.toggle("active", button.dataset.experiment === experiment);
  });
  if (experiment === "mutagenesis" && !state.mutations.length) {
    mutationResidue.value = "A";
  }
  if (experiment === "stability" && !state.mutations.length) {
    mutationResidue.value = "L";
  }
  logEvent("Experiment selected.", buttonLabel(experiment));
  runLocalPreview();
}

function buttonLabel(experiment) {
  const match = [...presetButtons].find((button) => button.dataset.experiment === experiment);
  return match ? match.querySelector("strong").textContent : experiment;
}

function buildAf3Json() {
  const sequence = getSequence();
  const name = (proteinName.value || "foldbench_job").trim().replace(/[^a-zA-Z0-9_-]/g, "_");
  const mutated = sequence.split("");
  state.mutations.forEach((mutation) => {
    mutated[mutation.position - 1] = mutation.to;
  });
  const copies = Math.max(1, Math.min(4, Number.parseInt(chainCopies.value, 10) || 1));
  const chain = {
    protein: {
      id: "A",
      sequence: mutated.join("")
    }
  };
  const payload = {
    name,
    modelSeeds: [1],
    sequences: [chain],
    dialect: "alphafold3",
    version: 2,
    metadata: {
      source: "FoldBench local workbench",
      experiment: state.experiment,
      previewMode: previewMode.value,
      moleculeType: moleculeType.value,
      ligandCode: ligandCode.value || null,
      chainCopies: copies,
      conditions: {
        pH: Number.parseFloat(phValue.value) || 7.4,
        temperature: temperaturePreset.value
      },
      variants: state.mutations.map(({ position, from, to }) => ({ position, from, to }))
    }
  };
  if (moleculeType.value === "protein_ligand" && ligandCode.value) {
    payload.sequences.push({
      ligand: {
        id: "L",
        ccdCodes: [ligandCode.value]
      }
    });
  }
  if (moleculeType.value === "complex" || copies > 1) {
    const ids = ["B", "C", "D"].slice(0, copies - 1);
    ids.forEach((id) => {
      payload.sequences.push({
        protein: {
          id,
          sequence: mutated.join("")
        }
      });
    });
  }
  return payload;
}

function scorePreview(stats) {
  let confidence = 58 + stats.hydrophobicity * 0.24 - stats.disorderRisk * 0.18 + stats.pocketSignal * 0.1;
  let ptm = 52 + stats.length * 0.14 - stats.chargeDensity * 0.18;
  let pae = 18 - stats.hydrophobicity * 0.05 + stats.disorderRisk * 0.06;
  const pH = Number.parseFloat(phValue.value) || 7.4;
  const copies = Math.max(1, Math.min(4, Number.parseInt(chainCopies.value, 10) || 1));

  state.mutations.forEach((mutation) => {
    if (charged.has(mutation.from) && hydrophobic.has(mutation.to)) confidence += 3.5;
    if (special.has(mutation.to)) confidence -= 4;
    if (mutation.from === mutation.to) confidence -= 1;
    ptm += hydrophobic.has(mutation.to) ? 1.4 : 0.2;
    pae += special.has(mutation.to) ? 1.5 : -0.3;
  });

  if (previewMode.value === "binding") confidence += stats.pocketSignal * 0.05;
  if (previewMode.value === "mutation") confidence += state.mutations.length * 1.5;
  if (ligandCode.value) {
    confidence += stats.pocketSignal > 45 ? 4 : -2;
    ptm += 2;
  }
  if (moleculeType.value === "complex" || copies > 1) {
    confidence -= 3 + copies;
    ptm += copies * 1.5;
    pae += 1.4 + copies * 0.6;
  }
  if (Math.abs(pH - stats.pI) > 2.8) {
    confidence -= 5;
    pae += 1.7;
  }
  if (temperaturePreset.value === "cold") confidence += 1;
  if (temperaturePreset.value === "physio") ptm += 1.5;
  if (temperaturePreset.value === "heat") {
    confidence -= 4;
    pae += 1.5;
  }
  if (state.experiment === "disorder") {
    confidence -= stats.disorderRisk * 0.12;
    pae += stats.disorderRisk * 0.03;
  }
  if (state.experiment === "stability") {
    confidence += state.mutations.length ? 2 : -2;
  }

  return {
    confidence: Math.max(28, Math.min(96, Math.round(confidence))),
    ptm: Math.max(20, Math.min(92, Math.round(ptm))),
    pae: Math.max(3, Math.min(28, Number(pae.toFixed(1))))
  };
}

function runLocalPreview() {
  const stats = analyze(clean(sequenceInput.value));
  const sequence = stats.sequence;
  if (!sequence.length) {
    predictionCard.innerHTML = "<strong>No sequence yet.</strong><p>Add a protein sequence to run the local preview.</p>";
    return;
  }
  renderMetrics(stats);
  renderSequenceMap(sequence);
  if (!state.imported) state.residues = buildSyntheticStructure(sequence);
  const scores = scorePreview(stats);
  const verdict = scores.confidence >= 82 ? "High-confidence candidate" : scores.confidence >= 66 ? "Promising candidate" : "Needs more design";
  const protocol = buttonLabel(state.experiment);
  const conditionSummary = `${protocol}; ${ligandCode.value || "no ligand"}; ${chainCopies.value} chain copy/copies; pH ${phValue.value}; ${temperaturePreset.options[temperaturePreset.selectedIndex].text}.`;
  predictionCard.innerHTML = `
    <strong>${verdict}</strong>
    <p>${conditionSummary} This is a fast local estimate for triage. Use the exported JSON with an AlphaFold 3 runner, then import the returned CIF/PDB here for inspection.</p>
    <div class="result-score">
      <div class="score-box"><span>pLDDT-like</span><strong>${scores.confidence}</strong></div>
      <div class="score-box"><span>pTM-like</span><strong>${scores.ptm}</strong></div>
      <div class="score-box"><span>PAE-like</span><strong>${scores.pae} Å</strong></div>
    </div>
  `;
  runStatus.textContent = `Preview complete: ${scores.confidence}%`;
  viewerReadout.textContent = `${sequence.length} residues rendered. Import an AF3 result to replace this synthetic preview.`;
  logEvent("Preview run.", `${verdict} with ${scores.confidence} pLDDT-like confidence.`);
}

function downloadAf3Json() {
  const payload = buildAf3Json();
  const blob = new Blob([JSON.stringify(payload, null, 2)], { type: "application/json" });
  const link = document.createElement("a");
  link.href = URL.createObjectURL(blob);
  link.download = `${payload.name || "foldbench_job"}.af3.json`;
  link.click();
  URL.revokeObjectURL(link.href);
  logEvent("AF3 JSON exported.", "Input package downloaded for a backend or local AlphaFold 3 run.");
}

function logEvent(title, detail) {
  const item = document.createElement("li");
  item.innerHTML = `<strong>${title}</strong> ${detail}`;
  eventLog.prepend(item);
  while (eventLog.children.length > 8) eventLog.lastElementChild.remove();
}

function resizeCanvas() {
  const rect = canvas.getBoundingClientRect();
  const ratio = window.devicePixelRatio || 1;
  canvas.width = Math.max(640, Math.floor(rect.width * ratio));
  canvas.height = Math.max(420, Math.floor(rect.height * ratio));
  ctx.setTransform(ratio, 0, 0, ratio, 0, 0);
}

function resizeIntroCanvas() {
  const rect = introCanvas.getBoundingClientRect();
  const ratio = window.devicePixelRatio || 1;
  introCanvas.width = Math.max(640, Math.floor(rect.width * ratio));
  introCanvas.height = Math.max(420, Math.floor(rect.height * ratio));
  introCtx.setTransform(ratio, 0, 0, ratio, 0, 0);
}

function rotate(point) {
  const sinX = Math.sin(state.rotX);
  const cosX = Math.cos(state.rotX);
  const sinY = Math.sin(state.rotY);
  const cosY = Math.cos(state.rotY);
  const y1 = point.y * cosX - point.z * sinX;
  const z1 = point.y * sinX + point.z * cosX;
  const x2 = point.x * cosY + z1 * sinY;
  const z2 = -point.x * sinY + z1 * cosY;
  return { x: x2, y: y1, z: z2 };
}

function project(point) {
  const rect = canvas.getBoundingClientRect();
  const rotated = rotate(point);
  const depth = 9 / (9 + rotated.z);
  return {
    residue: point,
    x: rect.width / 2 + rotated.x * state.zoom * depth,
    y: rect.height / 2 + rotated.y * state.zoom * depth,
    z: rotated.z,
    size: Math.max(4.5, 8.5 * depth)
  };
}

function drawBackbone(points) {
  if (points.length < 2) return;
  ctx.save();
  ctx.lineCap = "round";
  ctx.lineJoin = "round";
  ctx.strokeStyle = "rgba(0,0,0,0.48)";
  ctx.lineWidth = 15;
  ctx.beginPath();
  ctx.moveTo(points[0].x, points[0].y);
  for (let i = 1; i < points.length - 1; i += 1) {
    const midX = (points[i].x + points[i + 1].x) / 2;
    const midY = (points[i].y + points[i + 1].y) / 2;
    ctx.quadraticCurveTo(points[i].x, points[i].y, midX, midY);
  }
  ctx.lineTo(points[points.length - 1].x, points[points.length - 1].y);
  ctx.stroke();
  ctx.strokeStyle = "rgba(47, 130, 196, 0.94)";
  ctx.lineWidth = 8;
  ctx.stroke();
  ctx.strokeStyle = "rgba(247, 240, 228, 0.75)";
  ctx.lineWidth = 2;
  ctx.stroke();
  ctx.restore();
}

function drawResidue(item) {
  const residue = item.residue;
  const isMutation = state.mutations.some((mutation) => mutation.position === residue.id);
  const isHover = state.hovered && state.hovered.id === residue.id;
  const colors = {
    hydrophobic: "#73bd85",
    charged: "#e87970",
    polar: "#79b8e8",
    special: "#eac25d"
  };
  const radius = item.size + (isMutation ? 4 : 0) + (isHover ? 3 : 0);
  ctx.beginPath();
  ctx.arc(item.x, item.y, radius, 0, Math.PI * 2);
  ctx.fillStyle = colors[residue.type] || "#edf4f2";
  ctx.fill();
  ctx.lineWidth = isMutation || isHover ? 3 : 1.2;
  ctx.strokeStyle = isMutation ? "#ffffff" : "rgba(255,255,255,0.55)";
  ctx.stroke();
  if (isMutation || isHover) {
    ctx.fillStyle = "#edf4f2";
    ctx.font = "700 12px Inter, sans-serif";
    ctx.textAlign = "center";
    ctx.fillText(`${residue.code || residue.name}${residue.id}`, item.x, item.y - radius - 9);
  }
}

function render() {
  const rect = canvas.getBoundingClientRect();
  ctx.clearRect(0, 0, rect.width, rect.height);
  if (state.spin && !state.dragging) state.rotY += 0.0022;
  const projected = state.residues.map(project).sort((a, b) => a.z - b.z);
  state.projected = projected;
  drawBackbone([...projected].sort((a, b) => a.residue.id - b.residue.id));
  projected.forEach(drawResidue);
  requestAnimationFrame(render);
}

function drawIntro() {
  const rect = introCanvas.getBoundingClientRect();
  const width = rect.width;
  const height = rect.height;
  const time = performance.now() * 0.00045;
  introCtx.clearRect(0, 0, width, height);

  const centerY = height * 0.54;
  introCtx.lineCap = "round";

  for (let row = 0; row < 5; row += 1) {
    const y = centerY + (row - 2) * 34;
    const start = width * 0.16 + Math.sin(time * 1.8 + row) * 24;
    const end = width * 0.84 + Math.cos(time * 1.5 + row) * 18;
    introCtx.beginPath();
    introCtx.moveTo(start, y);
    introCtx.lineTo(end, y);
    introCtx.strokeStyle = row === 2 ? "rgba(121, 184, 232, 0.36)" : "rgba(240, 236, 226, 0.09)";
    introCtx.lineWidth = row === 2 ? 3 : 1.5;
    introCtx.stroke();
  }

  introCtx.beginPath();
  for (let i = 0; i <= 120; i += 1) {
    const x = width * 0.18 + (width * 0.64 * i) / 120;
    const y = centerY + Math.sin(i * 0.12 + time * 3) * 18 + Math.cos(i * 0.05) * 26;
    if (i === 0) introCtx.moveTo(x, y);
    else introCtx.lineTo(x, y);
  }
  introCtx.strokeStyle = "rgba(232, 162, 59, 0.72)";
  introCtx.lineWidth = 4;
  introCtx.stroke();

  ["#79b8e8", "#e8a23b", "#73bd85"].forEach((color, index) => {
    introCtx.beginPath();
    introCtx.arc(width * (0.28 + index * 0.22), centerY + Math.sin(time * 2 + index) * 18, 5, 0, Math.PI * 2);
    introCtx.fillStyle = color;
    introCtx.fill();
  });

  requestAnimationFrame(drawIntro);
}

function canvasPoint(event) {
  const rect = canvas.getBoundingClientRect();
  return { x: event.clientX - rect.left, y: event.clientY - rect.top };
}

function findResidue(point) {
  return [...state.projected]
    .reverse()
    .find((item) => Math.hypot(item.x - point.x, item.y - point.y) < item.size + 7)?.residue || null;
}

canvas.addEventListener("pointerdown", (event) => {
  canvas.setPointerCapture(event.pointerId);
  state.dragging = true;
  state.lastX = event.clientX;
  state.lastY = event.clientY;
  const hit = findResidue(canvasPoint(event));
  if (hit) {
    state.hovered = hit;
    mutationPosition.value = hit.id;
    viewerReadout.textContent = `${hit.id} ${hit.name || residueNames[hit.code] || hit.code} (${hit.type}) selected.`;
  }
});

canvas.addEventListener("pointermove", (event) => {
  const hit = findResidue(canvasPoint(event));
  if (hit && !state.dragging) state.hovered = hit;
  if (!state.dragging) return;
  state.rotY += (event.clientX - state.lastX) * 0.006;
  state.rotX += (event.clientY - state.lastY) * 0.006;
  state.lastX = event.clientX;
  state.lastY = event.clientY;
});

canvas.addEventListener("pointerup", () => {
  state.dragging = false;
});

canvas.addEventListener("wheel", (event) => {
  event.preventDefault();
  state.zoom = Math.max(22, Math.min(84, state.zoom - event.deltaY * 0.035));
}, { passive: false });

structureFile.addEventListener("change", async (event) => {
  const [file] = event.target.files;
  if (!file) return;
  const text = await file.text();
  const parsed = file.name.toLowerCase().endsWith(".pdb") ? parsePdb(text) : parseCif(text);
  if (!parsed.length) {
    logEvent("Import failed.", "No CA/P atoms were found in that structure file.");
    return;
  }
  state.imported = true;
  state.residues = parsed;
  viewerReadout.textContent = `${parsed.length} residues imported from ${file.name}.`;
  runStatus.textContent = "Imported structure loaded";
  logEvent("Structure imported.", `${parsed.length} residues loaded into the viewer.`);
});

analyzeSequence.addEventListener("click", () => {
  state.imported = false;
  runLocalPreview();
});

cleanSequence.addEventListener("click", () => {
  sequenceInput.value = getSequence();
  state.imported = false;
  runLocalPreview();
});

loadExample.addEventListener("click", () => {
  const sample = examples[Math.floor(Math.random() * examples.length)];
  proteinName.value = sample.name;
  sequenceInput.value = sample.sequence;
  state.mutations = [];
  state.imported = false;
  renderMutations();
  runLocalPreview();
});

presetButtons.forEach((button) => {
  button.addEventListener("click", () => selectExperiment(button.dataset.experiment));
});
[
  moleculeType,
  previewMode,
  ligandCode,
  chainCopies,
  phValue,
  temperaturePreset
].forEach((control) => {
  control.addEventListener("change", runLocalPreview);
});
addMutation.addEventListener("click", addMutationToQueue);
runPreview.addEventListener("click", runLocalPreview);
downloadJson.addEventListener("click", downloadAf3Json);
resetView.addEventListener("click", () => {
  state.rotX = -0.55;
  state.rotY = 0.74;
  state.zoom = 42;
});
toggleSpin.addEventListener("click", () => {
  state.spin = !state.spin;
  toggleSpin.textContent = state.spin ? "⟳" : "Ⅱ";
});
clearLog.addEventListener("click", () => {
  eventLog.innerHTML = "";
});
window.addEventListener("resize", resizeCanvas);
window.addEventListener("resize", resizeIntroCanvas);
enterLab.addEventListener("click", () => {
  introScreen.classList.add("hidden");
  document.body.classList.remove("intro-active");
  resizeCanvas();
});

document.body.classList.add("intro-active");
resizeCanvas();
resizeIntroCanvas();
renderMutations();
runLocalPreview();
logEvent("Workbench online.", "Prepare AF3 inputs, test variants, export JSON, and import returned structures.");
drawIntro();
render();
