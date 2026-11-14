(() => {
  "use strict";

  const canvas = document.getElementById("space");
  const gpuCanvas = document.getElementById("gpuSpace");
  const ctx = canvas.getContext("2d", { alpha: true });
  const els = {
    count: document.getElementById("count"),
    speed: document.getElementById("speed"),
    zoom: document.getElementById("zoom"),
    accuracy: document.getElementById("accuracy"),
    boundary: document.getElementById("boundary"),
    renderer: document.getElementById("renderer"),
    rendererStatus: document.getElementById("rendererStatus"),
    tree: document.getElementById("tree"),
    density: document.getElementById("density"),
    massMin: document.getElementById("massMin"),
    massMax: document.getElementById("massMax"),
    densityStatus: document.getElementById("densityStatus"),
    expansion: document.getElementById("expansion"),
    expansionStatus: document.getElementById("expansionStatus"),
    gravity: document.getElementById("gravity"),
    gravityValue: document.getElementById("gravityValue"),
    gravityPower: document.getElementById("gravityPower"),
    gravityPowerValue: document.getElementById("gravityPowerValue"),
    softening: document.getElementById("softening"),
    softeningValue: document.getElementById("softeningValue"),
    gravityStatus: document.getElementById("gravityStatus"),
    countValue: document.getElementById("countValue"),
    speedValue: document.getElementById("speedValue"),
    zoomValue: document.getElementById("zoomValue"),
    accuracyValue: document.getElementById("accuracyValue"),
    bodyStat: document.getElementById("bodyStat"),
    fpsStat: document.getElementById("fpsStat"),
    cellStat: document.getElementById("cellStat"),
    timingStat: document.getElementById("timingStat"),
    nlogStat: document.getElementById("nlogStat"),
    gainStat: document.getElementById("gainStat"),
    pause: document.getElementById("pause"),
    reset: document.getElementById("reset")
  };

  const cpuFMM = createOrbitFMMField();
  const CPU_MAX_BODIES = 262144;
  const MAX_DPR = 2;
  const BOX_MIN = -1;
  const BOX_MAX = 1;
  const BOX_SIZE = BOX_MAX - BOX_MIN;
  const palette = [
    [98, 217, 255],
    [157, 114, 255],
    [215, 255, 82],
    [255, 157, 85]
  ];
  const particleCanvas = document.createElement("canvas");
  const particleCtx = particleCanvas.getContext("2d");

  let width = 0;
  let height = 0;
  let dpr = 1;
  let bodies = [];
  let running = true;
  let lastTime = performance.now();
  let fpsSmooth = 60;
  let treeSnapshot = null;
  let treeCache = null;
  let operationCount = 1;
  let impulse = null;
  let particleImage = null;
  let fmmMsSmooth = 0;
  let renderMsSmooth = 0;
  let gpuRenderer = null;
  let gpuInitPromise = null;
  let gpuParticleData = null;
  let gpuComputeEngine = null;
  let gpuSeedConfig = null;
  let simulationCount = 512;
  let rendererRequest = 0;
  let boundaryMode = els.boundary.value;
  let torusActive = boundaryMode === "torus";
  let cameraZoom = 1;
  let cosmicTime = 0;
  let expansionFunction = () => 1;
  let expansionValid = true;
  let lastExpansionStatusTick = -1;
  let gravityStrength = 0.48;
  let gravityExponent = 2;
  let gravitySoftening2 = 0.042 ** 2;

  const randn = () => {
    const u = Math.max(1e-9, Math.random());
    const v = Math.random();
    return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
  };

  function resize() {
    const rect = canvas.getBoundingClientRect();
    dpr = Math.min(MAX_DPR, window.devicePixelRatio || 1);
    width = Math.max(320, rect.width);
    height = Math.max(320, rect.height);
    canvas.width = Math.round(width * dpr);
    canvas.height = Math.round(height * dpr);
    gpuCanvas.width = Math.round(width * dpr);
    gpuCanvas.height = Math.round(height * dpr);
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    particleCanvas.width = Math.round(width);
    particleCanvas.height = Math.round(height);
    particleImage = particleCtx.createImageData(
      particleCanvas.width,
      particleCanvas.height
    );
    if (gpuRenderer) gpuRenderer.configure();
    if (gpuComputeEngine) {
      gpuComputeEngine.configure(width, height, dpr);
    }
  }

  function selectedBodyCount() {
    const exponent = Number(els.count.value);
    return exponent >= 25.574
      ? 50000000
      : Math.round(2 ** exponent);
  }

  // Construit numériquement F(m), puis échantillonne m = F⁻¹(u).
  // Le maillage de 4096 segments est fixe : l'initialisation reste O(N).
  function buildInverseMassSampler() {
    try {
      const expression = els.density.value.trim();
      const minMass = Number(els.massMin.value);
      const maxMass = Number(els.massMax.value);
      if (!(minMass > 0) || !(maxMass > minMass)) {
        throw new Error("les bornes doivent vérifier 0 < min < max");
      }
      if (!expression || !/^[0-9A-Za-z_+\-*/().,\s]+$/.test(expression)) {
        throw new Error("expression vide ou caractère non autorisé");
      }

      const allowed = new Set([
        "m", "e", "pi", "exp", "log", "pow", "sqrt",
        "abs", "min", "max", "sin", "cos"
      ]);
      const identifiers = expression.match(/[A-Za-z_][A-Za-z_0-9]*/g) || [];
      if (identifiers.some((name) => !allowed.has(name))) {
        throw new Error("fonction inconnue dans f(m)");
      }

      const density = new Function(
        "m",
        `"use strict";
         const {exp,log,pow,sqrt,abs,min,max,sin,cos}=Math;
         const pi=Math.PI,e=Math.E;
         return (${expression});`
      );
      const bins = 4096;
      const xs = new Float64Array(bins + 1);
      const cdf = new Float64Array(bins + 1);
      const dx = (maxMass - minMass) / bins;
      let previous = Number(density(minMass));
      if (!Number.isFinite(previous) || previous < 0) {
        throw new Error("f(m) doit être finie et positive");
      }
      xs[0] = minMass;

      for (let i = 1; i <= bins; i++) {
        const x = minMass + i * dx;
        const value = Number(density(x));
        if (!Number.isFinite(value) || value < 0) {
          throw new Error("f(m) doit être finie et positive");
        }
        xs[i] = x;
        cdf[i] = cdf[i - 1] + (previous + value) * 0.5 * dx;
        previous = value;
      }

      const area = cdf[bins];
      if (!(area > 0)) throw new Error("l’intégrale de f(m) est nulle");
      for (let i = 1; i <= bins; i++) cdf[i] /= area;

      const invert = (u) => {
        let low = 0;
        let high = bins;
        while (high - low > 1) {
          const middle = (low + high) >> 1;
          if (cdf[middle] < u) low = middle;
          else high = middle;
        }
        const span = cdf[high] - cdf[low];
        const t = span > 1e-14 ? (u - cdf[low]) / span : 0;
        return xs[low] + (xs[high] - xs[low]) * t;
      };
      const inverseLut = new Float32Array(bins + 1);
      let inverseMean = 0;
      for (let index = 0; index <= bins; index++) {
        const mass = invert(index / bins);
        inverseLut[index] = mass;
        inverseMean += mass;
      }
      inverseMean /= inverseLut.length;

      els.densityStatus.classList.remove("error");
      els.densityStatus.value =
        "CDF numérique prête · inversion sur 4 096 segments";
      const sampler = () => invert(Math.random());
      sampler.inverseLut = inverseLut;
      sampler.inverseMean = inverseMean;
      sampler.minMass = minMass;
      sampler.maxMass = maxMass;
      return sampler;
    } catch (error) {
      els.densityStatus.classList.add("error");
      els.densityStatus.value = "Erreur · " + error.message;
      return null;
    }
  }

  function formatCosmicValue(value) {
    const absolute = Math.abs(value);
    if ((absolute > 0 && absolute < 0.001) || absolute >= 10000) {
      return value.toExponential(3);
    }
    return value.toFixed(3).replace(/\.?0+$/, "");
  }

  function updateExpansionStatus(force = false) {
    if (!expansionValid) return;
    const tick = Math.floor(cosmicTime * 10);
    if (!force && tick === lastExpansionStatusTick) return;
    lastExpansionStatusTick = tick;
    const scale = Number(expansionFunction(cosmicTime));
    if (!Number.isFinite(scale)) return;
    els.expansionStatus.classList.remove("error");
    els.expansionStatus.value =
      `Temps t = ${formatCosmicValue(cosmicTime)} · ` +
      `a(t) = ${formatCosmicValue(scale)}`;
  }

  function compileExpansionFunction() {
    try {
      const expression = els.expansion.value.trim();
      if (!expression || !/^[0-9A-Za-z_+\-*/().,\s]+$/.test(expression)) {
        throw new Error("expression vide ou caractère non autorisé");
      }
      const allowed = new Set([
        "t", "e", "pi", "exp", "log", "pow", "sqrt", "abs",
        "min", "max", "sin", "cos", "sinh", "cosh"
      ]);
      const identifiers = expression.match(/[A-Za-z_][A-Za-z_0-9]*/g) || [];
      if (identifiers.some((name) => !allowed.has(name))) {
        throw new Error("fonction inconnue dans a(t)");
      }
      const candidate = new Function(
        "t",
        `"use strict";
         const {exp,log,pow,sqrt,abs,min,max,sin,cos,sinh,cosh}=Math;
         const pi=Math.PI,e=Math.E;
         return (${expression});`
      );
      const probes = [0, cosmicTime, cosmicTime + 1, cosmicTime + 10];
      for (const time of probes) {
        const value = Number(candidate(time));
        if (!Number.isFinite(value)) {
          throw new Error("a(t) doit rester finie");
        }
      }
      expansionFunction = candidate;
      expansionValid = true;
      lastExpansionStatusTick = -1;
      updateExpansionStatus(true);
      return true;
    } catch (error) {
      // Une saisie invalide suspend l'expansion au lieu de continuer
      // silencieusement avec l'ancienne loi.
      expansionFunction = () => 1;
      expansionValid = false;
      els.expansionStatus.classList.add("error");
      els.expansionStatus.value =
        "Expansion suspendue · " + error.message;
      return false;
    }
  }

  function simulationStepSize(dt) {
    return Math.min(0.012, dt * 0.00038) * Number(els.speed.value);
  }

  function expansionScaleForStep(stepSize) {
    const nextTime = cosmicTime + stepSize;
    if (!expansionValid) {
      cosmicTime = nextTime;
      return 1;
    }
    const currentScale = Number(expansionFunction(cosmicTime));
    const nextScale = Number(expansionFunction(nextTime));
    if (
      !Number.isFinite(currentScale) ||
      !Number.isFinite(nextScale)
    ) {
      expansionFunction = () => 1;
      expansionValid = false;
      els.expansionStatus.classList.add("error");
      els.expansionStatus.value =
        "Expansion suspendue · a(t) n’est plus finie";
      return 1;
    }
    cosmicTime = nextTime;
    updateExpansionStatus();
    // a(t)=0 et le changement de signe sont singuliers pour le rapport
    // a(t+dt)/a(t). Une traversée est représentée par une inversion de
    // parité bornée, sans projection numérique vers l'infini.
    if (currentScale === 0 && nextScale === 0) return 1;
    if (currentScale === 0) return Math.sign(nextScale) || 1;
    if (nextScale === 0) return 1;
    if (Math.sign(currentScale) !== Math.sign(nextScale)) return -1;
    return nextScale / currentScale;
  }

  function updateGravityParameters() {
    gravityStrength = Number(els.gravity.value);
    gravityExponent = Number(els.gravityPower.value);
    const softening = Number(els.softening.value);
    gravitySoftening2 = softening * softening;
    els.gravityValue.value = gravityStrength.toFixed(2);
    els.gravityPowerValue.value =
      `1 / r${gravityExponent === 2 ? "²" : "^" + gravityExponent.toFixed(2)}`;
    els.softeningValue.value = softening.toFixed(3);
    els.gravityStatus.value =
      gravityStrength < 0
        ? "Gravité répulsive · paramètres appliqués en direct"
        : gravityStrength === 0
          ? "Gravité désactivée · mouvement inertiel"
          : "Gravité attractive · paramètres appliqués en direct";
  }

  function gravityKernel(q) {
    if (gravityExponent === 2) {
      const inverse = 1 / Math.sqrt(q);
      return inverse * inverse * inverse;
    }
    return Math.pow(q, -(gravityExponent + 1) * 0.5);
  }

  function resetBodies() {
    const sampler = buildInverseMassSampler();
    if (!sampler) return;
    cosmicTime = 0;
    lastExpansionStatusTick = -1;
    updateExpansionStatus(true);
    let n = selectedBodyCount();
    if (els.renderer.value !== "gpu-compute" && n > CPU_MAX_BODIES) {
      n = CPU_MAX_BODIES;
      els.count.value = "18";
      els.rendererStatus.value =
        "Limite CPU · 262 144 corps · choisissez GPU Compute pour monter";
    }
    simulationCount = n;
    if (els.renderer.value === "gpu-compute" && n > CPU_MAX_BODIES) {
      bodies = [];
      gpuSeedConfig = {
        massLut: sampler.inverseLut,
        inverseMean: sampler.inverseMean,
        seed: Math.floor(Math.random() * 0xffffffff),
        gravitation: gravityStrength,
        gravityPower: gravityExponent
      };
      els.bodyStat.value = n.toLocaleString("fr-FR");
      els.countValue.value = n.toLocaleString("fr-FR");
      els.densityStatus.value =
        `Échantillonnage inverse GPU · [${sampler.minMass.toPrecision(3)}, ` +
        `${sampler.maxMass.toPrecision(3)}] · génération sans objets JS`;
      impulse = null;
      fmmMsSmooth = 0;
      renderMsSmooth = 0;
      fpsSmooth = 60;
      treeCache = null;
      void setRendererMode();
      return;
    }

    gpuSeedConfig = null;
    bodies = new Array(n);
    const armCount = 3;
    let totalMass = 0;
    let sampledMin = Infinity;
    let sampledMax = -Infinity;

    for (let i = 0; i < n; i++) {
      const core = i < Math.max(8, n * 0.025);
      const radius = core
        ? Math.random() * 0.06
        : Math.pow(Math.random(), 0.58) * 0.83 + 0.035;
      const arm = i % armCount;
      const angle =
        arm * Math.PI * 2 / armCount +
        radius * 7.8 +
        randn() * (0.13 + radius * 0.13);
      const flatten = 0.73;
      const x = Math.cos(angle) * radius + randn() * 0.012;
      const y = Math.sin(angle) * radius * flatten + randn() * 0.012;
      const enclosed = 0.16 + Math.min(1, radius * 1.35);
      const orbital =
        Math.sqrt(
          Math.abs(gravityStrength) *
          enclosed *
          Math.pow(Math.max(radius, 0.035), 1 - gravityExponent)
        ) * 0.49;
      const vx = -Math.sin(angle) * orbital + randn() * 0.018;
      const vy = Math.cos(angle) * orbital * flatten + randn() * 0.018;
      const mass = sampler();
      totalMass += mass;
      if (mass < sampledMin) sampledMin = mass;
      if (mass > sampledMax) sampledMax = mass;

      bodies[i] = {
        x, y, vx, vy,
        ax: 0, ay: 0,
        mass,
        color: i % palette.length,
        cellIndex: 0
      };
    }

    // La masse totale reste 1 pour que la dynamique ne dépende que de la
    // forme de la distribution, pas de son unité choisie.
    for (const body of bodies) body.mass /= totalMass;

    els.bodyStat.value = n.toLocaleString("fr-FR");
    els.countValue.value = n.toLocaleString("fr-FR");
    els.densityStatus.value =
      `Échantillonnage inverse · brut [${sampledMin.toPrecision(3)}, ` +
      `${sampledMax.toPrecision(3)}] · Σm normalisée à 1`;
    impulse = null;
    fmmMsSmooth = 0;
    renderMsSmooth = 0;
    fpsSmooth = 60;
    treeCache = null;
    if (els.renderer.value === "gpu-compute") {
      if (!gpuComputeEngine || !gpuComputeEngine.reset(bodies)) {
        void setRendererMode();
      }
    }
  }

  function isTorus() {
    return torusActive;
  }

  function minimumImage(delta) {
    if (!isTorus()) return delta;
    if (delta > BOX_SIZE * 0.5) return delta - BOX_SIZE;
    if (delta < -BOX_SIZE * 0.5) return delta + BOX_SIZE;
    return delta;
  }

  function computeFMM() {
    const startedAt = performance.now();
    const n = bodies.length;
    boundaryMode = els.boundary.value;
    torusActive = boundaryMode === 'torus';
    const result = cpuFMM.compute(bodies, { boundary:boundaryMode,order:Number(els.accuracy.value),gravity:gravityStrength,exponent:gravityExponent,softening:Math.sqrt(gravitySoftening2) });
    const {levels, depth} = result;
    let ops = result.operations;
    if (impulse) {
      for (const b of bodies) {
        const rx = minimumImage(impulse.x - b.x);
        const ry = minimumImage(impulse.y - b.y);
        const q = rx * rx + ry * ry + 0.012;
        const f = impulse.force / (q * Math.sqrt(q));
        b.ax += rx * f;
        b.ay += ry * f;
      }
      impulse.force *= 0.91;
      if (impulse.force < 0.003) impulse = null;
    }

    treeSnapshot = { levels, depth };
    operationCount = Math.max(1, ops);
    let activeCells = 0;
    for (const level of levels) {
      for (const cell of level.cells) {
        if (cell.mass > 0) activeCells++;
      }
    }
    els.cellStat.value = activeCells.toLocaleString("fr-FR");
    const nlogOperations = n * Math.log2(Math.max(2, n));
    const gain = n * Math.max(1, n - 1) / operationCount;
    const nlogRatio = nlogOperations / operationCount;
    els.nlogStat.value =
      nlogRatio < 1
        ? nlogRatio.toFixed(2) + "×"
        : nlogRatio < 10
          ? nlogRatio.toFixed(1) + "×"
          : Math.round(nlogRatio) + "×";
    els.gainStat.value =
      gain < 10 ? gain.toFixed(1) + "×" : Math.round(gain) + "×";
    const elapsed = performance.now() - startedAt;
    fmmMsSmooth = fmmMsSmooth
      ? fmmMsSmooth * 0.9 + elapsed * 0.1
      : elapsed;
  }

  function applyBoundary(body) {
    const mode = boundaryMode;
    if (mode === "open") return;
    if (!Number.isFinite(body.x) || !Number.isFinite(body.vx)) {
      body.x = 0;
      body.vx = 0;
    }
    if (!Number.isFinite(body.y) || !Number.isFinite(body.vy)) {
      body.y = 0;
      body.vy = 0;
    }
    if (mode === "torus") {
      body.x =
        BOX_MIN + ((body.x - BOX_MIN) % BOX_SIZE + BOX_SIZE) % BOX_SIZE;
      body.y =
        BOX_MIN + ((body.y - BOX_MIN) % BOX_SIZE + BOX_SIZE) % BOX_SIZE;
      return;
    }

    while (body.x < BOX_MIN || body.x > BOX_MAX) {
      if (body.x < BOX_MIN) body.x = 2 * BOX_MIN - body.x;
      else body.x = 2 * BOX_MAX - body.x;
      body.vx *= -1;
    }
    while (body.y < BOX_MIN || body.y > BOX_MAX) {
      if (body.y < BOX_MIN) body.y = 2 * BOX_MIN - body.y;
      else body.y = 2 * BOX_MAX - body.y;
      body.vy *= -1;
    }
  }

  function step(dt, expansionScale) {
    computeFMM();
    const h = simulationStepSize(dt);
    for (const b of bodies) {
      b.vx += b.ax * h;
      b.vy += b.ay * h;
      b.vx *= 0.99996;
      b.vy *= 0.99996;
      b.x = (b.x + b.vx * h) * expansionScale;
      b.y = (b.y + b.vy * h) * expansionScale;
      applyBoundary(b);
    }
  }

  function worldToScreen(x, y) {
    const scale = Math.min(width, height) * 0.48 * cameraZoom;
    return [width * 0.5 + x * scale, height * 0.5 + y * scale];
  }

  function adaptiveViewBrightness() {
    return Math.max(0.55, Math.min(8, cameraZoom ** 0.75));
  }

  function adaptivePointScale() {
    return Math.max(0.8, Math.min(2, cameraZoom ** 0.2));
  }

  function drawBoundary() {
    if (els.boundary.value === "open") return;
    const topLeft = worldToScreen(BOX_MIN, BOX_MIN);
    const bottomRight = worldToScreen(BOX_MAX, BOX_MAX);
    ctx.save();
    ctx.strokeStyle =
      els.boundary.value === "torus"
        ? "rgba(157,114,255,.72)"
        : "rgba(215,255,82,.62)";
    ctx.lineWidth = 1.2;
    if (els.boundary.value === "torus") ctx.setLineDash([7, 5]);
    ctx.strokeRect(
      topLeft[0],
      topLeft[1],
      bottomRight[0] - topLeft[0],
      bottomRight[1] - topLeft[1]
    );
    ctx.restore();
  }

  function drawImpulse() {
    if (!impulse) return;
    const [x, y] = worldToScreen(impulse.x, impulse.y);
    ctx.strokeStyle =
      `rgba(215,255,82,${Math.min(0.8, impulse.force * 5)})`;
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.arc(x, y, 12 + impulse.force * 18, 0, Math.PI * 2);
    ctx.stroke();
  }

  function drawTree() {
    if (!els.tree.checked) return;
    if (els.renderer.value === "gpu-compute" && gpuComputeEngine) {
      ctx.save();
      ctx.lineWidth = 0.7;
      const minLevel = Math.max(1, gpuComputeEngine.depth - 3);
      for (
        let level = minLevel;
        level <= gpuComputeEngine.depth;
        level++
      ) {
        const side = 1 << level;
        ctx.strokeStyle =
          `rgba(98,217,255,${0.08 + (level - minLevel) * 0.035})`;
        ctx.beginPath();
        for (let index = 0; index <= side; index++) {
          const coordinate = BOX_MIN + BOX_SIZE * index / side;
          const verticalA = worldToScreen(coordinate, BOX_MIN);
          const verticalB = worldToScreen(coordinate, BOX_MAX);
          const horizontalA = worldToScreen(BOX_MIN, coordinate);
          const horizontalB = worldToScreen(BOX_MAX, coordinate);
          ctx.moveTo(verticalA[0], verticalA[1]);
          ctx.lineTo(verticalB[0], verticalB[1]);
          ctx.moveTo(horizontalA[0], horizontalA[1]);
          ctx.lineTo(horizontalB[0], horizontalB[1]);
        }
        ctx.stroke();
      }
      ctx.restore();
      return;
    }
    if (!treeSnapshot) return;
    const { levels, depth } = treeSnapshot;
    ctx.save();
    ctx.lineWidth = 0.85;
    for (let level = Math.max(1, depth - 3); level <= depth; level++) {
      ctx.strokeStyle =
        `rgba(98,217,255,${0.11 + (level - depth + 3) * 0.045})`;
      for (const cell of levels[level].cells) {
        if (cell.mass <= 0) continue;
        const a = worldToScreen(
          cell.cx - cell.size / 2,
          cell.cy - cell.size / 2
        );
        const b = worldToScreen(
          cell.cx + cell.size / 2,
          cell.cy + cell.size / 2
        );
        ctx.strokeRect(a[0], a[1], b[0] - a[0], b[1] - a[1]);
      }
    }
    ctx.restore();
  }

  function draw() {
    const startedAt = performance.now();
    ctx.fillStyle = "#05070b";
    ctx.fillRect(0, 0, width, height);

    const gradient = ctx.createRadialGradient(
      width * 0.5,
      height * 0.48,
      0,
      width * 0.5,
      height * 0.48,
      Math.max(width, height) * 0.58
    );
    gradient.addColorStop(0, "rgba(27,24,45,.34)");
    gradient.addColorStop(0.5, "rgba(9,12,20,.12)");
    gradient.addColorStop(1, "rgba(0,0,0,.3)");
    ctx.fillStyle = gradient;
    ctx.fillRect(0, 0, width, height);

    drawTree();
    drawBoundary();
    const viewBrightness = adaptiveViewBrightness();
    const pointScale = adaptivePointScale();
    if (bodies.length > 4096) {
      const pixels = particleImage.data;
      const rasterWidth = particleCanvas.width;
      const rasterHeight = particleCanvas.height;
      pixels.fill(0);
      for (const body of bodies) {
        const point = worldToScreen(body.x, body.y);
        const x = Math.round(point[0]);
        const y = Math.round(point[1]);
        if (x < 0 || y < 0 || x >= rasterWidth || y >= rasterHeight) continue;
        const index = (y * rasterWidth + x) * 4;
        const color = palette[body.color];
        pixels[index] = Math.min(255, pixels[index] + color[0]);
        pixels[index + 1] = Math.min(255, pixels[index + 1] + color[1]);
        pixels[index + 2] = Math.min(255, pixels[index + 2] + color[2]);
        pixels[index + 3] = 235;
      }
      particleCtx.putImageData(particleImage, 0, 0);
      ctx.save();
      ctx.globalCompositeOperation = "lighter";
      ctx.globalAlpha = Math.min(1, 0.3 * viewBrightness);
      ctx.filter = `blur(${(3 * pointScale).toFixed(2)}px)`;
      ctx.drawImage(particleCanvas, 0, 0, width, height);
      ctx.filter = "none";
      ctx.globalAlpha = Math.min(1, 0.92 * Math.sqrt(viewBrightness));
      ctx.drawImage(particleCanvas, 0, 0, width, height);
      ctx.restore();
    } else {
      ctx.save();
      ctx.globalCompositeOperation = "lighter";
      const radius =
        (bodies.length > 1200 ? 0.72 : bodies.length > 600 ? 0.9 : 1.15) *
        pointScale;
      for (const body of bodies) {
        const [x, y] = worldToScreen(body.x, body.y);
        if (x < -5 || y < -5 || x > width + 5 || y > height + 5) continue;
        const color = palette[body.color];
        ctx.fillStyle =
          `rgba(${color[0]},${color[1]},${color[2]},` +
          `${Math.min(0.9, 0.22 * viewBrightness)})`;
        ctx.beginPath();
        ctx.arc(x, y, radius * 3.4, 0, Math.PI * 2);
        ctx.fill();
        ctx.fillStyle =
          `rgba(${color[0]},${color[1]},${color[2]},` +
          `${Math.min(1, 0.88 * Math.sqrt(viewBrightness))})`;
        ctx.beginPath();
        ctx.arc(x, y, radius, 0, Math.PI * 2);
        ctx.fill();
      }
      ctx.restore();
    }

    drawImpulse();
    const elapsed = performance.now() - startedAt;
    renderMsSmooth = renderMsSmooth
      ? renderMsSmooth * 0.9 + elapsed * 0.1
      : elapsed;
    els.timingStat.value =
      `${fmmMsSmooth.toFixed(1)} / ${renderMsSmooth.toFixed(1)} ms`;
  }

  function drawOverlay(baseRenderMs = 0) {
    const startedAt = performance.now();
    ctx.clearRect(0, 0, width, height);
    drawTree();
    drawBoundary();
    drawImpulse();
    const elapsed = baseRenderMs + performance.now() - startedAt;
    renderMsSmooth = renderMsSmooth
      ? renderMsSmooth * 0.9 + elapsed * 0.1
      : elapsed;
    els.timingStat.value =
      `${fmmMsSmooth.toFixed(1)} / ${renderMsSmooth.toFixed(1)} ms`;
  }

  async function createWebGPURenderer() {
    if (!navigator.gpu) {
      throw new Error("WebGPU indisponible dans ce navigateur ou ce contexte");
    }
    const adapter = await navigator.gpu.requestAdapter({
      powerPreference: "high-performance"
    });
    if (!adapter) throw new Error("aucun adaptateur WebGPU disponible");
    const device = await adapter.requestDevice();
    const context = gpuCanvas.getContext("webgpu");
    if (!context) throw new Error("canvas WebGPU indisponible");
    const format = navigator.gpu.getPreferredCanvasFormat();
    const shader = device.createShaderModule({
      label: "Orbit particle shader",
      code: `
        struct Uniforms {
          worldScale: vec2<f32>,
          pointScale: vec2<f32>,
          brightness: f32,
          pad0: f32,
          pad1: f32,
          pad2: f32,
        }
        @group(0) @binding(0) var<uniform> uniforms: Uniforms;

        struct VertexOutput {
          @builtin(position) position: vec4<f32>,
          @location(0) uv: vec2<f32>,
          @location(1) color: vec3<f32>,
          @location(2) intensity: f32,
        }

        fn palette(index: u32) -> vec3<f32> {
          switch index {
            case 0u: { return vec3<f32>(0.243, 0.851, 1.0); }
            case 1u: { return vec3<f32>(0.616, 0.447, 1.0); }
            case 2u: { return vec3<f32>(0.843, 1.0, 0.322); }
            default: { return vec3<f32>(1.0, 0.616, 0.333); }
          }
        }

        @vertex
        fn vertexMain(
          @builtin(vertex_index) vertexIndex: u32,
          @location(0) particle: vec4<f32>
        ) -> VertexOutput {
          let corners = array<vec2<f32>, 6>(
            vec2<f32>(-1.0, -1.0),
            vec2<f32>( 1.0, -1.0),
            vec2<f32>(-1.0,  1.0),
            vec2<f32>(-1.0,  1.0),
            vec2<f32>( 1.0, -1.0),
            vec2<f32>( 1.0,  1.0)
          );
          let corner = corners[vertexIndex];
          let center = vec2<f32>(
            particle.x * uniforms.worldScale.x,
            -particle.y * uniforms.worldScale.y
          );
          var output: VertexOutput;
          output.position = vec4<f32>(
            center + corner * uniforms.pointScale,
            0.0,
            1.0
          );
          output.uv = corner;
          output.color = palette(u32(particle.z + 0.5));
          output.intensity =
            clamp(0.72 + 0.18 * particle.w, 0.65, 1.35) *
            uniforms.brightness;
          return output;
        }

        @fragment
        fn fragmentMain(input: VertexOutput) -> @location(0) vec4<f32> {
          let distance = length(input.uv);
          if (distance > 1.0) { discard; }
          let halo = pow(max(0.0, 1.0 - distance), 2.2);
          let core = smoothstep(0.32, 0.0, distance);
          let alpha = clamp((halo * 0.42 + core * 0.92) * input.intensity, 0.0, 1.0);
          return vec4<f32>(input.color, alpha);
        }
      `
    });

    const pipeline = await device.createRenderPipelineAsync({
      label: "Orbit WebGPU particle pipeline",
      layout: "auto",
      vertex: {
        module: shader,
        entryPoint: "vertexMain",
        buffers: [{
          arrayStride: 16,
          stepMode: "instance",
          attributes: [{
            shaderLocation: 0,
            offset: 0,
            format: "float32x4"
          }]
        }]
      },
      fragment: {
        module: shader,
        entryPoint: "fragmentMain",
        targets: [{
          format,
          blend: {
            color: {
              srcFactor: "src-alpha",
              dstFactor: "one",
              operation: "add"
            },
            alpha: {
              srcFactor: "one",
              dstFactor: "one",
              operation: "add"
            }
          }
        }]
      },
      primitive: { topology: "triangle-list" }
    });

    const uniformBuffer = device.createBuffer({
      label: "Orbit viewport uniforms",
      size: 32,
      usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST
    });
    const bindGroup = device.createBindGroup({
      layout: pipeline.getBindGroupLayout(0),
      entries: [{ binding: 0, resource: { buffer: uniformBuffer } }]
    });
    let vertexBuffer = null;
    let capacity = 0;

    const renderer = {
      device,
      configure() {
        context.configure({
          device,
          format,
          alphaMode: "opaque"
        });
      },
      ensureCapacity(count) {
        if (count <= capacity) return;
        capacity = 1;
        while (capacity < count) capacity *= 2;
        if (vertexBuffer) vertexBuffer.destroy();
        vertexBuffer = device.createBuffer({
          label: "Orbit particle instances",
          size: capacity * 16,
          usage: GPUBufferUsage.VERTEX | GPUBufferUsage.COPY_DST
        });
        gpuParticleData = new Float32Array(capacity * 4);
      },
      render(particles) {
        const startedAt = performance.now();
        this.ensureCapacity(particles.length);
        for (let index = 0; index < particles.length; index++) {
          const body = particles[index];
          const offset = index * 4;
          gpuParticleData[offset] = body.x;
          gpuParticleData[offset + 1] = body.y;
          gpuParticleData[offset + 2] = body.color;
          gpuParticleData[offset + 3] =
            Math.min(3.5, body.mass * particles.length);
        }
        device.queue.writeBuffer(
          vertexBuffer,
          0,
          gpuParticleData.buffer,
          0,
          particles.length * 16
        );
        const minDimension = Math.min(width, height);
        const radius =
          (particles.length > 65536 ? 2.2 : 3.2) *
          adaptivePointScale();
        device.queue.writeBuffer(
          uniformBuffer,
          0,
          new Float32Array([
            0.96 * minDimension / width * cameraZoom,
            0.96 * minDimension / height * cameraZoom,
            2 * radius * dpr / gpuCanvas.width,
            2 * radius * dpr / gpuCanvas.height,
            adaptiveViewBrightness(),
            0,
            0,
            0
          ])
        );

        const encoder = device.createCommandEncoder({
          label: "Orbit frame encoder"
        });
        const pass = encoder.beginRenderPass({
          colorAttachments: [{
            view: context.getCurrentTexture().createView(),
            clearValue: { r: 0.019, g: 0.027, b: 0.043, a: 1 },
            loadOp: "clear",
            storeOp: "store"
          }]
        });
        pass.setPipeline(pipeline);
        pass.setBindGroup(0, bindGroup);
        pass.setVertexBuffer(0, vertexBuffer);
        pass.draw(6, particles.length);
        pass.end();
        device.queue.submit([encoder.finish()]);
        return performance.now() - startedAt;
      }
    };
    renderer.configure();
    device.lost.then((info) => {
      if (gpuRenderer !== renderer) return;
      gpuRenderer = null;
      els.renderer.value = "cpu";
      gpuCanvas.classList.remove("active");
      els.rendererStatus.classList.add("error");
      els.rendererStatus.value =
        "GPU perdu · retour CPU" + (info.message ? " · " + info.message : "");
    });
    return renderer;
  }

  async function setRendererMode() {
    const request = ++rendererRequest;
    const mode = els.renderer.value;
    els.rendererStatus.classList.remove("error");

    if (mode !== "gpu-compute" && gpuComputeEngine) {
      gpuComputeEngine.destroy();
      gpuComputeEngine = null;
    }

    if (mode !== "gpu-compute" && simulationCount > CPU_MAX_BODIES) {
      els.count.value = "18";
      resetBodies();
    }

    if (mode === "cpu") {
      gpuCanvas.classList.remove("active");
      els.rendererStatus.value = "Canvas 2D actif";
      return;
    }

    if (mode === "gpu-compute") {
      if (
        gpuComputeEngine?.count === simulationCount &&
        !gpuSeedConfig
      ) {
        gpuCanvas.classList.add("active");
        els.rendererStatus.value =
          "WebGPU Compute actif · FMM, intégration et rendu sur GPU";
        return;
      }
      if (gpuComputeEngine) {
        gpuComputeEngine.destroy();
        gpuComputeEngine = null;
      }

      els.rendererStatus.value = "Compilation du FMM WebGPU…";
      try {
        if (typeof window.createOrbitGPUComputeEngine !== "function") {
          throw new Error("module GPU Compute introuvable");
        }
        const initialization = gpuSeedConfig;
        const engine = await window.createOrbitGPUComputeEngine({
          canvas: gpuCanvas,
          particles: bodies,
          count: simulationCount,
          initialization,
          width,
          height,
          dpr,
          onDeviceLost(info) {
            if (request !== rendererRequest) return;
            gpuComputeEngine = null;
            els.renderer.value = "cpu";
            if (simulationCount > CPU_MAX_BODIES) {
              els.count.value = "18";
              resetBodies();
            }
            gpuCanvas.classList.remove("active");
            els.rendererStatus.classList.add("error");
            els.rendererStatus.value =
              "GPU perdu · retour CPU" +
              (info.message ? " · " + info.message : "");
          }
        });
        if (
          request !== rendererRequest ||
          els.renderer.value !== "gpu-compute"
        ) {
          engine.destroy();
          return;
        }
        gpuComputeEngine = engine;
        if (gpuSeedConfig === initialization) gpuSeedConfig = null;
        gpuCanvas.classList.add("active");
        els.cellStat.value =
          "≤ " + engine.totalCells.toLocaleString("fr-FR");
        els.rendererStatus.value =
          "WebGPU Compute actif · FMM, intégration et rendu sur GPU";
      } catch (error) {
        if (request !== rendererRequest) return;
        gpuComputeEngine = null;
        if (error.maximumBodies > CPU_MAX_BODIES) {
          const supportedExponent =
            Math.floor(Math.log2(error.maximumBodies) * 4) / 4;
          els.count.value = String(supportedExponent);
          els.countValue.value =
            selectedBodyCount().toLocaleString("fr-FR");
          els.rendererStatus.value =
            "Ajustement à la limite mémoire de cette carte…";
          resetBodies();
          return;
        }
        els.renderer.value = "cpu";
        if (simulationCount > CPU_MAX_BODIES) {
          els.count.value = "18";
          resetBodies();
        }
        gpuCanvas.classList.remove("active");
        els.rendererStatus.classList.add("error");
        els.rendererStatus.value = "Repli CPU · " + error.message;
      }
      return;
    }

    els.rendererStatus.value = "Initialisation WebGPU…";
    try {
      if (!gpuRenderer) {
        gpuInitPromise ||= createWebGPURenderer();
        gpuRenderer = await gpuInitPromise;
      }
      if (request !== rendererRequest || els.renderer.value !== "gpu") return;
      gpuCanvas.classList.add("active");
      els.rendererStatus.value =
        "WebGPU actif · solveur FMM toujours sur CPU";
    } catch (error) {
      if (request !== rendererRequest) return;
      gpuInitPromise = null;
      gpuRenderer = null;
      els.renderer.value = "cpu";
      gpuCanvas.classList.remove("active");
      els.rendererStatus.classList.add("error");
      els.rendererStatus.value = "Repli CPU · " + error.message;
    }
  }

  function frame(now) {
    const rawDt = Math.max(1, now - lastTime);
    const dt = Math.min(40, rawDt);
    lastTime = now;
    fpsSmooth =
      fpsSmooth * 0.92 + 1000 / rawDt * 0.08;

    const rendererMode = els.renderer.value;
    if (rendererMode === "gpu-compute") {
      if (gpuComputeEngine) {
        try {
          let performanceSample = {
            submitMs: 0,
            gpuMs: gpuComputeEngine.lastGpuMs
          };
          if (running) {
            const expansionScale =
              expansionScaleForStep(simulationStepSize(dt));
            performanceSample = gpuComputeEngine.stepAndRender({
              dt,
              speed: Number(els.speed.value),
              boundary: els.boundary.value,
              order: Number(els.accuracy.value),
              impulse,
              expansionScale,
              zoom: cameraZoom,
              gravitation: gravityStrength,
              gravityPower: gravityExponent,
              softening2: gravitySoftening2
            });
            if (impulse) {
              impulse.force *= 0.91;
              if (impulse.force < 0.003) impulse = null;
            }
          }
          drawOverlay(performanceSample.submitMs);

          const n = gpuComputeEngine.count;
          const estimatedLinearWork =
            n * 10 + gpuComputeEngine.totalCells * 27;
          const nlogRatio =
            n * Math.log2(Math.max(2, n)) / estimatedLinearWork;
          const quadraticRatio =
            n * Math.max(1, n - 1) / estimatedLinearWork;
          els.cellStat.value =
            "≤ " +
            gpuComputeEngine.totalCells.toLocaleString("fr-FR");
          els.nlogStat.value =
            (nlogRatio < 10
              ? nlogRatio.toFixed(2)
              : Math.round(nlogRatio)) + "×";
          els.gainStat.value =
            (quadraticRatio < 10
              ? quadraticRatio.toFixed(1)
              : Math.round(quadraticRatio)) + "×";
          const frameInterval = 1000 / Math.max(1, fpsSmooth);
          els.timingStat.value = performanceSample.gpuMs > frameInterval * 2
            ? `${frameInterval.toFixed(1)} ms/image · ` +
              `file ${performanceSample.gpuMs.toFixed(0)} ms`
            : performanceSample.gpuMs > 0
              ? `${performanceSample.gpuMs.toFixed(1)} ms GPU · ` +
                `${performanceSample.submitMs.toFixed(2)} ms pilotage`
              : `mesure GPU… · ` +
                `${performanceSample.submitMs.toFixed(2)} ms pilotage`;
        } catch (error) {
          rendererRequest++;
          gpuComputeEngine.destroy();
          gpuComputeEngine = null;
          els.renderer.value = "cpu";
          if (simulationCount > CPU_MAX_BODIES) {
            els.count.value = "18";
            resetBodies();
          }
          gpuCanvas.classList.remove("active");
          els.rendererStatus.classList.add("error");
          els.rendererStatus.value = "Repli CPU · " + error.message;
          draw();
        }
      } else {
        draw();
      }
    } else {
      if (running) {
        const expansionScale =
          expansionScaleForStep(simulationStepSize(dt));
        step(dt, expansionScale);
      }
    }

    if (rendererMode === "gpu" && gpuRenderer) {
      try {
        const gpuMs = gpuRenderer.render(bodies);
        drawOverlay(gpuMs);
      } catch (error) {
        els.renderer.value = "cpu";
        gpuCanvas.classList.remove("active");
        els.rendererStatus.classList.add("error");
        els.rendererStatus.value = "Repli CPU · " + error.message;
        draw();
      }
    } else if (rendererMode !== "gpu-compute") {
      draw();
    }
    els.fpsStat.value = Math.round(fpsSmooth).toString();
    requestAnimationFrame(frame);
  }

  function setImpulse(event) {
    const rect = canvas.getBoundingClientRect();
    const scale = Math.min(width, height) * 0.48 * cameraZoom;
    impulse = {
      x: (event.clientX - rect.left - width * 0.5) / scale,
      y: (event.clientY - rect.top - height * 0.5) / scale,
      force: 0.16
    };
  }

  els.count.addEventListener("input", () => {
    els.countValue.value =
      selectedBodyCount().toLocaleString("fr-FR");
  });
  els.count.addEventListener("change", resetBodies);
  els.speed.addEventListener("input", () => {
    els.speedValue.value =
      Number(els.speed.value).toFixed(1) + "×";
  });
  els.zoom.addEventListener("input", () => {
    cameraZoom = 2 ** Number(els.zoom.value);
    els.zoomValue.value = cameraZoom.toFixed(2) + "×";
  });
  els.accuracy.addEventListener("input", () => {
    els.accuracyValue.value = "p = " + els.accuracy.value;
  });
  els.boundary.addEventListener("change", () => {
    boundaryMode = els.boundary.value;
    torusActive = boundaryMode === "torus";
    treeCache = null;
    for (const body of bodies) applyBoundary(body);
  });
  els.renderer.addEventListener("change", () => {
    if (
      els.renderer.value !== "gpu-compute" &&
      selectedBodyCount() > CPU_MAX_BODIES
    ) {
      els.count.value = "18";
      els.countValue.value =
        CPU_MAX_BODIES.toLocaleString("fr-FR");
    }
    void setRendererMode();
  });
  for (const input of [els.density, els.massMin, els.massMax]) {
    input.addEventListener("input", () => {
      els.densityStatus.classList.remove("error");
      els.densityStatus.value =
        "Paramètres modifiés · cliquez sur Ré-échantillonner";
    });
  }
  els.expansion.addEventListener("input", compileExpansionFunction);
  for (const input of [
    els.gravity,
    els.gravityPower,
    els.softening
  ]) {
    input.addEventListener("input", updateGravityParameters);
  }
  els.reset.addEventListener("click", resetBodies);
  els.pause.addEventListener("click", () => {
    running = !running;
    els.pause.textContent = running ? "Pause" : "Reprendre";
    els.pause.setAttribute("aria-pressed", String(!running));
  });
  canvas.addEventListener("pointerdown", setImpulse);
  canvas.addEventListener("wheel", (event) => {
    event.preventDefault();
    const direction = event.deltaY < 0 ? 1 : -1;
    const exponent = Math.max(
      Number(els.zoom.min),
      Math.min(Number(els.zoom.max), Number(els.zoom.value) + direction * 0.12)
    );
    els.zoom.value = String(exponent);
    cameraZoom = 2 ** exponent;
    els.zoomValue.value = cameraZoom.toFixed(2) + "×";
  }, { passive: false });
  window.addEventListener("resize", resize);
  document.addEventListener("visibilitychange", () => {
    lastTime = performance.now();
  });

  resize();
  updateGravityParameters();
  compileExpansionFunction();
  resetBodies();
  requestAnimationFrame(frame);
})();
