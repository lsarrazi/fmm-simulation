(() => {
  "use strict";
  const $ = id => document.getElementById(id);
  const canvas = $("matter"), ctx = canvas.getContext("2d");
  const imageCanvas = document.createElement("canvas"), imageContext = imageCanvas.getContext("2d");
  let solver, image, running = !new URLSearchParams(location.search).has("paused"), last = performance.now(), reportAt = 0, ms = 0, referenceDensity;
  let identities, candidateIds = [], visible = [], selected = null, selectedState = null;
  let observation = { x: 0, y: 0, locked: false }, viewport, pendingCount = 0, replayCursor = 0;
  const initialControls = ["resolution", "initial", "sigma", "fluctuation", "radius", "vmax", "population"];
  function reset() {
    const [n, nv] = $("resolution").value.split(",").map(Number);
    solver = new KineticMatter({ n, nv, vmax: Number($("vmax").value), sigma: $("initial").value === "cold" ? 0 : Number($("sigma").value), fluctuation: Number($("fluctuation").value), radius: Number($("radius").value) });
    identities = new LazyMatterIdentities(solver, { population: 10n ** BigInt($("population").value) });
    solver.onTransport = record => identities.record(record);
    candidateIds = identities.candidates(Number($("probes").value));
    selected = null; selectedState = null; visible = []; replayCursor = 0;
    $("identityDetail").textContent = "Cliquez sur un point ou recherchez un identifiant.";
    $("pause").disabled = false;
    $("pause").textContent = running ? "Pause" : "Reprendre";
    $("pause").setAttribute("aria-pressed", String(!running));
    referenceDensity = Math.max(...solver.mass) / solver.dx ** 2;
    imageCanvas.width = imageCanvas.height = n;
    image = imageContext.createImageData(n, n);
    $("pending").hidden = true;
    last = performance.now(); reportAt = 0;
    sample(); report(); draw();
  }
  function draw() {
    const rect = canvas.getBoundingClientRect(), dpr = Math.min(2, devicePixelRatio || 1);
    const w = Math.round(rect.width * dpr), h = Math.round(rect.height * dpr);
    if (canvas.width !== w || canvas.height !== h) { canvas.width = w; canvas.height = h; }
    for (let cell = 0; cell < solver.mass.length; cell++) {
      const density = solver.mass[cell] / solver.dx ** 2;
      const t = Math.max(0, Math.min(1, Math.log1p(density / referenceDensity) / Math.log(21)));
      const stops = [[7, 10, 20], [24, 76, 145], [67, 185, 203], [170, 240, 208], [255, 236, 170]];
      const v = t * 4, k = Math.min(3, Math.floor(v)), u = v - k;
      const x = cell % solver.n, y = Math.floor(cell / solver.n);
      const out = ((solver.n - 1 - y) * solver.n + x) * 4;
      for (let c = 0; c < 3; c++) image.data[out + c] = stops[k][c] * (1 - u) + stops[k + 1][c] * u;
      image.data[out + 3] = 255;
    }
    imageContext.putImageData(image, 0, 0);
    ctx.fillStyle = "#07090e"; ctx.fillRect(0, 0, w, h);
    const size = Math.min(w, h) * 0.9, left = (w - size) / 2, top = (h - size) / 2;
    viewport = { size, left, top, dpr };
    ctx.imageSmoothingEnabled = true;
    ctx.drawImage(imageCanvas, left, top, size, size);
    ctx.strokeStyle = "#374153"; ctx.lineWidth = dpr; ctx.strokeRect(left, top, size, size);
    ctx.fillStyle = "#8490a5"; ctx.font = `${11 * dpr}px monospace`;
    ctx.fillText("−1.5", left, top + size + 17 * dpr); ctx.fillText("+1.5", left + size - 30 * dpr, top + size + 17 * dpr);
    const point = p => [left + (p.x + solver.extent) / (2 * solver.extent) * size, top + (solver.extent - p.y) / (2 * solver.extent) * size];
    ctx.save(); ctx.beginPath(); ctx.rect(left, top, size, size); ctx.clip();
    if ($("sampleEnabled").checked) {
      const [ox, oy] = point(observation);
      ctx.strokeStyle = observation.locked ? "#d7ff52" : "#8bafc9";
      ctx.setLineDash([5 * dpr, 4 * dpr]); ctx.beginPath(); ctx.arc(ox, oy, Number($("observationRadius").value) / (2 * solver.extent) * size, 0, 2 * Math.PI); ctx.stroke(); ctx.setLineDash([]);
      ctx.fillStyle = "#d7ff52";
      for (const p of visible) { const [x, y] = point(p); ctx.beginPath(); ctx.arc(x, y, 2 * dpr, 0, 2 * Math.PI); ctx.fill(); }
    }
    if (selectedState?.ready && selectedState.alive) {
      const [x, y] = point(selectedState);
      ctx.strokeStyle = "#fff4d0"; ctx.lineWidth = 2 * dpr; ctx.beginPath(); ctx.arc(x, y, 7 * dpr, 0, 2 * Math.PI); ctx.stroke();
      ctx.fillStyle = "#fff4d0"; ctx.fillText("ID " + selectedState.id.slice(0, 15) + (selectedState.id.length > 15 ? "…" : ""), x + 10 * dpr, y - 10 * dpr);
    }
    ctx.restore();
  }
  function sample() {
    visible = []; pendingCount = 0;
    let budget = 32768;
    if (selected !== null) {
      selectedState = identities.resolve(selected, Math.min(4096, budget)); budget -= selectedState.work;
      $("identityDetail").textContent = selectedState.ready
        ? `ID ${selectedState.id}\n` + (selectedState.alive ? `x=${selectedState.x.toFixed(4)}, y=${selectedState.y.toFixed(4)} · vx=${selectedState.vx.toFixed(3)}, vy=${selectedState.vy.toFixed(3)}` : `Sorti du domaine (${selectedState.exit}).`) + `\nÉtape ${selectedState.tick} · identité conservée.`
        : `ID ${selected} · reconstruction ${selectedState.tick} / ${identities.history.length} étapes…`;
    }
    if ($("sampleEnabled").checked) {
      const r2 = Number($("observationRadius").value) ** 2;
      for (let j = 0; j < candidateIds.length; j++) {
        const id = candidateIds[(j + replayCursor) % candidateIds.length];
        const p = identities.resolve(id, Math.min(128, Math.max(0, budget)));
        budget -= p.work;
        if (!p.ready) pendingCount++;
        else if (p.alive && (p.x - observation.x) ** 2 + (p.y - observation.y) ** 2 <= r2) visible.push(p);
      }
      replayCursor = (replayCursor + 37) % candidateIds.length;
    }
  }
  function report() {
    const d = solver.diagnostics();
    $("massStat").value = (d.mass * 100).toFixed(2) + " %";
    $("timeStat").value = d.time.toFixed(3);
    $("spaceLoss").value = (d.spatialLoss * 100).toFixed(3) + " %";
    $("velocityLoss").value = (d.velocityLoss * 100).toFixed(3) + " %";
    $("angularStat").value = d.angularMomentum.toExponential(2);
    $("timing").value = ms.toFixed(1) + " ms";
    const error = Math.abs(d.mass + d.spatialLoss + d.velocityLoss - 1);
    $("status").textContent = `${solver.n}² positions × ${solver.nv}² vitesses · ${solver.size.toLocaleString("fr-FR")} cases. ` +
      (identities.full ? "Historique plein : expérience mise en pause. Réinitialisez pour recommencer. " : "") +
      (d.velocityLoss > 0.01 ? "Plus de 1 % de matière a quitté le domaine des vitesses : élargissez-le et réinitialisez. " : "") +
      `Bilan de masse, sorties incluses : erreur ${error.toExponential(1)}.`;
    $("sampleStatus").textContent = `N = ${identities.population.toString()} · ${visible.length} témoins dans la zone sur ${candidateIds.length} candidats. ` +
      `${identities.cache.size} états en cache · historique ${(identities.historyBytes / 1048576).toFixed(1)} / 32 Mio.` +
      (pendingCount ? ` Reconstruction de ${pendingCount} témoins en cours.` : "") +
      (observation.locked ? " Zone fixée." : " Zone mobile.");
  }
  for (const id of ["sigma", "fluctuation", "radius", "gravity", "softening", "speed", "observationRadius"]) {
    $(id).addEventListener("input", () => {
      const value = Number($(id).value);
      $(id + "Value").value = id === "fluctuation" ? Math.round(value * 100) + " %" : id === "speed" ? value.toFixed(1) + "×" : value.toFixed(2);
    });
  }
  for (const id of initialControls) $(id).addEventListener("change", () => { $("pending").hidden = false; $("sigma").disabled = $("initial").value === "cold"; });
  $("reset").addEventListener("click", reset);
  $("probes").addEventListener("change", () => { candidateIds = identities.candidates(Number($("probes").value)); });
  $("forget").addEventListener("click", () => { identities.forget(); });
  $("unlockObservation").addEventListener("click", () => { observation.locked = false; });
  function find(rebuild) {
    try {
      selected = identities.canonical($("identityInput").value);
      if (rebuild) identities.forgetIdentity(selected);
      sample(); draw();
    } catch (error) { $("identityDetail").textContent = error.message; selected = null; selectedState = null; }
  }
  $("findIdentity").addEventListener("click", () => find(false));
  $("rebuildIdentity").addEventListener("click", () => find(true));
  $("identityInput").addEventListener("keydown", event => { if (event.key === "Enter") find(false); });
  function worldPoint(event) {
    const rect = canvas.getBoundingClientRect();
    return { x: ((event.clientX - rect.left) * viewport.dpr - viewport.left) / viewport.size * 2 * solver.extent - solver.extent,
      y: solver.extent - ((event.clientY - rect.top) * viewport.dpr - viewport.top) / viewport.size * 2 * solver.extent };
  }
  canvas.addEventListener("pointermove", event => { if (!observation.locked) Object.assign(observation, worldPoint(event)); });
  canvas.addEventListener("pointerdown", event => {
    const point = worldPoint(event);
    const distance = 12 * viewport.dpr / viewport.size * 2 * solver.extent;
    let closest = null, best = distance ** 2;
    for (const p of visible) { const d = (p.x - point.x) ** 2 + (p.y - point.y) ** 2; if (d < best) { closest = p; best = d; } }
    if (closest) { $("identityInput").value = closest.id; find(false); }
    else Object.assign(observation, point, { locked: true });
  });
  $("pause").addEventListener("click", () => { running = !running; $("pause").textContent = running ? "Pause" : "Reprendre"; $("pause").setAttribute("aria-pressed", String(!running)); });
  document.addEventListener("visibilitychange", () => { last = performance.now(); });
  function frame(now) {
    const elapsed = Math.min(50, Math.max(0, now - last)); last = now;
    if (running && !document.hidden) {
      const start = performance.now();
      solver.advance(elapsed / 1000 * 0.25 * Number($("speed").value), { gravity: Number($("gravity").value), softening: Number($("softening").value) });
      if (identities.full) { running = false; $("pause").textContent = "Historique plein"; $("pause").disabled = true; $("pause").setAttribute("aria-pressed", "true"); }
      ms = ms * 0.85 + (performance.now() - start) * 0.15;
    }
    sample(); draw();
    if (now - reportAt > 350) { report(); reportAt = now; }
    requestAnimationFrame(frame);
  }
  reset(); requestAnimationFrame(frame);
})();
