/* Separate POC UI. A single committed render drives navigation and picking. */
(() => {
  "use strict";
  const $ = (id) => document.getElementById(id),
    canvas = $("scene"),
    camera = { x: 0, y: 0, zoom: 1 };
  let engine,
    population,
    renderer,
    ready = false,
    busy = false,
    running = false,
    target = null,
    dirty = true,
    generation = 0,
    selected = null,
    tracked = null,
    points = [],
    last = performance.now(),
    lastStats = 0,
    searchMs = 0,
    lastSample = null,
    lastLight = null,
    observationMemo = null;
  const meter = new FrameRateMeter(),
    pan = { x: 0, y: 0 };
  let drag = null;
  const number = (x) =>
    Number(x).toLocaleString("fr-FR", { maximumFractionDigits: 3 });
  function shape() {
    const b = canvas.getBoundingClientRect();
    return renderContext(
      viewportGeometry(b.width, b.height, {
        pixelRatio: devicePixelRatio,
        maxDimension: engine?.device.limits.maxTextureDimension2D ?? 8192,
      }),
      { worldShortSpan: 4 / camera.zoom, sigmaCss: 1, lengthPc: 10000 },
    );
  }
  function updateObject() {
    if (!selected || !population) {
      $("object").textContent = "Aucun objet sélectionné.";
      return;
    }
    const p = population.resolve(selected),
      props = stellarProperties(selected);
    $("object").textContent =
      stellarName(selected) +
      "\n" +
      selected +
      "\nPosition : " +
      number(p.x * 10000) +
      " ; " +
      number(p.y * 10000) +
      " pc\nVitesse : " +
      number(Math.hypot(p.vx, p.vy) * MilkyWay.velocityKmS) +
      " km/s\nMasse : " +
      number(props.mass) +
      " M☉\n" +
      (tracked === selected ? "Référentiel suivi" : "Sélection seule");
  }
  async function reset() {
    if (busy) return;
    const token = ++generation;
    busy = true;
    ready = false;
    running = false;
    target = null;
    $("play").disabled = $("advance").disabled = true;
    $("status").textContent = "Création des ancres et de la FMM GPU…";
    renderer?.destroy();
    engine?.destroy();
    renderer = engine = null;
    observationMemo = null;
    selected = tracked = null;
    points = [];
    camera.x = camera.y = 0;
    camera.zoom = 1;
    pan.x = pan.y = 0;
    try {
      const next = await AnchorEngine.create({
        count: Number($("anchorCount").value),
      });
      if (token !== generation) {
        next.destroy();
        return;
      }
      engine = next;
      population = new AnchorPopulation(engine);
      renderer = await GPUStellarLightRenderer.create({
        device: engine.device,
      });
      ready = true;
      dirty = true;
      $("status").textContent = "POC prêt · aucune grille cinétique.";
      canvas.dataset.ready = "true";
    } catch (error) {
      $("status").textContent = error.message;
      canvas.dataset.error = error.message;
    } finally {
      busy = false;
      $("play").disabled = $("advance").disabled = !ready;
      last = performance.now();
    }
  }
  function draw() {
    if (tracked) {
      const p = population.resolve(tracked);
      camera.x = p.x + pan.x;
      camera.y = p.y + pan.y;
    }
    const context = shape(),
      view = {
        x: camera.x,
        y: camera.y,
        width: context.worldSpan,
        aspectX: context.aspectX,
        aspectY: context.aspectY,
      },
      projection = population.project({ ...view, resolution: 96 });
    projection.cameraWidth = context.worldShortSpan;
    projection.zoom = camera.zoom;
    projection.fovDegrees = 50;
    const budget = Math.max(
        0,
        Math.min(65536, Math.floor(Number($("sourceBudget").value) || 0)),
      ),
      onset = Math.max(1, Number($("sourceOnset").value) || 160),
      detail = $("densityOnly").checked
        ? 0
        : stellarDetailLevel(projection.expected, budget, onset);
    const started = performance.now();
    let threshold =
        minimumLuminosity(projection, {
          population: population.population,
          gain: 1,
          exposureEV: Number($("exposure").value),
          context,
          lengthPc: 10000,
        }) / Math.max(detail, 1e-30),
      filter = StellarPopulation.luminosityFilter(threshold);
    for (
      let i = 0;
      i < 32 &&
      filter.fraction > 0 &&
      filter.fraction * projection.expected > budget * 0.6;
      i++
    ) {
      threshold *= 2;
      const next = StellarPopulation.luminosityFilter(threshold);
      if (!next.fraction) break;
      filter = next;
    }
    const signature = [
      budget,
      filter.key,
      detail > 0,
      context.worldSpan,
      context.aspectX,
      context.aspectY,
    ].join(":");
    const delta = observationMemo
      ? Math.abs(engine.time - observationMemo.time)
      : Infinity;
    const movement = observationMemo
      ? Math.hypot(view.x - observationMemo.x, view.y - observationMemo.y) +
        observationMemo.speed * delta
      : Infinity;
    const reuse =
      observationMemo &&
      signature === observationMemo.signature &&
      delta < 0.00001 &&
      movement < context.worldPerCssPixel * 0.25;
    if (reuse) {
      lastSample = {
        ...observationMemo.sample,
        points: observationMemo.ids
          .map((id) => population.resolve(id))
          .filter(
            (p) =>
              Math.abs(p.x - view.x) <= context.worldWidth / 2 &&
              Math.abs(p.y - view.y) <= context.worldHeight / 2,
          ),
        reused: true,
      };
    } else {
      lastSample =
        detail > 0 && budget
          ? population.sample({ ...view, budget, rankFilter: filter })
          : { points: [], visited: 0, approximate: false, shortfall: 0 };
      let speed = 0;
      for (const p of engine.anchors) {
        let b = 0;
        for (const v of p.B) b += Math.abs(v);
        speed = Math.max(speed, Math.hypot(p.vx, p.vy) + 8 * b);
      }
      observationMemo = {
        signature,
        time: engine.time,
        x: view.x,
        y: view.y,
        speed,
        sample: lastSample,
        ids: lastSample.points.map((p) => p.id),
      };
    }
    points = lastSample.points;
    if (selected && budget > 0 && !$("densityOnly").checked) {
      const p = population.resolve(selected);
      if (
        Math.abs(p.x - camera.x) < context.worldWidth / 2 &&
        Math.abs(p.y - camera.y) < context.worldHeight / 2 &&
        !points.some((q) => q.id === p.id)
      )
        points = [...points.slice(0, Math.max(0, budget - 1)), p];
    }
    searchMs = performance.now() - started;
    renderer.resize(context.pixelWidth, context.pixelHeight);
    lastLight = renderer.render(null, projection, points, {
      context,
      population: population.population,
      gain: 1,
      exposureEV: Number($("exposure").value),
      detailWeight: detail,
      selected,
      complete: false,
      lengthPc: 10000,
      opticalShape: DEFAULT_OPTICAL_SHAPE,
      stellarFade: DEFAULT_STELLAR_FADE,
    });
    renderer.present(canvas, {
      width: context.pixelWidth,
      height: context.pixelHeight,
      context,
    });
    // Picking uses this exact frame, never a live camera awaiting a new image.
    canvas._frame = { camera: { ...camera }, context, points: lastLight.drawn };
    canvas.dataset.time = String(engine.time);
    canvas.dataset.sources = String(lastLight.sources.length);
    canvas.dataset.followed = tracked || "";
    canvas.dataset.visited = String(lastSample.visited);
    canvas.dataset.approximate = String(lastSample.approximate);
    canvas.dataset.backend = "GPU FMM anchors";
    $("time").textContent = number(engine.time) + " My";
    $("hud").textContent =
      number(engine.count) +
      " ancres · " +
      number(lastLight.sources.length) +
      " sources\n×" +
      number(camera.zoom) +
      " · champ " +
      number(context.worldWidth * 10000) +
      " × " +
      number(context.worldHeight * 10000) +
      " pc";
    updateObject();
    dirty = false;
  }
  async function loop(now) {
    const elapsed = Math.min(0.1, Math.max(0, (now - last) / 1000));
    last = now;
    if (ready && !busy && !document.hidden) {
      busy = true;
      try {
        if (running || target !== null) {
          let duration =
            target !== null
              ? Math.min(0.8, Math.max(0, target - engine.time))
              : elapsed *
                Math.min(
                  lastLight?.sources.length ? 0.0001 : 1000,
                  Math.max(0, Number($("rate").value)) *
                    Number($("rateUnit").value),
                );
          if (duration > 0) await engine.advance(duration);
          if (target !== null && engine.time >= target - 1e-8) {
            target = null;
            $("advance").disabled = false;
            $("status").textContent = "Époque atteinte.";
          }
          dirty = true;
        }
        if (dirty) draw();
        meter.tick(performance.now());
        if (now - lastStats > 500) {
          lastStats = now;
          const d = engine.diagnostics();
          $("fps").textContent = meter.rate(now).toFixed(0) + " FPS";
          $("stats").textContent =
            "FMM + transport : " +
            number(engine.lastMs) +
            " ms/lot\nRecherche : " +
            number(searchMs) +
            " ms\nDécouverte : " +
            (lastSample?.reused
              ? "cache < 0,25 px / 10 ans"
              : "nouvelle recherche") +
            "\nNœuds visités : " +
            number(lastSample?.visited || 0) +
            "\nRecherche : " +
            (lastSample?.approximate
              ? "approximation bornée"
              : "résolution demandée atteinte") +
            "\nSources manquantes au budget : " +
            number(lastSample?.shortfall || 0) +
            "\nΔ moment cinétique : " +
            number(
              ((d.angularMomentum - engine.initialAngular) /
                Math.max(Math.abs(engine.initialAngular), 1e-20)) *
                100,
            ) +
            " %";
        }
      } catch (error) {
        running = false;
        target = null;
        $("status").textContent = error.stack;
        canvas.dataset.error = error.message;
      } finally {
        busy = false;
      }
    }
    requestAnimationFrame(loop);
  }
  $("play").onclick = () => {
    running = !running;
    target = null;
    $("advance").disabled = false;
    $("play").textContent = running ? "Pause" : "Démarrer";
  };
  $("advance").onclick = () => {
    const value = Number($("targetTime").value);
    if (!Number.isFinite(value) || value < engine.time) {
      $("status").textContent = "Choisir une époque future, ou recommencer.";
      return;
    }
    running = false;
    target = value;
    $("play").textContent = "Démarrer";
    $("advance").disabled = true;
    $("status").textContent = "Avancée FMM jusqu’à " + number(value) + " My…";
  };
  $("reset").onclick = () => {
    running = false;
    target = null;
    const retry = () => {
      if (busy) setTimeout(retry, 50);
      else reset();
    };
    retry();
  };
  for (const id of ["sourceBudget", "sourceOnset", "exposure", "densityOnly"])
    $(id).oninput = () => {
      dirty = true;
    };
  $("diagnostic").onclick = () => {
    tracked = null;
    camera.x = -0.208141793;
    camera.y = 0.120666986;
    camera.zoom = 80;
    dirty = true;
    $("status").textContent =
      "Région physique du diagnostic. Atteindre 26,5 My pour comparer ; les identités de ce POC sont distinctes.";
  };
  function follow() {
    try {
      const p = population.resolve($("reference").value.trim());
      selected = tracked = p.id;
      pan.x = pan.y = 0;
      dirty = true;
      $("status").textContent = "Objet suivi.";
    } catch (error) {
      $("status").textContent = error.message;
    }
  }
  $("follow").onclick = follow;
  $("reference").onkeydown = (e) => {
    if (e.key === "Enter") follow();
  };
  $("release").onclick = () => {
    tracked = null;
    dirty = true;
  };
  canvas.addEventListener(
    "wheel",
    (e) => {
      e.preventDefault();
      camera.zoom = Math.max(
        1,
        Math.min(1e9, camera.zoom * Math.exp(-e.deltaY * 0.002)),
      );
      dirty = true;
    },
    { passive: false },
  );
  canvas.addEventListener("pointerdown", (e) => {
    if (!ready || ![0, 1].includes(e.button)) return;
    e.preventDefault();
    drag = {
      x: e.clientX,
      y: e.clientY,
      camera: { ...camera },
      pan: { ...pan },
      button: e.button,
      moved: false,
      scale:
        canvas._frame?.context.worldPerCssPixel ?? shape().worldPerCssPixel,
    };
    canvas.setPointerCapture(e.pointerId);
  });
  canvas.addEventListener("pointermove", (e) => {
    if (!drag) return;
    const dx = e.clientX - drag.x,
      dy = e.clientY - drag.y;
    drag.moved ||= Math.abs(dx) + Math.abs(dy) > 3;
    if (!drag.moved) return;
    if (drag.button === 1 && tracked) {
      pan.x = drag.pan.x - dx * drag.scale;
      pan.y = drag.pan.y + dy * drag.scale;
    } else {
      tracked = null;
      camera.x = drag.camera.x - dx * drag.scale;
      camera.y = drag.camera.y + dy * drag.scale;
    }
    dirty = true;
  });
  canvas.addEventListener("pointerup", (e) => {
    if (!drag) return;
    const saved = drag;
    drag = null;
    if (saved.moved || !canvas._frame) return;
    const b = canvas.getBoundingClientRect(),
      f = canvas._frame,
      x =
        f.camera.x +
        (e.clientX - b.left - b.width / 2) * f.context.worldPerCssPixel,
      y =
        f.camera.y -
        (e.clientY - b.top - b.height / 2) * f.context.worldPerCssPixel;
    let best = 12 * f.context.worldPerCssPixel,
      pick = null;
    for (const p of f.points) {
      const d = Math.hypot(p.x - x, p.y - y);
      if (d < best) {
        best = d;
        pick = p;
      }
    }
    if (pick) {
      selected = pick.id;
      $("reference").value = selected;
      if (saved.button === 1) follow();
      dirty = true;
    }
  });
  canvas.addEventListener("pointercancel", () => {
    drag = null;
  });
  canvas.addEventListener("auxclick", (e) => {
    if (e.button === 1) e.preventDefault();
  });
  new ResizeObserver(() => {
    dirty = true;
  }).observe($("viewport"));
  window.addEventListener("beforeunload", () => {
    renderer?.destroy();
    engine?.destroy();
  });
  reset();
  requestAnimationFrame(loop);
})();
