/* Separate POC UI. A single committed render drives navigation and picking. */
(() => {
  "use strict";
  const $ = (id) => document.getElementById(id),
    canvas = $("scene"),
    camera = { x: 0, y: 0, zoom: 1 };
  let engine,
    population,
    renderer,
    diffuse,
    discovery,
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
    observationMemo = null,
    previousZoom = 1;
  const meter = new FrameRateMeter(),
    pan = { x: 0, y: 0 };
  let drag = null;
  let revision = 0;
  function invalidate() {
    dirty = true;
    revision++;
  }
  const number = (x) =>
    Number(x).toLocaleString("fr-FR", { maximumFractionDigits: 3 });
  function shape(view = camera) {
    const b = canvas.getBoundingClientRect();
    return renderContext(
      viewportGeometry(b.width, b.height, {
        pixelRatio: devicePixelRatio,
        maxDimension: engine?.device.limits.maxTextureDimension2D ?? 8192,
      }),
      { worldShortSpan: 4 / view.zoom, sigmaCss: 1, lengthPc: 10000 },
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
    discovery?.destroy();
    diffuse?.destroy();
    renderer?.destroy();
    engine?.destroy();
    discovery = diffuse = renderer = engine = null;
    observationMemo = null;
    lastLight = null;
    previousZoom = 1;
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
      diffuse = await GaussianDiffuseGPU.create({
        device: engine.device,
        population,
      });
      discovery = await AnchorDiscoveryGPU.create({
        device: engine.device,
        population,
      });
      ready = true;
      invalidate();
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
  async function draw() {
    if (tracked) {
      const p = population.resolve(tracked);
      camera.x = p.x + pan.x;
      camera.y = p.y + pan.y;
    }
    const drawRevision = revision,
      frameCamera = { ...camera },
      frameSelected = selected,
      frameTracked = tracked;
    const context = shape(frameCamera),
      view = {
        x: frameCamera.x,
        y: frameCamera.y,
        width: context.worldSpan,
        aspectX: context.aspectX,
        aspectY: context.aspectY,
      },
      projection = population.windowProjection(view);
    projection.cameraWidth = context.worldShortSpan;
    projection.zoom = frameCamera.zoom;
    projection.fovDegrees = 50;
    const budget = Math.max(
        0,
        Math.min(65536, Math.floor(Number($("sourceBudget").value) || 0)),
      ),
      onset = Math.max(1, Number($("sourceOnset").value) || 160),
      enabled = !$("densityOnly").checked && budget > 0,
      optics = POCStellarVisibility.options({
        population: population.population,
        context,
        exposureEV: Number($("exposure").value),
        onset,
      });
    const started = performance.now();
    const threshold = minimumLuminosity(projection, optics),
      filter = StellarPopulation.luminosityFilter(threshold);
    const signature = [
      budget,
      filter.key,
      enabled,
      optics.minimumPeak,
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
        enabled && filter.fraction > 0
          ? await (discovery || population).sample({
              ...view,
              budget,
              rankFilter: filter,
            })
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
    const { priorIds } = POCStellarVisibility.history(
      lastLight?.sources || [],
      {
        enabled,
        zoom: frameCamera.zoom,
        previousZoom,
      },
    );
    points = lastSample.points;
    if (frameSelected && enabled) {
      const p = population.resolve(frameSelected);
      if (
        Math.abs(p.x - frameCamera.x) < context.worldWidth / 2 &&
        Math.abs(p.y - frameCamera.y) < context.worldHeight / 2 &&
        !points.some((q) => q.id === p.id)
      )
        points = [...points, p];
    }
    points = enabled
      ? POCStellarVisibility.select(projection, points, {
          ...optics,
          budget,
          priorIds,
          selectedId: frameSelected,
        })
      : [];
    previousZoom = frameCamera.zoom;
    searchMs = performance.now() - started;
    renderer.resize(context.pixelWidth, context.pixelHeight);
    lastLight = renderer.render(null, projection, points, {
      ...optics,
      detailWeight: enabled ? 1 : 0,
      selected: frameSelected,
      complete: false,
      lengthPc: 10000,
      opticalShape: DEFAULT_OPTICAL_SHAPE,
      stellarFade: DEFAULT_STELLAR_FADE,
      diffuseLayer: diffuse,
    });
    renderer.present(canvas, {
      width: context.pixelWidth,
      height: context.pixelHeight,
      context,
    });
    // Picking uses this exact frame, never a live camera awaiting a new image.
    canvas._frame = {
      camera: { ...frameCamera },
      context,
      points: lastLight.drawn,
    };
    canvas.dataset.time = String(engine.time);
    canvas.dataset.sources = String(lastLight.sources.length);
    canvas.dataset.selectedPeak = String(
      lastLight.sources.find((p) => p.id === frameSelected)?.peak ?? 0,
    );
    canvas.dataset.followed = frameTracked || "";
    canvas.dataset.visited = String(lastSample.visited);
    canvas.dataset.approximate = String(lastSample.approximate);
    canvas.dataset.discoveryBackend = lastSample.backend || "none";
    canvas.dataset.discoveryThinning = String(lastSample.thinning || 1);
    canvas.dataset.searchMs = String(searchMs);
    canvas.dataset.backend = "GPU FMM anchors";
    canvas.dataset.diffuseBackend = diffuse.work.backend;
    $("time").textContent = number(engine.time) + " My";
    $("hud").textContent =
      number(engine.count) +
      " ancres · " +
      number(lastLight.sources.length) +
      " sources\n×" +
      number(frameCamera.zoom) +
      " · champ " +
      number(context.worldWidth * 10000) +
      " × " +
      number(context.worldHeight * 10000) +
      " pc";
    updateObject();
    dirty = revision !== drawRevision;
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
          invalidate();
        }
        if (dirty) await draw();
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
            "\nCatalogue : " +
            (lastSample?.backend || "désactivé") +
            (lastSample?.complete
              ? " · parcours terminé"
              : " · recherche approchée") +
            "\nFraction parcourue : 1/" +
            number(lastSample?.thinning || 1) +
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
      invalidate();
    };
  $("diagnostic").onclick = () => {
    tracked = null;
    camera.x = -0.208141793;
    camera.y = 0.120666986;
    camera.zoom = 80;
    invalidate();
    $("status").textContent =
      "Région physique du diagnostic. Atteindre 26,5 My pour comparer ; les identités de ce POC sont distinctes.";
  };
  function follow() {
    try {
      const p = population.resolve($("reference").value.trim());
      selected = tracked = p.id;
      pan.x = pan.y = 0;
      invalidate();
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
    invalidate();
  };
  canvas.addEventListener(
    "wheel",
    (e) => {
      e.preventDefault();
      camera.zoom = Math.max(
        1,
        Math.min(1e9, camera.zoom * Math.exp(-e.deltaY * 0.002)),
      );
      invalidate();
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
    invalidate();
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
      invalidate();
    }
  });
  canvas.addEventListener("pointercancel", () => {
    drag = null;
  });
  canvas.addEventListener("auxclick", (e) => {
    if (e.button === 1) e.preventDefault();
  });
  new ResizeObserver(() => {
    invalidate();
  }).observe($("viewport"));
  window.addEventListener("beforeunload", () => {
    discovery?.destroy();
    diffuse?.destroy();
    renderer?.destroy();
    engine?.destroy();
  });
  reset();
  requestAnimationFrame(loop);
})();
