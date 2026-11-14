/* Shared CPU FMM extracted from Orbit's existing implementation.
 * P2M / M2M / M2L / L2L / L2P; local order 0 or 1, bounded near field.
 * Each instance owns its tree cache. No DOM, rendering or integration.
 */
(function(root){
'use strict';
function createOrbitFMMField(){
 const NEAR_CELL_LIMIT=28, BOX_MIN=-1, BOX_SIZE=2;
 let bodies=[], treeCache=null, boundaryMode='open', torusActive=false;
 let gravityStrength=.48,gravityExponent=2,gravitySoftening2=.042**2,nearCellLimit=NEAR_CELL_LIMIT,minimumDepth=3;
 const els={boundary:{value:'open'},accuracy:{value:1}};
  function gravityKernel(q) {
    if (gravityExponent === 2) {
      const inverse = 1 / Math.sqrt(q);
      return inverse * inverse * inverse;
    }
    return Math.pow(q, -(gravityExponent + 1) * 0.5);
  }


  function createLevel(level, minX, minY, size) {
    const side = 1 << level;
    const cellSize = size / side;
    const cells = new Array(side * side);

    for (let y = 0; y < side; y++) {
      for (let x = 0; x < side; x++) {
        const index = y * side + x;
        cells[index] = {
          index,
          x,
          y,
          cx: minX + (x + 0.5) * cellSize,
          cy: minY + (y + 0.5) * cellSize,
          size: cellSize,
          mass: 0,
          mx: 0,
          my: 0,
          ax: 0,
          ay: 0,
          jxx: 0,
          jxy: 0,
          jyy: 0,
          bodies: null
        };
      }
    }
    return {
      side,
      cells,
      cellSize,
      seen: new Int32Array(cells.length)
    };
  }

  function prepareLevels(depth, minX, minY, size) {
    if (!treeCache || treeCache.depth !== depth) {
      const levels = [];
      for (let level = 0; level <= depth; level++) {
        levels.push(createLevel(level, minX, minY, size));
      }
      for (const cell of levels[depth].cells) cell.bodies = [];
      treeCache = { depth, levels };
    }

    for (let levelIndex = 0; levelIndex <= depth; levelIndex++) {
      const level = treeCache.levels[levelIndex];
      const cellSize = size / level.side;
      level.cellSize = cellSize;
      for (const cell of level.cells) {
        cell.cx = minX + (cell.x + 0.5) * cellSize;
        cell.cy = minY + (cell.y + 0.5) * cellSize;
        cell.size = cellSize;
        cell.mass = 0;
        cell.mx = 0;
        cell.my = 0;
        cell.ax = 0;
        cell.ay = 0;
        cell.jxx = 0;
        cell.jxy = 0;
        cell.jyy = 0;
        if (levelIndex === depth) cell.bodies.length = 0;
      }
    }
    return treeCache.levels;
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

  function periodicCellDistance(a, b, side) {
    const distance = Math.abs(a - b);
    return isTorus() ? Math.min(distance, side - distance) : distance;
  }

  // M2L d'ordre 1 : le monopôle source devient un champ local
  // (accélération + jacobienne) au centre de la cellule cible.
  function addSourceField(target, source, order) {
    if (source.mass <= 0) return 0;

    const rx = minimumImage(source.mx / source.mass - target.cx);
    const ry = minimumImage(source.my / source.mass - target.cy);
    const q = rx * rx + ry * ry + gravitySoftening2;
    const kernel = gravityKernel(q);
    const strength = gravityStrength * source.mass;
    target.ax += strength * rx * kernel;
    target.ay += strength * ry * kernel;

    if (order > 0) {
      const derivative = kernel / q;
      const coefficient = gravityExponent + 1;
      target.jxx +=
        strength * (-kernel + coefficient * rx * rx * derivative);
      target.jxy +=
        strength * (coefficient * rx * ry * derivative);
      target.jyy +=
        strength * (-kernel + coefficient * ry * ry * derivative);
    }
    return 1;
  }

  function directField(body, source) {
    const rx = minimumImage(source.x - body.x);
    const ry = minimumImage(source.y - body.y);
    const q = rx * rx + ry * ry + gravitySoftening2;
    const f = gravityStrength * source.mass * gravityKernel(q);
    body.ax += rx * f;
    body.ay += ry * f;
  }

  // Garde-fou linéaire : une feuille pathologiquement dense est condensée.
  // Ainsi, le voisinage direct reste borné même si toutes les particules
  // se concentrent dans la même cellule.
  function aggregateNearField(body, cell) {
    const ownMass = cell.index === body.cellIndex ? body.mass : 0;
    const mass = cell.mass - ownMass;
    if (mass <= 1e-12) return;

    const sx = (cell.mx - body.x * ownMass) / mass;
    const sy = (cell.my - body.y * ownMass) / mass;
    const rx = minimumImage(sx - body.x);
    const ry = minimumImage(sy - body.y);
    const q = rx * rx + ry * ry + gravitySoftening2;
    const f = gravityStrength * mass * gravityKernel(q);
    body.ax += rx * f;
    body.ay += ry * f;
  }


  function computeFMM() {
    const startedAt = performance.now();
    const n = bodies.length;
    boundaryMode = els.boundary.value;
    torusActive = boundaryMode === "torus";
    // 4^L est proportionnel à N : le nombre total de cellules est O(N).
    const depth = Math.max(
      minimumDepth,
      Math.min(6, Math.ceil(Math.log(Math.max(1, n / 4)) / Math.log(4)))
    );

    let minX = BOX_MIN;
    let minY = BOX_MIN;
    let span = BOX_SIZE;
    if (els.boundary.value === "open") {
      let maxX = -Infinity;
      let maxY = -Infinity;
      minX = Infinity;
      minY = Infinity;
      for (const body of bodies) {
        if (body.x < minX) minX = body.x;
        if (body.x > maxX) maxX = body.x;
        if (body.y < minY) minY = body.y;
        if (body.y > maxY) maxY = body.y;
      }
      span = Math.max(maxX - minX, maxY - minY, 1.8) * 1.08;
      const centerX = (minX + maxX) * 0.5;
      const centerY = (minY + maxY) * 0.5;
      minX = centerX - span * 0.5;
      minY = centerY - span * 0.5;
    }

    const levels = prepareLevels(depth, minX, minY, span);

    // P2M : chaque particule est insérée une seule fois.
    const leaf = levels[depth];
    for (const b of bodies) {
      const ix = Math.max(
        0,
        Math.min(leaf.side - 1, Math.floor((b.x - minX) / leaf.cellSize))
      );
      const iy = Math.max(
        0,
        Math.min(leaf.side - 1, Math.floor((b.y - minY) / leaf.cellSize))
      );
      const index = iy * leaf.side + ix;
      const cell = leaf.cells[index];
      cell.mass += b.mass;
      cell.mx += b.x * b.mass;
      cell.my += b.y * b.mass;
      cell.bodies.push(b);
      b.cellIndex = index;
      b.ax = 0;
      b.ay = 0;
    }

    // M2M : remontée des moments, quatre enfants par cellule.
    for (let level = depth - 1; level >= 0; level--) {
      const current = levels[level];
      const child = levels[level + 1];
      for (const cell of current.cells) {
        for (let oy = 0; oy < 2; oy++) {
          for (let ox = 0; ox < 2; ox++) {
            const c =
              child.cells[
                (cell.y * 2 + oy) * child.side + cell.x * 2 + ox
              ];
            cell.mass += c.mass;
            cell.mx += c.mx;
            cell.my += c.my;
          }
        }
      }
    }

    // M2L : chaque cellule a une liste d'interaction de taille bornée (≤ 27).
    let ops = n;
    const order = Number(els.accuracy.value);
    for (let levelIndex = 2; levelIndex <= depth; levelIndex++) {
      const level = levels[levelIndex];
      const parent = levels[levelIndex - 1];
      if (isTorus()) level.seen.fill(0);

      for (const target of level.cells) {
        if (target.mass <= 0) continue;
        const stamp = target.index + 1;
        const px = target.x >> 1;
        const py = target.y >> 1;
        for (let noy = -1; noy <= 1; noy++) {
          for (let nox = -1; nox <= 1; nox++) {
            let npx = px + nox;
            let npy = py + noy;
            if (isTorus()) {
              npx = (npx + parent.side) % parent.side;
              npy = (npy + parent.side) % parent.side;
            } else if (
              npx < 0 ||
              npy < 0 ||
              npx >= parent.side ||
              npy >= parent.side
            ) {
              continue;
            }

            for (let cy = 0; cy < 2; cy++) {
              for (let cx = 0; cx < 2; cx++) {
                const sx = npx * 2 + cx;
                const sy = npy * 2 + cy;
                const sourceIndex = sy * level.side + sx;
                if (isTorus()) {
                  if (level.seen[sourceIndex] === stamp) continue;
                  level.seen[sourceIndex] = stamp;
                }
                if (
                  periodicCellDistance(sx, target.x, level.side) <= 1 &&
                  periodicCellDistance(sy, target.y, level.side) <= 1
                ) continue;
                ops += addSourceField(
                  target,
                  level.cells[sourceIndex],
                  order
                );
              }
            }
          }
        }
      }
    }

    // L2L : le champ local du parent est traduit vers chaque enfant.
    for (let levelIndex = 1; levelIndex <= depth; levelIndex++) {
      const level = levels[levelIndex];
      const parent = levels[levelIndex - 1];
      for (const cell of level.cells) {
        if (cell.mass <= 0) continue;
        const p =
          parent.cells[(cell.y >> 1) * parent.side + (cell.x >> 1)];
        const dx = cell.cx - p.cx;
        const dy = cell.cy - p.cy;
        cell.ax += p.ax + p.jxx * dx + p.jxy * dy;
        cell.ay += p.ay + p.jxy * dx + p.jyy * dy;
        if (order > 0) {
          cell.jxx += p.jxx;
          cell.jxy += p.jxy;
          cell.jyy += p.jyy;
        }
      }
    }

    // L2P + P2P : champ local puis neuf cellules voisines au maximum.
    for (const cell of leaf.cells) {
      if (!cell.bodies.length) continue;

      for (const body of cell.bodies) {
        const dx = body.x - cell.cx;
        const dy = body.y - cell.cy;
        body.ax +=
          cell.ax + (order ? cell.jxx * dx + cell.jxy * dy : 0);
        body.ay +=
          cell.ay + (order ? cell.jxy * dx + cell.jyy * dy : 0);

        for (let oy = -1; oy <= 1; oy++) {
          for (let ox = -1; ox <= 1; ox++) {
            let nx = cell.x + ox;
            let ny = cell.y + oy;
            if (isTorus()) {
              nx = (nx + leaf.side) % leaf.side;
              ny = (ny + leaf.side) % leaf.side;
            } else if (
              nx < 0 ||
              ny < 0 ||
              nx >= leaf.side ||
              ny >= leaf.side
            ) {
              continue;
            }

            const near = leaf.cells[ny * leaf.side + nx];
            if (near.bodies.length <= nearCellLimit) {
              for (const source of near.bodies) {
                if (source !== body) {
                  directField(body, source);
                  ops++;
                }
              }
            } else {
              aggregateNearField(body, near);
              ops++;
            }
          }
        }
      }
    }


 return { levels, depth, operations:ops };
 }
 return {compute(particles,options={}){
 bodies=particles;els.boundary.value=options.boundary??'open';els.accuracy.value=options.order??1;
 gravityStrength=options.gravity??.48;gravityExponent=options.exponent??2;gravitySoftening2=(options.softening??.042)**2;
 nearCellLimit=options.nearLimit??NEAR_CELL_LIMIT;
 minimumDepth=options.minDepth??3;
 return computeFMM();
 }};
}
root.createOrbitFMMField=createOrbitFMMField;
if(typeof module!=='undefined')module.exports=createOrbitFMMField;
})(typeof window!=='undefined'?window:globalThis);
