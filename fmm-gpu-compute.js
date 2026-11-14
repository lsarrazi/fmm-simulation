(() => {
  "use strict";

  const WORKGROUP_SIZE = 256;
  const MAX_DISPATCH_WORKGROUPS = 65535;
  const DIRECT_LIMIT = 28;
  const MASS_SCALE = 268435456;
  const PARTICLE_STRIDE = 32;
  const CELL_STRIDE = 48;
  const TOPOLOGY_STRIDE = 32;
  const INTERACTIONS_PER_CELL = 27;

  function levelOffset(level) {
    return ((1 << (2 * level)) - 1) / 3;
  }

  function depthForCount(count) {
    return Math.max(
      3,
      Math.min(6, Math.ceil(Math.log(Math.max(1, count / 4)) / Math.log(4)))
    );
  }

  function dispatchParticleWorkgroups(pass, count) {
    const workgroups = Math.ceil(count / WORKGROUP_SIZE);
    const x = Math.min(workgroups, MAX_DISPATCH_WORKGROUPS);
    const y = Math.ceil(workgroups / MAX_DISPATCH_WORKGROUPS);
    pass.dispatchWorkgroups(x, y);
  }

  function buildTopology(depth) {
    const totalCells = levelOffset(depth + 1);
    const leafOffset = levelOffset(depth);
    const leafSide = 1 << depth;
    const topologyData = new ArrayBuffer(totalCells * TOPOLOGY_STRIDE);
    const topologyFloats = new Float32Array(topologyData);
    const topologyInts = new Uint32Array(topologyData);
    const interactionData = new Uint32Array(
      totalCells * INTERACTIONS_PER_CELL * 2
    );

    for (let level = 0; level <= depth; level++) {
      const side = 1 << level;
      const size = 2 / side;
      const offset = levelOffset(level);
      for (let y = 0; y < side; y++) {
        for (let x = 0; x < side; x++) {
          const local = y * side + x;
          const flat = offset + local;
          const word = flat * 8;
          topologyFloats[word] = -1 + (x + 0.5) * size;
          topologyFloats[word + 1] = -1 + (y + 0.5) * size;
          topologyFloats[word + 2] = size;
          topologyInts[word + 3] =
            level === 0
              ? 0
              : levelOffset(level - 1) + (y >> 1) * (side >> 1) + (x >> 1);
        }
      }
    }

    function writeInteractions(torus) {
      const interactionBase = torus
        ? totalCells * INTERACTIONS_PER_CELL
        : 0;
      for (let level = 2; level <= depth; level++) {
        const side = 1 << level;
        const parentSide = side >> 1;
        const offset = levelOffset(level);
        for (let ty = 0; ty < side; ty++) {
          for (let tx = 0; tx < side; tx++) {
            const targetFlat = offset + ty * side + tx;
            const seen = new Set();
            const sources = [];
            const px = tx >> 1;
            const py = ty >> 1;
            for (let oy = -1; oy <= 1; oy++) {
              for (let ox = -1; ox <= 1; ox++) {
                let npx = px + ox;
                let npy = py + oy;
                if (torus) {
                  npx = (npx + parentSide) % parentSide;
                  npy = (npy + parentSide) % parentSide;
                } else if (
                  npx < 0 ||
                  npy < 0 ||
                  npx >= parentSide ||
                  npy >= parentSide
                ) {
                  continue;
                }
                for (let cy = 0; cy < 2; cy++) {
                  for (let cx = 0; cx < 2; cx++) {
                    const sx = npx * 2 + cx;
                    const sy = npy * 2 + cy;
                    const dx = Math.abs(sx - tx);
                    const dy = Math.abs(sy - ty);
                    const periodicDx = torus ? Math.min(dx, side - dx) : dx;
                    const periodicDy = torus ? Math.min(dy, side - dy) : dy;
                    if (periodicDx <= 1 && periodicDy <= 1) continue;
                    const sourceFlat = offset + sy * side + sx;
                    if (seen.has(sourceFlat)) continue;
                    seen.add(sourceFlat);
                    sources.push(sourceFlat);
                  }
                }
              }
            }
            if (sources.length > INTERACTIONS_PER_CELL) {
              throw new Error("liste M2L supérieure à la borne prévue");
            }
            const topologyWord = targetFlat * 8;
            topologyInts[topologyWord + (torus ? 5 : 4)] = sources.length;
            const destination =
              interactionBase + targetFlat * INTERACTIONS_PER_CELL;
            for (let index = 0; index < sources.length; index++) {
              interactionData[destination + index] = sources[index];
            }
          }
        }
      }
    }

    writeInteractions(false);
    writeInteractions(true);
    return {
      depth,
      totalCells,
      leafOffset,
      leafSide,
      leafCount: leafSide * leafSide,
      topologyData,
      interactionData
    };
  }

  function packParticles(particles) {
    const data = new Float32Array(particles.length * 8);
    for (let index = 0; index < particles.length; index++) {
      const body = particles[index];
      const offset = index * 8;
      data[offset] = body.x;
      data[offset + 1] = body.y;
      data[offset + 2] = body.vx;
      data[offset + 3] = body.vy;
      data[offset + 4] = 0;
      data[offset + 5] = 0;
      data[offset + 6] = body.mass;
      data[offset + 7] = body.color;
    }
    return data;
  }

  const computeWGSL = `
    const WORKGROUP_SIZE: u32 = 256u;
    const DIRECT_LIMIT: u32 = 28u;
    const MASS_SCALE: f32 = 268435456.0;
    const BOX_SIZE: f32 = 2.0;

    struct Config {
      bodyCount: u32,
      leafSide: u32,
      leafOffset: u32,
      totalCells: u32,
      depth: u32,
      boundary: u32,
      order: u32,
      interactionBase: u32,
      dt: f32,
      softening2: f32,
      gravitation: f32,
      damping: f32,
      impulseX: f32,
      impulseY: f32,
      impulseForce: f32,
      expansionScale: f32,
      worldScaleX: f32,
      worldScaleY: f32,
      pointScaleX: f32,
      pointScaleY: f32,
      gravityPower: f32,
      viewBrightness: f32,
      pad1: f32,
      pad2: f32,
    }

    struct Particle {
      position: vec2<f32>,
      velocity: vec2<f32>,
      acceleration: vec2<f32>,
      mass: f32,
      color: f32,
    }

    struct LeafAccum {
      mass: atomic<i32>,
      mx: atomic<i32>,
      my: atomic<i32>,
      count: atomic<u32>,
    }

    struct Cell {
      mass: f32,
      mx: f32,
      my: f32,
      pad0: f32,
      ax: f32,
      ay: f32,
      jxx: f32,
      jxy: f32,
      jyy: f32,
      pad1: f32,
      pad2: f32,
      pad3: f32,
    }

    struct Topology {
      center: vec2<f32>,
      size: f32,
      parent: u32,
      normalCount: u32,
      torusCount: u32,
      pad0: u32,
      pad1: u32,
    }

    @group(0) @binding(0) var<uniform> config: Config;
    @group(0) @binding(1) var<storage, read> particlesIn: array<Particle>;
    @group(0) @binding(2) var<storage, read_write> particlesOut: array<Particle>;
    @group(0) @binding(3) var<storage, read_write> leafAccum: array<LeafAccum>;
    @group(0) @binding(4) var<storage, read_write> leafSlots: array<atomic<u32>>;
    @group(0) @binding(5) var<storage, read_write> particleCell: array<u32>;
    @group(0) @binding(6) var<storage, read_write> cells: array<Cell>;
    @group(0) @binding(7) var<storage, read> topology: array<Topology>;
    @group(0) @binding(8) var<storage, read> interactions: array<u32>;

    fn levelOffset(level: u32) -> u32 {
      return ((1u << (2u * level)) - 1u) / 3u;
    }

    fn particleIndex(globalId: vec3<u32>) -> u32 {
      return globalId.x + globalId.y * 16776960u;
    }

    fn minimumImage(delta: f32) -> f32 {
      if (config.boundary != 2u) { return delta; }
      if (delta > 1.0) { return delta - BOX_SIZE; }
      if (delta < -1.0) { return delta + BOX_SIZE; }
      return delta;
    }

    fn accelerationFrom(
      targetPosition: vec2<f32>,
      source: vec2<f32>,
      sourceMass: f32
    ) -> vec2<f32> {
      let displacement = vec2<f32>(
        minimumImage(source.x - targetPosition.x),
        minimumImage(source.y - targetPosition.y)
      );
      let q = dot(displacement, displacement) + config.softening2;
      var kernel: f32;
      if (abs(config.gravityPower - 2.0) < 0.0001) {
        let inverse = inverseSqrt(q);
        kernel = inverse * inverse * inverse;
      } else {
        kernel = pow(q, -(config.gravityPower + 1.0) * 0.5);
      }
      return displacement * (
        config.gravitation * sourceMass * kernel
      );
    }

    @compute @workgroup_size(256)
    fn clear(
      @builtin(global_invocation_id) globalId: vec3<u32>
    ) {
      let index = globalId.x;
      if (index < config.totalCells) {
        var emptyCell: Cell;
        cells[index] = emptyCell;
      }
      let leafCount = config.leafSide * config.leafSide;
      if (index < leafCount) {
        atomicStore(&leafAccum[index].mass, 0);
        atomicStore(&leafAccum[index].mx, 0);
        atomicStore(&leafAccum[index].my, 0);
        atomicStore(&leafAccum[index].count, 0u);
        if(config.pad2>0.0){atomicStore(&leafSlots[leafCount*DIRECT_LIMIT+index],0xffffffffu);}
      }
    }

    @compute @workgroup_size(256)
    fn binParticles(
      @builtin(global_invocation_id) globalId: vec3<u32>
    ) {
      let index = particleIndex(globalId);
      if (index >= config.bodyCount) { return; }
      let particle = particlesIn[index];
      let normalized = clamp(
        (particle.position + vec2<f32>(1.0)) * 0.5,
        vec2<f32>(0.0),
        vec2<f32>(0.999999)
      );
      let coordinate = min(
        vec2<u32>(normalized * f32(config.leafSide)),
        vec2<u32>(config.leafSide - 1u)
      );
      let cell = coordinate.y * config.leafSide + coordinate.x;
      particleCell[index] = cell;
      if(config.pad2>0.0){particleCell[config.bodyCount+index]=atomicExchange(&leafSlots[config.leafSide*config.leafSide*DIRECT_LIMIT+cell],index);}
      let slot = atomicAdd(&leafAccum[cell].count, 1u);
      if (slot < DIRECT_LIMIT) {
        atomicStore(&leafSlots[cell * DIRECT_LIMIT + slot],index);
      }

      let scaledMass = i32(round(max(0.0, particle.mass) * MASS_SCALE));
      let scaledMx = i32(round(
        clamp(particle.mass * particle.position.x, -2.0, 2.0) * MASS_SCALE
      ));
      let scaledMy = i32(round(
        clamp(particle.mass * particle.position.y, -2.0, 2.0) * MASS_SCALE
      ));
      atomicAdd(&leafAccum[cell].mass, scaledMass);
      atomicAdd(&leafAccum[cell].mx, scaledMx);
      atomicAdd(&leafAccum[cell].my, scaledMy);
    }

    @compute @workgroup_size(256)
    fn p2m(
      @builtin(global_invocation_id) globalId: vec3<u32>
    ) {
      let leaf = globalId.x;
      let leafCount = config.leafSide * config.leafSide;
      if (leaf >= leafCount) { return; }
      let targetCell = config.leafOffset + leaf;
      var result: Cell;
      result.mass = f32(atomicLoad(&leafAccum[leaf].mass)) / MASS_SCALE;
      result.mx = f32(atomicLoad(&leafAccum[leaf].mx)) / MASS_SCALE;
      result.my = f32(atomicLoad(&leafAccum[leaf].my)) / MASS_SCALE;
      cells[targetCell] = result;
    }

    @compute @workgroup_size(256)
    fn m2m(
      @builtin(local_invocation_index) localIndex: u32
    ) {
      var level = config.depth;
      loop {
        if (level == 0u) { break; }
        let parentLevel = level - 1u;
        let parentCount = 1u << (2u * parentLevel);
        let parentOffset = levelOffset(parentLevel);
        let childOffset = levelOffset(level);
        var index = localIndex;
        loop {
          if (index >= parentCount) { break; }
          let parentSide = 1u << parentLevel;
          let childSide = parentSide * 2u;
          let px = index % parentSide;
          let py = index / parentSide;
          var result: Cell;
          for (var oy: u32 = 0u; oy < 2u; oy = oy + 1u) {
            for (var ox: u32 = 0u; ox < 2u; ox = ox + 1u) {
              let child =
                childOffset + (py * 2u + oy) * childSide + px * 2u + ox;
              result.mass = result.mass + cells[child].mass;
              result.mx = result.mx + cells[child].mx;
              result.my = result.my + cells[child].my;
            }
          }
          cells[parentOffset + index] = result;
          index = index + WORKGROUP_SIZE;
        }
        storageBarrier();
        workgroupBarrier();
        level = parentLevel;
      }
    }

    @compute @workgroup_size(256)
    fn m2l(
      @builtin(global_invocation_id) globalId: vec3<u32>
    ) {
      let targetIndex = globalId.x;
      if (targetIndex >= config.totalCells) { return; }
      var targetCell = cells[targetIndex];
      if (targetCell.mass <= 0.0) { return; }
      let targetTopology = topology[targetIndex];
      let interactionCount = select(
        targetTopology.normalCount,
        targetTopology.torusCount,
        config.boundary == 2u
      );
      var ax = 0.0;
      var ay = 0.0;
      var jxx = 0.0;
      var jxy = 0.0;
      var jyy = 0.0;
      for (
        var interaction: u32 = 0u;
        interaction < interactionCount;
        interaction = interaction + 1u
      ) {
        let sourceIndex =
          interactions[
            config.interactionBase +
            targetIndex * 27u +
            interaction
          ];
        let source = cells[sourceIndex];
        if (source.mass <= 0.0) { continue; }
        let sourceCenter = vec2<f32>(
          source.mx / source.mass,
          source.my / source.mass
        );
        let displacement = vec2<f32>(
          minimumImage(sourceCenter.x - targetTopology.center.x),
          minimumImage(sourceCenter.y - targetTopology.center.y)
        );
        let q = dot(displacement, displacement) + config.softening2;
        var kernel: f32;
        if (abs(config.gravityPower - 2.0) < 0.0001) {
          let inverse = inverseSqrt(q);
          kernel = inverse * inverse * inverse;
        } else {
          kernel = pow(q, -(config.gravityPower + 1.0) * 0.5);
        }
        let strength = config.gravitation * source.mass;
        ax = ax + strength * displacement.x * kernel;
        ay = ay + strength * displacement.y * kernel;
        if (config.order > 0u) {
          let derivative = kernel / q;
          let coefficient = config.gravityPower + 1.0;
          jxx = jxx + strength * (
            -kernel +
            coefficient * displacement.x * displacement.x * derivative
          );
          jxy = jxy + strength * (
            coefficient * displacement.x * displacement.y * derivative
          );
          jyy = jyy + strength * (
            -kernel +
            coefficient * displacement.y * displacement.y * derivative
          );
        }
      }
      targetCell.ax = ax;
      targetCell.ay = ay;
      targetCell.jxx = jxx;
      targetCell.jxy = jxy;
      targetCell.jyy = jyy;
      cells[targetIndex] = targetCell;
    }

    @compute @workgroup_size(256)
    fn l2l(
      @builtin(local_invocation_index) localIndex: u32
    ) {
      var level = 1u;
      loop {
        if (level > config.depth) { break; }
        let offset = levelOffset(level);
        let count = 1u << (2u * level);
        var index = localIndex;
        loop {
          if (index >= count) { break; }
          let cellIndex = offset + index;
          var cell = cells[cellIndex];
          if (cell.mass > 0.0) {
            let parentIndex = topology[cellIndex].parent;
            let parent = cells[parentIndex];
            let delta = topology[cellIndex].center - topology[parentIndex].center;
            cell.ax = cell.ax + parent.ax + parent.jxx * delta.x + parent.jxy * delta.y;
            cell.ay = cell.ay + parent.ay + parent.jxy * delta.x + parent.jyy * delta.y;
            if (config.order > 0u) {
              cell.jxx = cell.jxx + parent.jxx;
              cell.jxy = cell.jxy + parent.jxy;
              cell.jyy = cell.jyy + parent.jyy;
            }
            cells[cellIndex] = cell;
          }
          index = index + WORKGROUP_SIZE;
        }
        storageBarrier();
        workgroupBarrier();
        level = level + 1u;
      }
    }

    fn wrappedCell(value: i32, side: i32) -> i32 {
      return ((value % side) + side) % side;
    }

    @compute @workgroup_size(256)
    fn evaluateAndIntegrate(
      @builtin(global_invocation_id) globalId: vec3<u32>
    ) {
      let index = particleIndex(globalId);
      if (index >= config.bodyCount) { return; }
      let input = particlesIn[index];
      let leaf = particleCell[index];
      let leafCell = cells[config.leafOffset + leaf];
      let leafCenter = topology[config.leafOffset + leaf].center;
      let localDelta = input.position - leafCenter;
      var acceleration = vec2<f32>(
        leafCell.ax + select(0.0, leafCell.jxx * localDelta.x + leafCell.jxy * localDelta.y, config.order > 0u),
        leafCell.ay + select(0.0, leafCell.jxy * localDelta.x + leafCell.jyy * localDelta.y, config.order > 0u)
      );

      let cellX = i32(leaf % config.leafSide);
      let cellY = i32(leaf / config.leafSide);
      let side = i32(config.leafSide);
      for (var oy: i32 = -1; oy <= 1; oy = oy + 1) {
        for (var ox: i32 = -1; ox <= 1; ox = ox + 1) {
          var nx = cellX + ox;
          var ny = cellY + oy;
          if (config.boundary == 2u) {
            nx = wrappedCell(nx, side);
            ny = wrappedCell(ny, side);
          } else if (nx < 0 || ny < 0 || nx >= side || ny >= side) {
            continue;
          }
          let near = u32(ny) * config.leafSide + u32(nx);
          let count = atomicLoad(&leafAccum[near].count);
          if (count <= DIRECT_LIMIT) {
            for (var slot: u32 = 0u; slot < count; slot = slot + 1u) {
              let sourceIndex = atomicLoad(&leafSlots[near * DIRECT_LIMIT + slot]);
              if (sourceIndex == index) { continue; }
              let source = particlesIn[sourceIndex];
              acceleration = acceleration + accelerationFrom(
                input.position,
                source.position,
                source.mass
              );
            }
          } else if(config.pad2>0.0){
            // Exact overflow traversal in the planetary patch, linear in local occupancy.
            var sourceIndex=atomicLoad(&leafSlots[config.leafSide*config.leafSide*DIRECT_LIMIT+near]);
            loop{if(sourceIndex==0xffffffffu){break;}if(sourceIndex!=index){let source=particlesIn[sourceIndex];acceleration+=accelerationFrom(input.position,source.position,source.mass);}sourceIndex=particleCell[config.bodyCount+sourceIndex];}
          } else if (config.pad1 > 0.0) {
            // Local physical refinement cannot replace dense neighboring stars
            // by a single aggregate. Reuse the same FMM, with exact overflow P2P.
            for (var sourceIndex: u32 = 0u; sourceIndex < config.bodyCount; sourceIndex++) {
              if (sourceIndex == index || particleCell[sourceIndex] != near) { continue; }
              let source = particlesIn[sourceIndex];
              acceleration += accelerationFrom(input.position, source.position, source.mass);
            }
          } else {
            let sourceCell = cells[config.leafOffset + near];
            let ownMass = select(0.0, input.mass, near == leaf);
            let sourceMass = sourceCell.mass - ownMass;
            if (sourceMass > 0.0000001) {
              let sourceCenter = vec2<f32>(
                (sourceCell.mx - input.position.x * ownMass) / sourceMass,
                (sourceCell.my - input.position.y * ownMass) / sourceMass
              );
              acceleration = acceleration + accelerationFrom(
                input.position,
                sourceCenter,
                sourceMass
              );
            }
          }
        }
      }

      // En frontière ouverte, la topologie GPU couvre [-1, 1]². Une
      // particule sortie de cette racine ne doit pas extrapoler le champ
      // local de la cellule du bord : cette extrapolation peut diverger.
      // Hors racine, le monopôle global est une approximation lointaine
      // stable, bornée et de coût constant par particule.
      if (
        config.boundary == 0u &&
        (abs(input.position.x) > 1.0 || abs(input.position.y) > 1.0)
      ) {
        let root = cells[0];
        let sourceMass = root.mass - input.mass;
        acceleration = vec2<f32>(0.0);
        if (sourceMass > 0.0000001) {
          let sourceCenter = vec2<f32>(
            (root.mx - input.position.x * input.mass) / sourceMass,
            (root.my - input.position.y * input.mass) / sourceMass
          );
          acceleration = accelerationFrom(
            input.position,
            sourceCenter,
            sourceMass
          );
        }
      }

      // An explicitly split dominant star is evaluated exactly in the shared
      // force kernel; all remaining interactions still use the same FMM.
      if (config.pad2 > 0.0 && index > 0u) {
        acceleration += accelerationFrom(input.position, particlesIn[0].position, config.pad2);
      }

      if (config.impulseForce > 0.0) {
        let impulse = vec2<f32>(config.impulseX, config.impulseY);
        let displacement = vec2<f32>(
          minimumImage(impulse.x - input.position.x),
          minimumImage(impulse.y - input.position.y)
        );
        let q = dot(displacement, displacement) + 0.012;
        acceleration = acceleration + displacement * (
          config.impulseForce / (q * sqrt(q))
        );
      }

      var output = input;
      output.acceleration = acceleration;
      output.velocity = (input.velocity + acceleration * config.dt) * config.damping;
      output.position =
        (input.position + output.velocity * config.dt) * config.expansionScale;

      if (config.boundary == 2u) {
        output.position = vec2<f32>(
          -1.0 + fract((output.position.x + 1.0) * 0.5) * 2.0,
          -1.0 + fract((output.position.y + 1.0) * 0.5) * 2.0
        );
      } else if (config.boundary == 1u) {
        if (output.position.x < -1.0) {
          output.position.x = -2.0 - output.position.x;
          output.velocity.x = -output.velocity.x;
        } else if (output.position.x > 1.0) {
          output.position.x = 2.0 - output.position.x;
          output.velocity.x = -output.velocity.x;
        }
        if (output.position.y < -1.0) {
          output.position.y = -2.0 - output.position.y;
          output.velocity.y = -output.velocity.y;
        } else if (output.position.y > 1.0) {
          output.position.y = 2.0 - output.position.y;
          output.velocity.y = -output.velocity.y;
        }
      }
      particlesOut[index] = output;
    }
    // Heliocentric kick-drift-kick. The shared FMM writes the new absolute
    // accelerations between these two passes; subtract the star's reaction.
    @compute @workgroup_size(256) fn orbitDrift(@builtin(global_invocation_id) id:vec3<u32>){
      let i=particleIndex(id);if(i>=config.bodyCount){return;}var p=particlesIn[i];
      if(i>0u){p.velocity+=0.5*config.worldScaleX*(p.acceleration-particlesIn[0].acceleration);p.position+=config.worldScaleX*p.velocity;}
      particlesOut[i]=p;
    }
    @compute @workgroup_size(256) fn orbitKick(@builtin(global_invocation_id) id:vec3<u32>){
      let i=particleIndex(id);if(i>=config.bodyCount){return;}var p=particlesIn[i];
      if(i>0u){p.velocity+=0.5*config.worldScaleX*(p.acceleration-particlesIn[0].acceleration);}particlesOut[i]=p;
    }

  `;

  const renderWGSL = `
    struct Config {
      bodyCount: u32,
      leafSide: u32,
      leafOffset: u32,
      totalCells: u32,
      depth: u32,
      boundary: u32,
      order: u32,
      interactionBase: u32,
      dt: f32,
      softening2: f32,
      gravitation: f32,
      damping: f32,
      impulseX: f32,
      impulseY: f32,
      impulseForce: f32,
      expansionScale: f32,
      worldScaleX: f32,
      worldScaleY: f32,
      pointScaleX: f32,
      pointScaleY: f32,
      gravityPower: f32,
      viewBrightness: f32,
      pad1: f32,
      pad2: f32,
    }

    struct Particle {
      position: vec2<f32>,
      velocity: vec2<f32>,
      acceleration: vec2<f32>,
      mass: f32,
      color: f32,
    }

    struct VertexOutput {
      @builtin(position) position: vec4<f32>,
      @location(0) uv: vec2<f32>,
      @location(1) color: vec3<f32>,
      @location(2) intensity: f32,
    }

    @group(0) @binding(0) var<uniform> config: Config;
    @group(0) @binding(1) var<storage, read> particles: array<Particle>;

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
      @builtin(instance_index) instanceIndex: u32
    ) -> VertexOutput {
      let corners = array<vec2<f32>, 3>(
        vec2<f32>( 0.0,  2.0),
        vec2<f32>(-1.7320508, -1.0),
        vec2<f32>( 1.7320508, -1.0)
      );
      let particle = particles[instanceIndex];
      let corner = corners[vertexIndex];
      let center = vec2<f32>(
        particle.position.x * config.worldScaleX,
        -particle.position.y * config.worldScaleY
      );
      var output: VertexOutput;
      output.position = vec4<f32>(
        center + corner * vec2<f32>(config.pointScaleX, config.pointScaleY),
        0.0,
        1.0
      );
      output.uv = corner;
      output.color = palette(u32(particle.color + 0.5));
      output.intensity = clamp(
        0.72 + 0.18 * particle.mass * f32(config.bodyCount),
        0.65,
        1.35
      );
      return output;
    }

    @fragment
    fn fragmentMain(input: VertexOutput) -> @location(0) vec4<f32> {
      let distance = length(input.uv);
      if (distance > 1.0) { discard; }
      let halo = pow(max(0.0, 1.0 - distance), 2.2);
      let core = smoothstep(0.32, 0.0, distance);
      let populationScale = clamp(
        pow(8192.0 / f32(config.bodyCount), 0.55),
        0.015,
        1.0
      );
      let alpha = clamp(
        (halo * 0.36 + core * 0.9) *
        input.intensity *
        populationScale *
        config.viewBrightness,
        0.0,
        1.0
      );
      return vec4<f32>(input.color, alpha);
    }
  `;

  const initializeWGSL = `
    const PI: f32 = 3.141592653589793;

    struct InitConfig {
      bodyCount: u32,
      lutSize: u32,
      seed: u32,
      pad0: u32,
      massNormalization: f32,
      gravitation: f32,
      flattening: f32,
      gravityPower: f32,
    }

    struct Particle {
      position: vec2<f32>,
      velocity: vec2<f32>,
      acceleration: vec2<f32>,
      mass: f32,
      color: f32,
    }

    @group(0) @binding(0) var<uniform> config: InitConfig;
    @group(0) @binding(1) var<storage, read> inverseMass: array<f32>;
    @group(0) @binding(2) var<storage, read_write> particles: array<Particle>;

    fn hash(value: u32) -> u32 {
      var result = value;
      result = result ^ (result >> 16u);
      result = result * 0x7feb352du;
      result = result ^ (result >> 15u);
      result = result * 0x846ca68bu;
      result = result ^ (result >> 16u);
      return result;
    }

    fn particleIndex(globalId: vec3<u32>) -> u32 {
      return globalId.x + globalId.y * 16776960u;
    }

    fn random(index: u32, channel: u32) -> f32 {
      let bits = hash(
        index ^ config.seed ^ (channel * 0x9e3779b9u)
      );
      return (f32(bits) + 0.5) / 4294967296.0;
    }

    fn gaussian(index: u32, channel: u32) -> f32 {
      let first = max(0.0000001, random(index, channel));
      let second = random(index, channel + 1u);
      return sqrt(-2.0 * log(first)) * cos(2.0 * PI * second);
    }

    fn sampleMass(index: u32) -> f32 {
      let coordinate =
        random(index, 11u) * f32(config.lutSize - 1u);
      let lower = u32(floor(coordinate));
      let upper = min(lower + 1u, config.lutSize - 1u);
      return mix(
        inverseMass[lower],
        inverseMass[upper],
        fract(coordinate)
      );
    }

    @compute @workgroup_size(256)
    fn initialize(
      @builtin(global_invocation_id) globalId: vec3<u32>
    ) {
      let index = particleIndex(globalId);
      if (index >= config.bodyCount) { return; }

      let isCore = f32(index) < max(8.0, f32(config.bodyCount) * 0.025);
      let radius = select(
        pow(random(index, 1u), 0.58) * 0.83 + 0.035,
        random(index, 1u) * 0.06,
        isCore
      );
      let arm = f32(index % 3u);
      let angle =
        arm * 2.0 * PI / 3.0 +
        radius * 7.8 +
        gaussian(index, 3u) * (0.13 + radius * 0.13);
      let position = vec2<f32>(
        cos(angle) * radius + gaussian(index, 5u) * 0.012,
        sin(angle) * radius * config.flattening +
          gaussian(index, 7u) * 0.012
      );
      let enclosed = 0.16 + min(1.0, radius * 1.35);
      let orbital =
        sqrt(
          abs(config.gravitation) *
          enclosed *
          pow(max(radius, 0.035), 1.0 - config.gravityPower)
        ) * 0.49;
      let velocity = vec2<f32>(
        -sin(angle) * orbital + gaussian(index, 9u) * 0.018,
        cos(angle) * orbital * config.flattening +
          gaussian(index, 13u) * 0.018
      );

      var particle: Particle;
      particle.position = position;
      particle.velocity = velocity;
      particle.acceleration = vec2<f32>(0.0);
      particle.mass = sampleMass(index) * config.massNormalization;
      particle.color = f32(index % 4u);
      particles[index] = particle;
    }
  `;

  async function createOrbitGPUComputeEngine(options) {
    const {
      canvas,
      particles,
      count: requestedCount,
      initialization,
      width,
      height,
      dpr,
      onDeviceLost
    } = options;
    if (!navigator.gpu) {
      throw new Error("WebGPU indisponible dans ce navigateur ou ce contexte");
    }
    const adapter = await navigator.gpu.requestAdapter({
      powerPreference: "high-performance"
    });
    if (!adapter) throw new Error("aucun adaptateur WebGPU disponible");
    const count = requestedCount ?? particles.length;
    const particleBytes = count * PARTICLE_STRIDE;
    const maximumParticleBytes = Math.min(
      Number(adapter.limits.maxStorageBufferBindingSize),
      Number(adapter.limits.maxBufferSize)
    );
    if (particleBytes > maximumParticleBytes) {
      const supportedBodies = Math.floor(
        maximumParticleBytes / PARTICLE_STRIDE
      );
      const limitError = new Error(
        "cette carte limite un buffer à " +
        supportedBodies.toLocaleString("fr-FR") +
        " corps"
      );
      limitError.maximumBodies = supportedBodies;
      throw limitError;
    }
    const device = await adapter.requestDevice({
      requiredLimits: {
        maxStorageBufferBindingSize: particleBytes,
        maxBufferSize: particleBytes
      }
    });
    const context = canvas.getContext("webgpu");
    if (!context) throw new Error("canvas WebGPU indisponible");
    const format = navigator.gpu.getPreferredCanvasFormat();
    const topology = buildTopology(depthForCount(count));

    const computeBindGroupLayout = device.createBindGroupLayout({
      label: "Orbit FMM compute bindings",
      entries: [
        {
          binding: 0,
          visibility: GPUShaderStage.COMPUTE,
          buffer: { type: "uniform" }
        },
        {
          binding: 1,
          visibility: GPUShaderStage.COMPUTE,
          buffer: { type: "read-only-storage" }
        },
        {
          binding: 2,
          visibility: GPUShaderStage.COMPUTE,
          buffer: { type: "storage" }
        },
        {
          binding: 3,
          visibility: GPUShaderStage.COMPUTE,
          buffer: { type: "storage" }
        },
        {
          binding: 4,
          visibility: GPUShaderStage.COMPUTE,
          buffer: { type: "storage" }
        },
        {
          binding: 5,
          visibility: GPUShaderStage.COMPUTE,
          buffer: { type: "storage" }
        },
        {
          binding: 6,
          visibility: GPUShaderStage.COMPUTE,
          buffer: { type: "storage" }
        },
        {
          binding: 7,
          visibility: GPUShaderStage.COMPUTE,
          buffer: { type: "read-only-storage" }
        },
        {
          binding: 8,
          visibility: GPUShaderStage.COMPUTE,
          buffer: { type: "read-only-storage" }
        }
      ]
    });
    const renderBindGroupLayout = device.createBindGroupLayout({
      label: "Orbit FMM render bindings",
      entries: [
        {
          binding: 0,
          visibility: GPUShaderStage.VERTEX | GPUShaderStage.FRAGMENT,
          buffer: { type: "uniform" }
        },
        {
          binding: 1,
          visibility: GPUShaderStage.VERTEX,
          buffer: { type: "read-only-storage" }
        }
      ]
    });
    const computeLayout = device.createPipelineLayout({
      bindGroupLayouts: [computeBindGroupLayout]
    });
    const renderLayout = device.createPipelineLayout({
      bindGroupLayouts: [renderBindGroupLayout]
    });
    device.pushErrorScope("validation");
    const computeModule = device.createShaderModule({
      label: "Orbit FMM compute shaders",
      code: computeWGSL
    });
    const renderModule = device.createShaderModule({
      label: "Orbit FMM render shaders",
      code: renderWGSL
    });
    const initializeModule = device.createShaderModule({
      label: "Orbit particle initialization shader",
      code: initializeWGSL
    });
    const [
      computeCompilation,
      renderCompilation,
      initializeCompilation
    ] = await Promise.all([
      computeModule.getCompilationInfo(),
      renderModule.getCompilationInfo(),
      initializeModule.getCompilationInfo()
    ]);
    const moduleValidationError = await device.popErrorScope();
    const shaderErrors = [
      ...computeCompilation.messages,
      ...renderCompilation.messages,
      ...initializeCompilation.messages
    ].filter((message) => message.type === "error");
    if (moduleValidationError) {
      throw new Error("validation WGSL · " + moduleValidationError.message);
    }
    if (shaderErrors.length) {
      throw new Error(shaderErrors.map((message) =>
        `WGSL ${message.lineNum}:${message.linePos} · ${message.message}`
      ).join(" | "));
    }

    const initializePipeline = initialization
      ? await device.createComputePipelineAsync({
          layout: "auto",
          compute: {
            module: initializeModule,
            entryPoint: "initialize"
          }
        })
      : null;

    const [
      clearPipeline,
      binPipeline,
      p2mPipeline,
      m2mPipeline,
      m2lPipeline,
      l2lPipeline,
      integratePipeline,
      renderPipeline
    ] = await Promise.all([
      ["clear", computeModule],
      ["binParticles", computeModule],
      ["p2m", computeModule],
      ["m2m", computeModule],
      ["m2l", computeModule],
      ["l2l", computeModule],
      ["evaluateAndIntegrate", computeModule]
    ].map(([entryPoint, module]) =>
      device.createComputePipelineAsync({
        layout: computeLayout,
        compute: { module, entryPoint }
      })
    ).concat([
      device.createRenderPipelineAsync({
        layout: renderLayout,
        vertex: {
          module: renderModule,
          entryPoint: "vertexMain"
        },
        fragment: {
          module: renderModule,
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
      })
    ]));

    function createBuffer(label, size, usage, data) {
      const buffer = device.createBuffer({
        label,
        size: Math.max(4, size),
        usage,
        mappedAtCreation: Boolean(data)
      });
      if (data) {
        const destination = new Uint8Array(buffer.getMappedRange());
        const source = data instanceof ArrayBuffer
          ? new Uint8Array(data)
          : new Uint8Array(data.buffer, data.byteOffset, data.byteLength);
        destination.set(source);
        buffer.unmap();
      }
      return buffer;
    }

    const initialParticles = initialization
      ? null
      : packParticles(particles);
    const particleUsage =
      GPUBufferUsage.STORAGE |
      GPUBufferUsage.COPY_DST |
      GPUBufferUsage.COPY_SRC;
    const particleBuffers = [
      createBuffer(
        "Orbit particles A",
        count * PARTICLE_STRIDE,
        particleUsage,
        initialParticles
      ),
      createBuffer(
        "Orbit particles B",
        count * PARTICLE_STRIDE,
        particleUsage,
        initialParticles
      )
    ];
    if (initialization) {
      const massLut = initialization.massLut;
      const massLutBuffer = createBuffer(
        "Orbit inverse mass CDF",
        massLut.byteLength,
        GPUBufferUsage.STORAGE,
        massLut
      );
      const initializationData = new ArrayBuffer(32);
      const initializationInts = new Uint32Array(initializationData);
      const initializationFloats = new Float32Array(initializationData);
      initializationInts[0] = count;
      initializationInts[1] = massLut.length;
      initializationInts[2] = initialization.seed >>> 0;
      initializationFloats[4] =
        1 / (count * initialization.inverseMean);
      initializationFloats[5] = initialization.gravitation ?? 0.48;
      initializationFloats[6] = 0.73;
      initializationFloats[7] = initialization.gravityPower ?? 2;
      const initializationBuffer = createBuffer(
        "Orbit particle initialization parameters",
        initializationData.byteLength,
        GPUBufferUsage.UNIFORM,
        initializationData
      );
      const initializationBindGroup = device.createBindGroup({
        layout: initializePipeline.getBindGroupLayout(0),
        entries: [
          { binding: 0, resource: { buffer: initializationBuffer } },
          { binding: 1, resource: { buffer: massLutBuffer } },
          { binding: 2, resource: { buffer: particleBuffers[0] } }
        ]
      });
      const initializationEncoder = device.createCommandEncoder({
        label: "Orbit GPU particle initialization"
      });
      const initializationPass = initializationEncoder.beginComputePass();
      initializationPass.setPipeline(initializePipeline);
      initializationPass.setBindGroup(0, initializationBindGroup);
      dispatchParticleWorkgroups(initializationPass, count);
      initializationPass.end();
      initializationEncoder.copyBufferToBuffer(
        particleBuffers[0],
        0,
        particleBuffers[1],
        0,
        particleBytes
      );
      device.queue.submit([initializationEncoder.finish()]);
      await device.queue.onSubmittedWorkDone();
      initializationBuffer.destroy();
      massLutBuffer.destroy();
    }
    const configBuffer = createBuffer(
      "Orbit FMM configuration",
      96,
      GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST
    );
    const leafAccumBuffer = createBuffer(
      "Orbit leaf atomic accumulators",
      topology.leafCount * 16,
      GPUBufferUsage.STORAGE
    );
    const leafSlotsBuffer = createBuffer(
      "Orbit bounded direct lists",
      topology.leafCount * DIRECT_LIMIT * 4,
      GPUBufferUsage.STORAGE
    );
    const particleCellBuffer = createBuffer(
      "Orbit particle cell ids",
      count * 4,
      GPUBufferUsage.STORAGE
    );
    const cellBuffer = createBuffer(
      "Orbit FMM cells",
      topology.totalCells * CELL_STRIDE,
      GPUBufferUsage.STORAGE
    );
    const topologyBuffer = createBuffer(
      "Orbit FMM topology",
      topology.topologyData.byteLength,
      GPUBufferUsage.STORAGE,
      new Uint8Array(topology.topologyData)
    );
    const interactionBuffer = createBuffer(
      "Orbit FMM interaction lists",
      topology.interactionData.byteLength,
      GPUBufferUsage.STORAGE,
      topology.interactionData
    );

    const computeBindGroups = particleBuffers.map((input, index) =>
      device.createBindGroup({
        layout: computeBindGroupLayout,
        entries: [
          { binding: 0, resource: { buffer: configBuffer } },
          { binding: 1, resource: { buffer: input } },
          { binding: 2, resource: { buffer: particleBuffers[1 - index] } },
          { binding: 3, resource: { buffer: leafAccumBuffer } },
          { binding: 4, resource: { buffer: leafSlotsBuffer } },
          { binding: 5, resource: { buffer: particleCellBuffer } },
          { binding: 6, resource: { buffer: cellBuffer } },
          { binding: 7, resource: { buffer: topologyBuffer } },
          { binding: 8, resource: { buffer: interactionBuffer } }
        ]
      })
    );
    const renderBindGroups = particleBuffers.map((buffer) =>
      device.createBindGroup({
        layout: renderBindGroupLayout,
        entries: [
          { binding: 0, resource: { buffer: configBuffer } },
          { binding: 1, resource: { buffer } }
        ]
      })
    );
    const configData = new ArrayBuffer(96);
    const configInts = new Uint32Array(configData);
    const configFloats = new Float32Array(configData);
    let current = 0;
    let frameIndex = 0;
    let lastGpuMs = 0;
    let disposed = false;
    let cssWidth = width;
    let cssHeight = height;
    let pixelRatio = dpr;

    function configure() {
      context.configure({
        device,
        format,
        alphaMode: "opaque"
      });
    }
    configure();

    function updateConfig(frame) {
      const boundary =
        frame.boundary === "torus"
          ? 2
          : frame.boundary === "reflect"
            ? 1
            : 0;
      const minDimension = Math.min(cssWidth, cssHeight);
      const zoom = frame.zoom ?? 1;
      const viewBrightness = Math.max(
        0.55,
        Math.min(8, zoom ** 0.75)
      );
      const pointZoom = Math.max(
        0.8,
        Math.min(2, zoom ** 0.2)
      );
      const radius =
        count > 4000000
          ? 0.85
          : count > 1000000
            ? 1.15
            : count > 65536
              ? 1.65
              : 2.4;
      configInts[0] = count;
      configInts[1] = topology.leafSide;
      configInts[2] = topology.leafOffset;
      configInts[3] = topology.totalCells;
      configInts[4] = topology.depth;
      configInts[5] = boundary;
      configInts[6] = frame.order;
      configInts[7] =
        boundary === 2
          ? topology.totalCells * INTERACTIONS_PER_CELL
          : 0;
      configFloats[8] =
        Math.min(0.012, frame.dt * 0.00038) * frame.speed;
      configFloats[9] = frame.softening2 ?? 0.0018;
      configFloats[10] = frame.gravitation ?? 0.48;
      configFloats[11] = 0.99996;
      configFloats[12] = frame.impulse ? frame.impulse.x : 0;
      configFloats[13] = frame.impulse ? frame.impulse.y : 0;
      configFloats[14] = frame.impulse ? frame.impulse.force : 0;
      configFloats[15] = frame.expansionScale ?? 1;
      configFloats[16] =
        0.96 * minDimension / cssWidth * (frame.zoom ?? 1);
      configFloats[17] =
        0.96 * minDimension / cssHeight * (frame.zoom ?? 1);
      configFloats[18] =
        2 * radius * pointZoom * pixelRatio / canvas.width;
      configFloats[19] =
        2 * radius * pointZoom * pixelRatio / canvas.height;
      configFloats[20] = frame.gravityPower ?? 2;
      configFloats[21] = viewBrightness;
      configFloats[22] = 0;
      configFloats[23] = 0;
      device.queue.writeBuffer(configBuffer, 0, configData);
    }

    const engine = {
      count,
      depth: topology.depth,
      totalCells: topology.totalCells,
      leafCount: topology.leafCount,
      lastGpuMs: 0,
      configure(nextWidth, nextHeight, nextDpr) {
        cssWidth = nextWidth;
        cssHeight = nextHeight;
        pixelRatio = nextDpr;
        configure();
      },
      reset(nextParticles) {
        if (nextParticles.length !== count) return false;
        const packed = packParticles(nextParticles);
        device.queue.writeBuffer(particleBuffers[0], 0, packed);
        device.queue.writeBuffer(particleBuffers[1], 0, packed);
        current = 0;
        return true;
      },
      stepAndRender(frame) {
        const submitStarted = performance.now();
        updateConfig(frame);
        const encoder = device.createCommandEncoder({
          label: "Orbit GPU FMM frame"
        });
        const computePass = encoder.beginComputePass({
          label: "Orbit GPU FMM passes"
        });
        computePass.setBindGroup(0, computeBindGroups[current]);
        computePass.setPipeline(clearPipeline);
        computePass.dispatchWorkgroups(
          Math.ceil(Math.max(topology.totalCells, topology.leafCount) / WORKGROUP_SIZE)
        );
        computePass.setPipeline(binPipeline);
        dispatchParticleWorkgroups(computePass, count);
        computePass.setPipeline(p2mPipeline);
        computePass.dispatchWorkgroups(Math.ceil(topology.leafCount / WORKGROUP_SIZE));
        computePass.setPipeline(m2mPipeline);
        computePass.dispatchWorkgroups(1);
        computePass.setPipeline(m2lPipeline);
        computePass.dispatchWorkgroups(Math.ceil(topology.totalCells / WORKGROUP_SIZE));
        computePass.setPipeline(l2lPipeline);
        computePass.dispatchWorkgroups(1);
        computePass.setPipeline(integratePipeline);
        dispatchParticleWorkgroups(computePass, count);
        computePass.end();

        const output = 1 - current;
        const renderPass = encoder.beginRenderPass({
          colorAttachments: [{
            view: context.getCurrentTexture().createView(),
            clearValue: { r: 0.019, g: 0.027, b: 0.043, a: 1 },
            loadOp: "clear",
            storeOp: "store"
          }]
        });
        renderPass.setPipeline(renderPipeline);
        renderPass.setBindGroup(0, renderBindGroups[output]);
        renderPass.draw(3, count);
        renderPass.end();
        device.queue.submit([encoder.finish()]);
        current = output;
        frameIndex++;

        if (frameIndex % 45 === 0) {
          const gpuStarted = performance.now();
          device.queue.onSubmittedWorkDone().then(() => {
            lastGpuMs = performance.now() - gpuStarted;
            engine.lastGpuMs = lastGpuMs;
          });
        }
        return {
          submitMs: performance.now() - submitStarted,
          gpuMs: lastGpuMs
        };
      },
      destroy() {
        if (disposed) return;
        disposed = true;
        for (const buffer of particleBuffers) buffer.destroy();
        configBuffer.destroy();
        leafAccumBuffer.destroy();
        leafSlotsBuffer.destroy();
        particleCellBuffer.destroy();
        cellBuffer.destroy();
        topologyBuffer.destroy();
        interactionBuffer.destroy();
        device.destroy();
      }
    };

    device.lost.then((info) => {
      if (!disposed && onDeviceLost) onDeviceLost(info);
    });
    return engine;
  }

  // Headless force adapter for local stellar dynamics. It deliberately shares
  // computeWGSL, topology construction and every FMM pass with the original app.
  // No renderer, damping or GPU device ownership is added by this adapter.
  async function createOrbitGPUForceField({device,capacity=8192}) {
    const U=GPUBufferUsage,buffers=[],topologies=new Map();let disposed=false,lastMs=0;
    function buffer(label,size,usage,data){const b=device.createBuffer({label,size:Math.max(16,size),usage});buffers.push(b);if(data)device.queue.writeBuffer(b,0,data);return b;}
    device.pushErrorScope('validation');
    try {
      const module=device.createShaderModule({label:'Orbit shared FMM · local forces',code:computeWGSL}),info=await module.getCompilationInfo(),errors=info.messages.filter(m=>m.type==='error');if(errors.length)throw Error(errors.map(m=>m.message).join('\n'));
      const layout=device.createBindGroupLayout({entries:Array.from({length:9},(_,binding)=>({binding,visibility:GPUShaderStage.COMPUTE,buffer:{type:binding===0?'uniform':[1,7,8].includes(binding)?'read-only-storage':'storage'}}))}),pipelineLayout=device.createPipelineLayout({bindGroupLayouts:[layout]}),pipelines={};
      for(const entryPoint of ['clear','binParticles','p2m','m2m','m2l','l2l','evaluateAndIntegrate'])pipelines[entryPoint]=await device.createComputePipelineAsync({layout:pipelineLayout,compute:{module,entryPoint}});
      const maxDepth=Math.min(6,Math.ceil(Math.log(Math.max(1,capacity/4))/Math.log(4))),maxCells=levelOffset(maxDepth+1),maxLeaves=4**maxDepth;
      const configs=[0,1].map(i=>buffer('Orbit local force config '+i,96,U.UNIFORM|U.COPY_DST)),input=buffer('Orbit local force particles',capacity*PARTICLE_STRIDE,U.STORAGE|U.COPY_DST),outputs=[0,1].map(i=>buffer('Orbit local force output '+i,capacity*PARTICLE_STRIDE,U.STORAGE|U.COPY_SRC)),accum=buffer('Orbit local force accumulators',maxLeaves*16,U.STORAGE),slots=buffer('Orbit local force near lists',maxLeaves*DIRECT_LIMIT*4,U.STORAGE),particleCells=buffer('Orbit local force cell IDs',capacity*4,U.STORAGE),cells=buffer('Orbit local force multipoles',maxCells*CELL_STRIDE,U.STORAGE),readback=buffer('Orbit local force readback',capacity*PARTICLE_STRIDE*2,U.MAP_READ|U.COPY_DST);
      function topologyFor(depth){if(topologies.has(depth))return topologies.get(depth);const topology=buildTopology(depth),tb=buffer('Orbit local topology '+depth,topology.topologyData.byteLength,U.STORAGE|U.COPY_DST,new Uint8Array(topology.topologyData)),ib=buffer('Orbit local interactions '+depth,topology.interactionData.byteLength,U.STORAGE|U.COPY_DST,topology.interactionData);const groups=outputs.map((output,k)=>device.createBindGroup({layout,entries:[configs[k],input,output,accum,slots,particleCells,cells,tb,ib].map((b,binding)=>({binding,resource:{buffer:b}}))}));const result={...topology,groups};topologies.set(depth,result);return result;}
      const validation=await device.popErrorScope();if(validation)throw Error(validation.message);
      return {
        capacity,get lastMs(){return lastMs;},backend:'Orbit GPU FMM',
        async compute(particles,{gravity,softening,coarseSoftening}){
          if(disposed)throw Error('FMM GPU locale détruite.');const n=particles.length;if(n>capacity)throw Error('Budget FMM GPU locale dépassé.');const passes=coarseSoftening===undefined?1:2;const sum=particles.reduce((s,p)=>s+p.mass,0);if(!n||!(sum>0))return new Float32Array(n*2);const started=performance.now(),depth=Math.max(0,Math.min(6,Math.ceil(Math.log(Math.max(1,n/4))/Math.log(4)))),topology=topologyFor(depth),packed=packParticles(particles.map(p=>({x:p.x/1.08,y:p.y/1.08,vx:0,vy:0,mass:p.mass/sum,color:0})));device.queue.writeBuffer(input,0,packed);
          for(let k=0;k<passes;k++){const config=new ArrayBuffer(96),ints=new Uint32Array(config),floats=new Float32Array(config);ints.set([n,topology.leafSide,topology.leafOffset,topology.totalCells,depth,0,1,0]);floats[9]=((k?coarseSoftening:softening)/1.08)**2;floats[10]=gravity*sum/1.08**3;floats[11]=1;floats[15]=1;floats[20]=2;floats[22]=1;device.queue.writeBuffer(configs[k],0,config);}
          const e=device.createCommandEncoder({label:'Orbit local fine/coarse FMM'});
          for(let k=0;k<passes;k++){const pass=e.beginComputePass();pass.setBindGroup(0,topology.groups[k]);if(k===0){pass.setPipeline(pipelines.clear);pass.dispatchWorkgroups(Math.ceil(Math.max(topology.totalCells,topology.leafCount)/WORKGROUP_SIZE));pass.setPipeline(pipelines.binParticles);dispatchParticleWorkgroups(pass,n);}pass.setPipeline(pipelines.p2m);pass.dispatchWorkgroups(Math.ceil(topology.leafCount/WORKGROUP_SIZE));pass.setPipeline(pipelines.m2m);pass.dispatchWorkgroups(1);pass.setPipeline(pipelines.m2l);pass.dispatchWorkgroups(Math.ceil(topology.totalCells/WORKGROUP_SIZE));pass.setPipeline(pipelines.l2l);pass.dispatchWorkgroups(1);pass.setPipeline(pipelines.evaluateAndIntegrate);dispatchParticleWorkgroups(pass,n);pass.end();e.copyBufferToBuffer(outputs[k],0,readback,k*n*PARTICLE_STRIDE,n*PARTICLE_STRIDE);}
          device.queue.submit([e.finish()]);await readback.mapAsync(GPUMapMode.READ,0,n*PARTICLE_STRIDE*passes);const values=new Float32Array(readback.getMappedRange(0,n*PARTICLE_STRIDE*passes)),result=new Float32Array(n*2);for(let i=0;i<n;i++){result[i*2]=(values[i*8+4]-(passes===2?values[n*8+i*8+4]:0))*1.08;result[i*2+1]=(values[i*8+5]-(passes===2?values[n*8+i*8+5]:0))*1.08;}readback.unmap();lastMs=performance.now()-started;return result;
        },
        destroy(){if(disposed)return;disposed=true;for(const b of buffers)b.destroy();}
      };
    } catch(error) {await device.popErrorScope().catch(()=>{});for(const b of buffers)b.destroy();throw error;}
  }

  async function createOrbitGPUOrbitSystem({device,bodies,lengthScale}){
    const U=GPUBufferUsage,buffers=[],n=bodies.length,depth=Math.max(0,Math.min(6,Math.ceil(Math.log(Math.max(1,n/4))/Math.log(4)))),topology=buildTopology(depth),scale=lengthScale||Math.max(1,...bodies.map(p=>Math.max(Math.abs(p.x),Math.abs(p.y))))*1.3;
    if(n===1)return{backend:'Orbit GPU FMM',scale,count:1,lastMs:0,async advance(dt,steps=1){return{bodies:bodies.map(p=>({...p})),advanced:dt*steps,steps};},destroy(){}};if(n<2||n>32768)throw Error('Le système local doit contenir entre 1 et 32768 corps.');
    const mass=bodies.slice(1).reduce((s,p)=>s+p.mass,0);if(!(mass>0&&bodies[0].mass>0))throw Error('Masses orbitales invalides.');
    const buffer=(label,size,usage,data)=>{const b=device.createBuffer({label,size:Math.max(16,size),usage});buffers.push(b);if(data)device.queue.writeBuffer(b,0,data);return b;};
    device.pushErrorScope('validation');let scoped=true;
    try{
      const module=device.createShaderModule({label:'Orbit shared GPU FMM · planetary KDK',code:computeWGSL}),info=await module.getCompilationInfo();if(info.messages.some(m=>m.type==='error'))throw Error(info.messages.filter(m=>m.type==='error').map(m=>m.message).join('\n'));
      const layout=device.createBindGroupLayout({entries:Array.from({length:9},(_,binding)=>({binding,visibility:GPUShaderStage.COMPUTE,buffer:{type:binding===0?'uniform':[1,7,8].includes(binding)?'read-only-storage':'storage'}}))}),pipelineLayout=device.createPipelineLayout({bindGroupLayouts:[layout]}),pipelines={};
      for(const entryPoint of ['clear','binParticles','p2m','m2m','m2l','l2l','evaluateAndIntegrate','orbitDrift','orbitKick'])pipelines[entryPoint]=await device.createComputePipelineAsync({layout:pipelineLayout,compute:{module,entryPoint}});
      const config=buffer('Planetary force settings',96,U.UNIFORM|U.COPY_DST),particles=[0,1].map(i=>buffer('Planetary GPU states '+i,n*PARTICLE_STRIDE,U.STORAGE|U.COPY_SRC|U.COPY_DST)),accum=buffer('Planetary FMM bins',topology.leafCount*16,U.STORAGE),slots=buffer('Planetary near lists',topology.leafCount*(DIRECT_LIMIT+1)*4,U.STORAGE),particleCells=buffer('Planetary cells',n*8,U.STORAGE),cells=buffer('Planetary multipoles',topology.totalCells*CELL_STRIDE,U.STORAGE),tb=buffer('Planetary shared topology',topology.topologyData.byteLength,U.STORAGE|U.COPY_DST,new Uint8Array(topology.topologyData)),ib=buffer('Planetary interactions',topology.interactionData.byteLength,U.STORAGE|U.COPY_DST,topology.interactionData),readback=buffer('Planetary snapshot',n*PARTICLE_STRIDE,U.MAP_READ|U.COPY_DST),groups=particles.map((p,i)=>device.createBindGroup({layout,entries:[config,p,particles[1-i],accum,slots,particleCells,cells,tb,ib].map((b,binding)=>({binding,resource:{buffer:b}}))}));
      device.queue.writeBuffer(particles[0],0,packParticles(bodies.map((p,i)=>({...p,x:p.x/scale,y:p.y/scale,vx:p.vx/scale,vy:p.vy/scale,mass:i?p.mass/mass:0,color:0}))));
      const error=await device.popErrorScope();scoped=false;if(error)throw Error(error.message);let current=0,initialized=false,disposed=false;
      const writeConfig=dt=>{const b=new ArrayBuffer(96),u=new Uint32Array(b),f=new Float32Array(b);u.set([n,topology.leafSide,topology.leafOffset,topology.totalCells,depth,0,1,0]);f[9]=(1e-7/scale)**2;f[10]=39.47841760435743*mass/scale**3;f[11]=f[15]=1;f[16]=dt;f[20]=2;f[22]=1;f[23]=bodies[0].mass/mass;device.queue.writeBuffer(config,0,b);};
      function force(e){const pass=e.beginComputePass();pass.setBindGroup(0,groups[current]);for(const [name,count]of [['clear',Math.ceil(Math.max(topology.totalCells,topology.leafCount)/256)],['binParticles',Math.ceil(n/256)],['p2m',Math.ceil(topology.leafCount/256)],['m2m',1],['m2l',Math.ceil(topology.totalCells/256)],['l2l',1],['evaluateAndIntegrate',Math.ceil(n/256)]]){pass.setPipeline(pipelines[name]);pass.dispatchWorkgroups(count);}pass.end();current=1-current;}
      function kick(e,name){const pass=e.beginComputePass();pass.setBindGroup(0,groups[current]);pass.setPipeline(pipelines[name]);pass.dispatchWorkgroups(Math.ceil(n/256));pass.end();current=1-current;}
      return{backend:'Orbit GPU FMM',scale,count:n,lastMs:0,async advance(dt,steps=1){if(disposed)throw Error('Système GPU détruit.');if(!(dt>0)||!Number.isInteger(steps)||steps<1||steps>128)throw Error('Pas orbital invalide.');const started=performance.now();writeConfig(dt);const e=device.createCommandEncoder({label:'Planetary shared FMM KDK batch'});if(!initialized){force(e);initialized=true;}for(let i=0;i<steps;i++){kick(e,'orbitDrift');force(e);kick(e,'orbitKick');}e.copyBufferToBuffer(particles[current],0,readback,0,n*PARTICLE_STRIDE);device.queue.submit([e.finish()]);await readback.mapAsync(GPUMapMode.READ);const data=new Float32Array(readback.getMappedRange()),result=bodies.map((p,i)=>({...p,x:data[i*8]*scale,y:data[i*8+1]*scale,vx:data[i*8+2]*scale,vy:data[i*8+3]*scale,ax:data[i*8+4]*scale,ay:data[i*8+5]*scale}));readback.unmap();this.lastMs=performance.now()-started;return{bodies:result,advanced:dt*steps,steps};},destroy(){if(disposed)return;disposed=true;for(const b of buffers)b.destroy();}};
    }catch(error){if(scoped)await device.popErrorScope().catch(()=>{});for(const b of buffers)b.destroy();throw error;}
  }
  window.createOrbitGPUOrbitSystem=createOrbitGPUOrbitSystem;
  window.createOrbitGPUComputeEngine = createOrbitGPUComputeEngine;
  window.createOrbitGPUForceField = createOrbitGPUForceField;
})();
