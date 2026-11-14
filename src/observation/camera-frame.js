/* Coordinate and publication boundary shared by CPU and GPU presentation. */
(function (root) {
  "use strict";
  // Each image owns its coordinate frame. Origins use model coordinates while
  // centres/anchors use the frame's unit. Rebase near the host for AU precision.
  // The physical field spans the longest viewport axis; pixels stay isotropic.
  function viewportGeometry(
    width = 600,
    height = 600,
    { pixelRatio = 1, size = null, maxDimension = 8192 } = {},
  ) {
    width = Math.max(1, Number(width) || 600);
    height = Math.max(1, Number(height) || 600);
    const span = Math.max(width, height),
      scale =
        size === null
          ? Math.min(Math.max(0.1, pixelRatio), maxDimension / span)
          : size / span,
      pixelWidth = Math.max(1, Math.round(width * scale)),
      pixelHeight = Math.max(1, Math.round(height * scale));
    return {
      size: Math.max(pixelWidth, pixelHeight),
      pixelSpan: Math.max(pixelWidth, pixelHeight),
      cssWidth: width,
      cssHeight: height,
      cssShort: Math.min(width, height),
      aspectX: width / span,
      aspectY: height / span,
      pixelWidth,
      pixelHeight,
      referenceWidth: (600 * width) / span,
      referenceHeight: (600 * height) / span,
      cssSpan: span,
      pixelRatio: scale,
    };
  }
  function renderContext(
    viewport,
    {
      worldShortSpan = 1,
      worldSpan = null,
      fovDegrees = 50,
      sigmaCss = 1,
      photometricReferenceCss = 600,
      units = "model",
      lengthPc = 1,
      observationErrorCss = 1,
    } = {},
  ) {
    const cssWidth = viewport.cssSpan * viewport.aspectX,
      cssHeight = viewport.cssSpan * viewport.aspectY,
      cssShort = Math.min(cssWidth, cssHeight),
      pixelSpan = viewport.size;
    worldSpan = worldSpan ?? (worldShortSpan * viewport.cssSpan) / cssShort;
    const sensorPixels = (pixelSpan * worldShortSpan) / worldSpan,
      pixelScale = viewport.pixelRatio;
    return Object.freeze({
      ...viewport,
      cssWidth,
      cssHeight,
      cssShort,
      pixelSpan,
      worldShortSpan,
      worldSpan,
      worldWidth: worldSpan * viewport.aspectX,
      worldHeight: worldSpan * viewport.aspectY,
      worldPerCssPixel: worldSpan / viewport.cssSpan,
      worldPerPixel: worldSpan / pixelSpan,
      sensorPixels,
      focalPixels: sensorPixels / (2 * Math.tan((fovDegrees * Math.PI) / 360)),
      sigmaCss,
      sigmaPixels: sigmaCss * pixelScale,
      exposurePixels: photometricReferenceCss * pixelScale,
      diffuseScale:
        ((photometricReferenceCss * pixelScale) / sensorPixels) ** 2,
      observationTolerance:
        (observationErrorCss * worldSpan) / viewport.cssSpan,
      units,
      lengthPc,
      fovDegrees,
    });
  }
  // Coverage is checked in destination pixels. This never enlarges the rendered
  // field, crops the viewport or requests another observation/render.
  function coveredCameraView(image, requested, imageContext, displayContext) {
    const covers = (view) => {
      const t = cameraPresentation(image, view, displayContext.pixelSpan),
        scale = (t.scale * displayContext.pixelSpan) / imageContext.pixelSpan,
        w = imageContext.pixelWidth * scale,
        h = imageContext.pixelHeight * scale,
        cx = displayContext.pixelWidth / 2 + t.x,
        cy = displayContext.pixelHeight / 2 + t.y,
        e = 1e-5;
      return (
        cx - w / 2 <= e &&
        cy - h / 2 <= e &&
        cx + w / 2 >= displayContext.pixelWidth - e &&
        cy + h / 2 >= displayContext.pixelHeight - e
      );
    };
    if (covers(requested)) return { view: requested, pending: false };
    return { view: image, pending: true };
  }

  function viewportProjection(projection, aspectX = 1, aspectY = 1) {
    const n = projection.resolution,
      loX = ((1 - aspectX) * n) / 2,
      hiX = ((1 + aspectX) * n) / 2,
      loY = ((1 - aspectY) * n) / 2,
      hiY = ((1 + aspectY) * n) / 2;
    let mass = 0;
    for (let y = Math.floor(loY); y < Math.ceil(hiY); y++)
      for (let x = Math.floor(loX); x < Math.ceil(hiX); x++)
        mass +=
          projection.image[y * n + x] *
          Math.max(0, Math.min(x + 1, hiX) - Math.max(x, loX)) *
          Math.max(0, Math.min(y + 1, hiY) - Math.max(y, loY));
    return { ...projection, aspectX, aspectY, viewMass: mass };
  }
  function cameraPresentation(image, view, size) {
    const delta = (a, b) => ({
      x:
        (a.originX - b.originX) * a.modelLengthPc +
        a.x * a.unitPc -
        b.x * b.unitPc,
      y:
        (a.originY - b.originY) * a.modelLengthPc +
        a.y * a.unitPc -
        b.y * b.unitPc,
    });
    const shift = delta(image, view),
      width = view.width * view.unitPc;
    if (image.follow && view.follow && image.follow.id === view.follow.id) {
      const movement = delta(
        { ...view, ...view.follow },
        { ...image, ...image.follow },
      );
      shift.x += movement.x;
      shift.y += movement.y;
    }
    return {
      scale: (image.width * image.unitPc) / width,
      x: (shift.x / width) * size,
      y: (-shift.y / width) * size,
    };
  }
  // Reproject a committed scene using navigation only. Physics may mutate live
  // state while a GPU calculation is pending; it must never move an older image.
  function committedCameraView(image, view, pan = { x: 0, y: 0 }) {
    if (image.space !== view.space || image.referenceId !== view.referenceId)
      return view;
    const result = { ...view, originX: image.originX, originY: image.originY };
    if (view.space !== "system") {
      result.x +=
        ((view.originX - image.originX) * view.modelLengthPc) / view.unitPc;
      result.y +=
        ((view.originY - image.originY) * view.modelLengthPc) / view.unitPc;
    }
    if (image.follow && view.follow?.id === image.follow.id) {
      result.follow = image.follow;
      result.x = image.follow.x + pan.x;
      result.y = image.follow.y + pan.y;
    } else result.follow = null;
    return result;
  }

  /**
   * Owns the last complete image description. Large buffers are borrowed, not
   * copied: producers must replace them rather than mutate a published image.
   * Live camera/physics objects are never retained as mutable frame metadata.
   */
  class FramePublication {
    constructor() {
      this.current = null;
    }
    clear() {
      this.current = null;
    }
    publish({ camera, context, points, projection, solarScene, revision }) {
      const frameCamera = Object.freeze({
        ...camera,
        follow: camera.follow ? Object.freeze({ ...camera.follow }) : null,
      });
      this.current = Object.freeze({
        camera: frameCamera,
        context,
        points,
        projection,
        solarScene,
        revision,
      });
      return this.current;
    }
  }
  const api = {
    viewportGeometry,
    renderContext,
    coveredCameraView,
    viewportProjection,
    cameraPresentation,
    committedCameraView,
    FramePublication,
  };
  Object.assign(root, api);
  if (typeof module !== "undefined") module.exports = api;
})(typeof window !== "undefined" ? window : globalThis);
