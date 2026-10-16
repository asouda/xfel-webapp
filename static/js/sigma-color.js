/* Colors for the single interactive isosurface.
   Signed: positive sigma is blue and negative sigma is red, darker farther from zero.
   Temperature: cold blue at the low end of the range, hot red at the high end. */
(function (root, factory) {
  const api = factory();
  if (typeof module === "object" && module.exports) module.exports = api;
  root.XFELSigma = api;
})(typeof globalThis !== "undefined" ? globalThis : this, function () {
  "use strict";

  const COLD_LOW = "#5aa7ff";
  const COLD_HIGH = "#031a4a";
  const WARM_LOW = "#ff8d78";
  const WARM_HIGH = "#6e0812";
  const ZERO = "#b7bec8";
  const TEMPERATURE = [
    [0, "#0a1f6e"],
    [0.18, "#1a62d0"],
    [0.36, "#1aa8c9"],
    [0.52, "#e6c84a"],
    [0.7, "#ef7a12"],
    [0.86, "#e23b1a"],
    [1, "#8d1010"],
  ];

  function clamp01(value) {
    return Math.min(1, Math.max(0, value));
  }

  function hexToRgb(hex) {
    const n = parseInt(hex.slice(1), 16);
    return [(n >> 16) & 255, (n >> 8) & 255, n & 255];
  }

  function rgbToHex(r, g, b) {
    const byte = (value) => Math.round(value).toString(16).padStart(2, "0");
    return `#${byte(r)}${byte(g)}${byte(b)}`;
  }

  function lerpColor(start, end, t) {
    const u = clamp01(t);
    const a = hexToRgb(start);
    const b = hexToRgb(end);
    return rgbToHex(
      a[0] + (b[0] - a[0]) * u,
      a[1] + (b[1] - a[1]) * u,
      a[2] + (b[2] - a[2]) * u,
    );
  }

  function lerpStops(stops, t) {
    const u = clamp01(t);
    let index = 0;
    while (index < stops.length - 2 && u > stops[index + 1][0]) index += 1;
    const [t0, c0] = stops[index];
    const [t1, c1] = stops[index + 1];
    const span = t1 - t0 || 1;
    return lerpColor(c0, c1, (u - t0) / span);
  }

  function signedColor(sigma, minimum, maximum) {
    if (sigma > 0) {
      const start = Math.max(minimum, 0);
      const end = Math.max(maximum, start);
      const t = end === start ? 1 : (sigma - start) / (end - start);
      return lerpColor(COLD_LOW, COLD_HIGH, t);
    }
    if (sigma < 0) {
      const start = Math.min(maximum, 0);
      const end = Math.min(minimum, start);
      const t = start === end ? 1 : (start - sigma) / (start - end);
      return lerpColor(WARM_LOW, WARM_HIGH, t);
    }
    return ZERO;
  }

  function temperatureColor(sigma, minimum, maximum) {
    const span = maximum - minimum;
    const t = span === 0 ? 0.5 : (sigma - minimum) / span;
    return lerpStops(TEMPERATURE, t);
  }

  function colorForSigma(sigma, minimum, maximum, scale) {
    if (scale === "temperature") return temperatureColor(sigma, minimum, maximum);
    return signedColor(sigma, minimum, maximum);
  }

  function formatSigma(value) {
    if (!Number.isFinite(value)) return "";
    return String(Math.round(value * 10000) / 10000);
  }

  return { colorForSigma, formatSigma };
});
