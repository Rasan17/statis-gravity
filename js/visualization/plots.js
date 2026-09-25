/**
 * Statis-Gravity - Statistical Plots Suite
 * Renders Box-and-Whisker, Distribution Histograms, Scatter + Regression CI, and ROC Curves.
 */

import { ChartEngine } from './chart-engine.js';
import { Distributions } from '../stats/distributions.js';
import { Teaching } from '../stats/teaching.js';

export const Plots = {
  /**
   * Renders Box & Whisker plot for one or more groups
   * @param {ChartEngine} engine
   * @param {Array<{name: string, stats: object, color?: string}>} groupStats
   * @param {string} title
   */
  renderBoxPlot(engine, groupStats, title = 'Distribution & Box-and-Whisker Plot') {
    engine.lastRenderFn = () => this.renderBoxPlot(engine, groupStats, title);
    engine.clear();
    const b = engine.getPlotBounds ? engine.getPlotBounds() : engine.getBounds();
    const ctx = engine.ctx;
    const pal = engine.palette;

    let groups = [];
    if (Array.isArray(groupStats)) {
      groups = groupStats.map(g => g.stats ? g : { name: g.name || 'Sample Data', stats: g, color: g.color });
    } else if (groupStats && typeof groupStats === 'object') {
      groups = [{ name: groupStats.name || 'Sample Data', stats: groupStats, color: pal.primary }];
    }

    if (groups.length === 0) return;

    // Find global min and max across all groups
    let globalMin = Infinity;
    let globalMax = -Infinity;
    for (const g of groups) {
      if (!g.stats || g.stats.n === 0) continue;
      if (g.stats.min < globalMin) globalMin = g.stats.min;
      if (g.stats.max > globalMax) globalMax = g.stats.max;
    }

    if (!isFinite(globalMin) || !isFinite(globalMax)) return;

    // Margin for Y axis
    const span = globalMax - globalMin || 1;
    const yMin = globalMin - 0.1 * span;
    const yMax = globalMax + 0.1 * span;
    const yRange = yMax - yMin;

    const yToPixel = (val) => b.y + b.height - ((val - yMin) / yRange) * b.height;

    // Build Y ticks
    const yTicks = [];
    const stepCount = 5;
    for (let i = 0; i <= stepCount; i++) {
      const v = yMin + (i / stepCount) * yRange;
      yTicks.push({ norm: (v - yMin) / yRange, label: v.toFixed(1) });
    }

    // Build X ticks
    const k = groups.length;
    const xTicks = groups.map((g, idx) => ({
      norm: (idx + 0.5) / k,
      label: g.name
    }));

    engine.drawAxes({ yTicks, xTicks, title, yLabel: 'Observed Value' });

    // Draw each box
    const slotWidth = b.width / k;
    const boxWidth = Math.min(55, slotWidth * 0.45);

    const colors = [pal.primary, pal.secondary, pal.accent, '#f59e0b'];

    groups.forEach((g, idx) => {
      const s = g.stats;
      if (!s || s.n === 0) return;

      const groupColor = g.color || colors[idx % colors.length];
      const centerX = b.x + (idx + 0.5) * slotWidth;

      const q1Y = yToPixel(s.q1);
      const q3Y = yToPixel(s.q3);
      const medY = yToPixel(s.median);
      const minFenceY = yToPixel(Math.max(s.min, s.lowerFence));
      const maxFenceY = yToPixel(Math.min(s.max, s.upperFence));

      // Whiskers (stem & end-caps)
      ctx.strokeStyle = pal.axis;
      ctx.lineWidth = 1.5;
      ctx.beginPath();
      // Lower whisker
      ctx.moveTo(centerX, q1Y);
      ctx.lineTo(centerX, minFenceY);
      ctx.moveTo(centerX - boxWidth * 0.25, minFenceY);
      ctx.lineTo(centerX + boxWidth * 0.25, minFenceY);
      // Upper whisker
      ctx.moveTo(centerX, q3Y);
      ctx.lineTo(centerX, maxFenceY);
      ctx.moveTo(centerX - boxWidth * 0.25, maxFenceY);
      ctx.lineTo(centerX + boxWidth * 0.25, maxFenceY);
      ctx.stroke();

      // Box body
      const boxTop = Math.min(q1Y, q3Y);
      const boxHeight = Math.abs(q1Y - q3Y) || 2;

      ctx.fillStyle = `${groupColor}22`; // 14% opacity
      ctx.strokeStyle = groupColor;
      ctx.lineWidth = 2;
      ctx.beginPath();
      ctx.roundRect?.(centerX - boxWidth / 2, boxTop, boxWidth, boxHeight, 4) ||
        ctx.rect(centerX - boxWidth / 2, boxTop, boxWidth, boxHeight);
      ctx.fill();
      ctx.stroke();

      // Median line
      ctx.strokeStyle = pal.textBold;
      ctx.lineWidth = 2.5;
      ctx.beginPath();
      ctx.moveTo(centerX - boxWidth / 2, medY);
      ctx.lineTo(centerX + boxWidth / 2, medY);
      ctx.stroke();

      // Mean Diamond Marker (◆)
      if (typeof s.mean === 'number' && isFinite(s.mean)) {
        const meanY = yToPixel(s.mean);
        ctx.fillStyle = pal.accent;
        ctx.strokeStyle = '#ffffff';
        ctx.lineWidth = 1.2;
        ctx.beginPath();
        const d = 5;
        ctx.moveTo(centerX, meanY - d);
        ctx.lineTo(centerX + d, meanY);
        ctx.lineTo(centerX, meanY + d);
        ctx.lineTo(centerX - d, meanY);
        ctx.closePath();
        ctx.fill();
        ctx.stroke();
      }

      // Outlier dots (distinguish mild vs extreme)
      if (s.outliers && s.outliers.length > 0) {
        for (const out of s.outliers) {
          const val = typeof out === 'object' && out !== null ? out.value : out;
          const isExtreme = typeof out === 'object' && out !== null ? out.type === 'Extreme' : false;
          const outY = yToPixel(val);

          ctx.fillStyle = isExtreme ? pal.danger : (pal.warning || '#f59e0b');
          ctx.strokeStyle = isExtreme ? '#ffffff' : pal.danger;
          ctx.lineWidth = 1.5;
          ctx.beginPath();
          ctx.arc(centerX, outY, isExtreme ? 5.5 : 4, 0, 2 * Math.PI);
          ctx.fill();
          ctx.stroke();

          // Callout for outlier
          ctx.fillStyle = isExtreme ? pal.danger : pal.text;
          ctx.font = `600 10px ${engine.options.fontFamily || 'sans-serif'}`;
          ctx.textAlign = 'left';
          ctx.fillText(`${val.toFixed(1)}${isExtreme ? ' (Extr)' : ''}`, centerX + boxWidth * 0.4, outY + 3);
        }
      }

      // Jittered raw points
      if (s.values && s.values.length <= 150) {
        ctx.fillStyle = `${groupColor}88`;
        for (let i = 0; i < s.values.length; i++) {
          const v = s.values[i];
          const ptY = yToPixel(v);
          const jitter = (Math.sin(i * 12.9898 + v) * 0.35) * (boxWidth / 2);
          ctx.beginPath();
          ctx.arc(centerX + jitter, ptY, 3, 0, 2 * Math.PI);
          ctx.fill();
        }
      }

      // Single cohort side annotations
      if (k === 1) {
        ctx.font = `500 10px ${engine.options.fontFamily || 'sans-serif'}`;
        ctx.textAlign = 'left';
        ctx.fillStyle = pal.textBold;
        ctx.fillText(`Median: ${s.median.toFixed(1)}`, centerX + boxWidth / 2 + 8, medY + 3);

        ctx.fillStyle = pal.text;
        ctx.fillText(`Q3: ${s.q3.toFixed(1)}`, centerX + boxWidth / 2 + 8, q3Y + 3);
        ctx.fillText(`Q1: ${s.q1.toFixed(1)}`, centerX + boxWidth / 2 + 8, q1Y + 3);

        if (typeof s.mean === 'number') {
          ctx.fillStyle = pal.accent;
          ctx.fillText(`Mean: ${s.mean.toFixed(1)}`, centerX - boxWidth / 2 - 62, yToPixel(s.mean) + 3);
        }
      }
    });
  },

  /**
   * Renders Comparative Cohort Error Bar / Dispersion Plot
   * Supports:
   *  - 'ci95': 95% Confidence Interval (Mean ± t* x SEM)
   *  - 'sem': Standard Error of Mean (Mean ± 1 SEM)
   *  - 'sd': Standard Deviation (Mean ± 1 SD)
   *  - 'iqr': Interquartile Range (Box & Whiskers / Median ± IQR via renderBoxPlot)
   * @param {ChartEngine} engine
   * @param {Array<{name: string, stats: object, color?: string}>} groupStats
   * @param {object|string} [options={}] - Options object { mode: 'ci95'|'sem'|'sd'|'iqr', title?: string } or title string
   */
  renderErrorBarPlot(engine, groupStats, options = {}) {
    const opts = typeof options === 'string' ? { title: options } : (options || {});
    const mode = (opts.mode || 'ci95').toLowerCase();
    const title = opts.title || 'Comparative Cohort Distribution';

    // If IQR is requested, delegate directly to the Box-and-Whisker plot
    if (mode === 'iqr') {
      return this.renderBoxPlot(engine, groupStats, title);
    }

    engine.lastRenderFn = () => this.renderErrorBarPlot(engine, groupStats, options);
    engine.lastRender = engine.lastRenderFn;
    engine.clear();

    const b = engine.getPlotBounds ? engine.getPlotBounds() : engine.getBounds();
    const ctx = engine.ctx;
    const pal = engine.palette;

    let groups = [];
    if (Array.isArray(groupStats)) {
      groups = groupStats.map(g => g.stats ? g : { name: g.name || 'Sample Data', stats: g, color: g.color });
    } else if (groupStats && typeof groupStats === 'object') {
      groups = [{ name: groupStats.name || 'Sample Data', stats: groupStats, color: pal.primary }];
    }

    if (groups.length === 0) return;

    // Determine error bounds for each group based on chosen mode
    const processedGroups = groups.map((g, idx) => {
      const s = g.stats || {};
      const mean = typeof s.mean === 'number' && isFinite(s.mean) ? s.mean : 0;
      const sem = typeof s.sem === 'number' && isFinite(s.sem) ? s.sem : (s.sd && s.n ? s.sd / Math.sqrt(s.n) : 0);
      const sd = typeof s.sd === 'number' && isFinite(s.sd) ? s.sd : 0;

      let lower = mean;
      let upper = mean;
      let label = '95% CI';
      let subLabel = '';

      if (mode === 'sem') {
        lower = mean - sem;
        upper = mean + sem;
        label = '±1 SEM';
        subLabel = `SEM: ±${sem.toFixed(2)}`;
      } else if (mode === 'sd') {
        lower = mean - sd;
        upper = mean + sd;
        label = '±1 SD';
        subLabel = `SD: ±${sd.toFixed(2)}`;
      } else {
        // 'ci95' default
        if (Array.isArray(s.ci95) && s.ci95.length === 2 && isFinite(s.ci95[0]) && isFinite(s.ci95[1])) {
          lower = s.ci95[0];
          upper = s.ci95[1];
        } else {
          const margin = sem * 1.96;
          lower = mean - margin;
          upper = mean + margin;
        }
        label = '95% CI';
        subLabel = `[${lower.toFixed(2)}, ${upper.toFixed(2)}]`;
      }

      return {
        ...g,
        mean,
        sem,
        sd,
        lower,
        upper,
        label,
        subLabel
      };
    });

    // Find global min and max across all groups, including raw points and error bounds
    let globalMin = Infinity;
    let globalMax = -Infinity;
    for (const g of processedGroups) {
      const s = g.stats;
      if (!s || s.n === 0) continue;
      if (g.lower < globalMin) globalMin = g.lower;
      if (g.upper > globalMax) globalMax = g.upper;
      if (typeof s.min === 'number' && isFinite(s.min) && s.min < globalMin) globalMin = s.min;
      if (typeof s.max === 'number' && isFinite(s.max) && s.max > globalMax) globalMax = s.max;
    }

    if (!isFinite(globalMin) || !isFinite(globalMax)) return;

    // Extra margin so error bar caps and callouts have breathing room
    const span = globalMax - globalMin || 1;
    const yMin = globalMin - 0.15 * span;
    const yMax = globalMax + 0.18 * span;
    const yRange = yMax - yMin;

    const yToPixel = (val) => b.y + b.height - ((val - yMin) / yRange) * b.height;

    // Build Y ticks
    const yTicks = [];
    const stepCount = 5;
    for (let i = 0; i <= stepCount; i++) {
      const v = yMin + (i / stepCount) * yRange;
      yTicks.push({ norm: (v - yMin) / yRange, label: v.toFixed(1) });
    }

    // Build X ticks
    const k = processedGroups.length;
    const xTicks = processedGroups.map((g, idx) => ({
      norm: (idx + 0.5) / k,
      label: g.name
    }));

    const modeHeaders = {
      ci95: '95% Confidence Interval (Mean ± 95% CI)',
      sem: 'Standard Error of Mean (Mean ± 1 SEM)',
      sd: 'Standard Deviation (Mean ± 1 SD)',
      iqr: 'Interquartile Range'
    };
    const modeTitle = modeHeaders[mode] || '95% Confidence Interval';

    engine.drawAxes({
      yTicks,
      xTicks,
      title: `${title}`,
      yLabel: 'Observed Value'
    });

    const slotWidth = b.width / k;
    const capWidth = Math.min(38, Math.max(22, slotWidth * 0.28));
    const colors = [pal.primary, pal.secondary, pal.accent, '#f59e0b'];

    // Draw connecting delta bridge between 2 groups if k === 2
    if (k === 2 && processedGroups[0].stats && processedGroups[1].stats) {
      const gA = processedGroups[0];
      const gB = processedGroups[1];
      const xA = b.x + 0.5 * slotWidth;
      const xB = b.x + 1.5 * slotWidth;
      const yA = yToPixel(gA.mean);
      const yB = yToPixel(gB.mean);

      ctx.save();
      ctx.strokeStyle = `${pal.axis}`;
      ctx.setLineDash([4, 4]);
      ctx.lineWidth = 1.5;
      ctx.beginPath();
      ctx.moveTo(xA, yA);
      ctx.lineTo(xB, yB);
      ctx.stroke();
      ctx.setLineDash([]);

      // Delta label pill in middle of bridge
      const midX = (xA + xB) / 2;
      const midY = (yA + yB) / 2;
      const deltaM = gA.mean - gB.mean;
      const deltaText = `ΔM = ${deltaM >= 0 ? '+' : ''}${deltaM.toFixed(2)}`;

      ctx.font = `600 10px ${engine.options.fontFamily || 'sans-serif'}`;
      const textW = ctx.measureText ? (ctx.measureText(deltaText)?.width || 50) : 50;
      ctx.fillStyle = pal.bgSurfaceElevated || pal.bgCard || '#1e293b';
      ctx.strokeStyle = pal.border || '#334155';
      ctx.lineWidth = 1;
      ctx.beginPath();
      ctx.roundRect?.(midX - textW / 2 - 6, midY - 10, textW + 12, 18, 4) ||
        ctx.rect(midX - textW / 2 - 6, midY - 10, textW + 12, 18);
      ctx.fill();
      ctx.stroke();

      ctx.fillStyle = pal.textBold;
      ctx.textAlign = 'center';
      ctx.textBaseline = 'middle';
      ctx.fillText(deltaText, midX, midY);
      ctx.restore();
    }

    // Render each group
    processedGroups.forEach((g, idx) => {
      const s = g.stats;
      if (!s || s.n === 0) return;

      const groupColor = g.color || colors[idx % colors.length];
      const centerX = b.x + (idx + 0.5) * slotWidth;

      const meanY = yToPixel(g.mean);
      const lowerY = yToPixel(g.lower);
      const upperY = yToPixel(g.upper);

      // 1. Raw points jittered in background
      if (s.values && s.values.length > 0) {
        ctx.save();
        ctx.fillStyle = `${groupColor}44`;
        for (let i = 0; i < s.values.length; i++) {
          const v = s.values[i];
          const ptY = yToPixel(v);
          const jitter = (Math.sin(i * 12.9898 + v) * 0.35) * (slotWidth * 0.3);
          ctx.beginPath();
          ctx.arc(centerX + jitter, ptY, 3.2, 0, 2 * Math.PI);
          ctx.fill();
        }
        ctx.restore();
      }

      // 2. Translucent pillar / shaded dispersion band
      const pillarWidth = Math.min(52, slotWidth * 0.35);
      const pillarTop = Math.min(lowerY, upperY);
      const pillarHeight = Math.abs(lowerY - upperY) || 2;
      ctx.save();
      ctx.fillStyle = `${groupColor}14`;
      ctx.beginPath();
      ctx.roundRect?.(centerX - pillarWidth / 2, pillarTop, pillarWidth, pillarHeight, 6) ||
        ctx.rect(centerX - pillarWidth / 2, pillarTop, pillarWidth, pillarHeight);
      ctx.fill();
      ctx.restore();

      // 3. Error bar vertical stem
      ctx.save();
      ctx.strokeStyle = groupColor;
      ctx.lineWidth = 2.5;
      ctx.beginPath();
      ctx.moveTo(centerX, lowerY);
      ctx.lineTo(centerX, upperY);
      ctx.stroke();

      // 4. Horizontal serif end-caps
      ctx.beginPath();
      // Upper cap
      ctx.moveTo(centerX - capWidth / 2, upperY);
      ctx.lineTo(centerX + capWidth / 2, upperY);
      // Lower cap
      ctx.moveTo(centerX - capWidth / 2, lowerY);
      ctx.lineTo(centerX + capWidth / 2, lowerY);
      ctx.stroke();
      ctx.restore();

      // 5. Mean Point Marker (Prominent circle with halo)
      ctx.save();
      ctx.fillStyle = groupColor;
      ctx.strokeStyle = '#ffffff';
      ctx.lineWidth = 2.5;
      ctx.beginPath();
      ctx.arc(centerX, meanY, 6.5, 0, 2 * Math.PI);
      ctx.fill();
      ctx.stroke();

      // Inner center dot
      ctx.fillStyle = '#ffffff';
      ctx.beginPath();
      ctx.arc(centerX, meanY, 2, 0, 2 * Math.PI);
      ctx.fill();
      ctx.restore();

      // 6. Text callout annotations
      ctx.save();
      ctx.font = `600 11px ${engine.options.fontFamily || 'sans-serif'}`;
      ctx.textAlign = 'left';
      ctx.textBaseline = 'middle';

      // Mean callout
      ctx.fillStyle = pal.textBold;
      ctx.fillText(`M = ${g.mean.toFixed(2)}`, centerX + capWidth / 2 + 8, meanY);

      // Sub-label for dispersion / error metric
      ctx.font = `600 10px ${engine.options.fontFamily || 'sans-serif'}`;
      ctx.fillStyle = pal.textMuted || pal.text;
      ctx.fillText(g.subLabel, centerX + capWidth / 2 + 8, meanY + 13);

      // Explicit error bar height / total span
      ctx.font = `500 9px ${engine.options.fontFamily || 'sans-serif'}`;
      ctx.fillStyle = pal.textDim || pal.textMuted || pal.text;
      ctx.fillText(`Span: ${(g.upper - g.lower).toFixed(2)}`, centerX + capWidth / 2 + 8, meanY + 25);

      // Sample size n
      ctx.fillText(`n = ${s.n}`, centerX + capWidth / 2 + 8, meanY + 36);

      // Top cap value
      ctx.textAlign = 'right';
      ctx.font = `500 9px ${engine.options.fontFamily || 'sans-serif'}`;
      ctx.fillStyle = pal.textDim || pal.textMuted || pal.text;
      ctx.fillText(g.upper.toFixed(2), centerX - capWidth / 2 - 6, upperY);

      // Bottom cap value
      ctx.fillText(g.lower.toFixed(2), centerX - capWidth / 2 - 6, lowerY);

      ctx.restore();
    });

    // 7. Header note in top-right corner
    ctx.save();
    ctx.font = `italic 10px ${engine.options.fontFamily || 'sans-serif'}`;
    ctx.textAlign = 'right';
    ctx.textBaseline = 'top';
    ctx.fillStyle = pal.textMuted || pal.text;
    ctx.fillText(`Mode: ${modeTitle}`, b.x + b.width, b.y - 18);
    ctx.restore();
  },

  /**
   * Renders Violin Density Plot combining Kernel Density Estimation (KDE) with internal Box/Quartile markers
   * @param {ChartEngine} engine
   * @param {object|Array} data - Stats object or array of group stats
   * @param {string} title
   */
  renderViolinPlot(engine, data, title = 'Violin Density Plot (KDE & Quartiles)') {
    engine.lastRenderFn = () => this.renderViolinPlot(engine, data, title);
    engine.clear();
    const b = engine.getPlotBounds ? engine.getPlotBounds() : engine.getBounds();
    const ctx = engine.ctx;
    const pal = engine.palette;

    let groups = [];
    if (Array.isArray(data)) {
      groups = data.map(g => g.stats ? g : { name: g.name || 'Sample Data', stats: g, color: g.color });
    } else if (data && typeof data === 'object') {
      groups = [{ name: data.name || 'Sample Data', stats: data, color: pal.secondary }];
    }

    if (groups.length === 0) return;

    let globalMin = Infinity;
    let globalMax = -Infinity;
    for (const g of groups) {
      const s = g.stats;
      if (!s || s.n === 0) continue;
      if (s.min < globalMin) globalMin = s.min;
      if (s.max > globalMax) globalMax = s.max;
    }

    if (!isFinite(globalMin) || !isFinite(globalMax)) return;

    const span = globalMax - globalMin || 1;
    const yMin = globalMin - 0.12 * span;
    const yMax = globalMax + 0.12 * span;
    const yRange = yMax - yMin;

    const yToPixel = (val) => b.y + b.height - ((val - yMin) / yRange) * b.height;

    // Build Y ticks
    const yTicks = [];
    const stepCount = 5;
    for (let i = 0; i <= stepCount; i++) {
      const v = yMin + (i / stepCount) * yRange;
      yTicks.push({ norm: (v - yMin) / yRange, label: v.toFixed(1) });
    }

    // Build X ticks
    const k = groups.length;
    const xTicks = groups.map((g, idx) => ({
      norm: (idx + 0.5) / k,
      label: g.name
    }));

    engine.drawAxes({ yTicks, xTicks, title, yLabel: 'Observed Value' });

    const slotWidth = b.width / k;
    const maxHalfWidth = Math.min(75, slotWidth * 0.42);

    groups.forEach((g, idx) => {
      const s = g.stats;
      if (!s || !s.values || s.values.length < 2) return;

      const groupColor = g.color || pal.secondary;
      const centerX = b.x + (idx + 0.5) * slotWidth;
      const values = s.values;
      const n = values.length;

      // Silverman rule of thumb bandwidth
      const iqr = s.iqr > 0 ? s.iqr : (s.sd || 1);
      const A = Math.min(s.sd || 1, (iqr / 1.34) || 1) || 1;
      const h = Math.max(1e-3 * span, 0.9 * A * Math.pow(n, -0.2));

      // Gaussian KDE function
      const kde = (y) => {
        let sum = 0;
        const invH = 1 / h;
        for (let i = 0; i < n; i++) {
          const u = (y - values[i]) * invH;
          sum += Math.exp(-0.5 * u * u);
        }
        return sum / (n * h * Math.sqrt(2 * Math.PI));
      };

      // Evaluate density on grid
      const gridSteps = 100;
      const gridY = [];
      const density = [];
      let maxDensity = 0;

      for (let i = 0; i <= gridSteps; i++) {
        const yVal = s.min + (i / gridSteps) * (s.max - s.min);
        const d = kde(yVal);
        gridY.push(yVal);
        density.push(d);
        if (d > maxDensity) maxDensity = d;
      }

      if (maxDensity <= 0) maxDensity = 1;

      // Draw mirrored violin silhouette
      ctx.beginPath();
      // Right side (ascending Y)
      for (let i = 0; i <= gridSteps; i++) {
        const py = yToPixel(gridY[i]);
        const px = centerX + (density[i] / maxDensity) * maxHalfWidth;
        if (i === 0) ctx.moveTo(centerX, py);
        else ctx.lineTo(px, py);
      }
      ctx.lineTo(centerX, yToPixel(gridY[gridSteps]));

      // Left side (descending Y)
      for (let i = gridSteps; i >= 0; i--) {
        const py = yToPixel(gridY[i]);
        const px = centerX - (density[i] / maxDensity) * maxHalfWidth;
        ctx.lineTo(px, py);
      }
      ctx.closePath();

      // Translucent gradient fill
      const grad = ctx.createLinearGradient(centerX - maxHalfWidth, 0, centerX + maxHalfWidth, 0);
      grad.addColorStop(0, `${groupColor}15`);
      grad.addColorStop(0.5, `${groupColor}44`);
      grad.addColorStop(1, `${groupColor}15`);
      ctx.fillStyle = grad;
      ctx.fill();

      ctx.strokeStyle = groupColor;
      ctx.lineWidth = 2;
      ctx.stroke();

      // Range whisker stem (min to max along center)
      ctx.strokeStyle = pal.axis;
      ctx.lineWidth = 2;
      ctx.beginPath();
      ctx.moveTo(centerX, yToPixel(s.min));
      ctx.lineTo(centerX, yToPixel(s.max));
      ctx.stroke();

      // Internal IQR Box (Q1 to Q3)
      const q1Y = yToPixel(s.q1);
      const q3Y = yToPixel(s.q3);
      const iqrBoxWidth = 10;
      ctx.fillStyle = pal.bg === '#ffffff' ? '#e2e8f0' : '#1e293b';
      ctx.strokeStyle = pal.textBold;
      ctx.lineWidth = 2;
      const iqrTop = Math.min(q1Y, q3Y);
      const iqrHeight = Math.abs(q1Y - q3Y) || 2;
      ctx.beginPath();
      ctx.roundRect?.(centerX - iqrBoxWidth / 2, iqrTop, iqrBoxWidth, iqrHeight, 3) ||
        ctx.rect(centerX - iqrBoxWidth / 2, iqrTop, iqrBoxWidth, iqrHeight);
      ctx.fill();
      ctx.stroke();

      // White solid dot for Median
      const medY = yToPixel(s.median);
      ctx.fillStyle = '#ffffff';
      ctx.strokeStyle = pal.textBold;
      ctx.lineWidth = 1.5;
      ctx.beginPath();
      ctx.arc(centerX, medY, 4, 0, 2 * Math.PI);
      ctx.fill();
      ctx.stroke();

      // Mean cyan diamond marker (◆)
      if (typeof s.mean === 'number' && isFinite(s.mean)) {
        const meanY = yToPixel(s.mean);
        ctx.fillStyle = pal.primary;
        ctx.strokeStyle = '#ffffff';
        ctx.lineWidth = 1.2;
        ctx.beginPath();
        const mSize = 4.5;
        ctx.moveTo(centerX, meanY - mSize);
        ctx.lineTo(centerX + mSize, meanY);
        ctx.lineTo(centerX, meanY + mSize);
        ctx.lineTo(centerX - mSize, meanY);
        ctx.closePath();
        ctx.fill();
        ctx.stroke();
      }

      // Single-cohort annotations
      if (k === 1) {
        ctx.font = `600 10px ${engine.options.fontFamily || 'sans-serif'}`;
        ctx.textAlign = 'left';
        ctx.fillStyle = pal.textBold;
        ctx.fillText(`Median: ${s.median.toFixed(1)}`, centerX + maxHalfWidth * 0.7 + 10, medY + 3);

        if (typeof s.mean === 'number') {
          ctx.fillStyle = pal.primary;
          ctx.fillText(`Mean: ${s.mean.toFixed(1)}`, centerX - maxHalfWidth * 0.7 - 65, yToPixel(s.mean) + 3);
        }

        ctx.fillStyle = pal.text;
        ctx.fillText(`IQR [${s.q1.toFixed(1)}, ${s.q3.toFixed(1)}]`, centerX + maxHalfWidth * 0.7 + 10, q3Y + 3);
      }
    });
  },

  /**
   * Renders Histogram with fitted Normal curve
   * @param {ChartEngine} engine
   * @param {object} stats - from Descriptive.calculate
   * @param {string} title
   */
  renderHistogram(engine, stats, title = 'Frequency Distribution & Normal Fit') {
    engine.lastRenderFn = () => this.renderHistogram(engine, stats, title);
    engine.clear();
    const b = engine.getPlotBounds();
    const ctx = engine.ctx;
    const pal = engine.palette;

    if (!stats || stats.n < 3) return;

    const values = stats.values;
    const min = stats.min;
    const max = stats.max;
    const span = max - min || 1;

    // Number of bins via Freedman-Diaconis or Sturges
    const numBins = Math.max(5, Math.min(25, Math.round(1 + 3.322 * Math.log10(stats.n))));
    const binWidth = span / numBins;

    const bins = new Array(numBins).fill(0);
    for (const v of values) {
      const idx = Math.min(numBins - 1, Math.floor((v - min) / binWidth));
      bins[idx]++;
    }

    const maxCount = Math.max(...bins, 1);
    const yMax = maxCount * 1.25;

    // Ticks
    const yTicks = [
      { norm: 0, label: '0' },
      { norm: 0.5, label: (yMax * 0.5).toFixed(0) },
      { norm: 1.0, label: yMax.toFixed(0) }
    ];

    const xTicks = [];
    for (let i = 0; i <= 4; i++) {
      const v = min + (i / 4) * span;
      xTicks.push({ norm: i / 4, label: v.toFixed(1) });
    }

    engine.drawAxes({ yTicks, xTicks, title, xLabel: 'Observation Value', yLabel: 'Frequency' });

    // Draw Histogram bars
    const barWidth = b.width / numBins;
    for (let i = 0; i < numBins; i++) {
      const count = bins[i];
      const h = (count / yMax) * b.height;
      const x = b.x + i * barWidth;
      const y = b.y + b.height - h;

      ctx.fillStyle = `${pal.primary}33`;
      ctx.strokeStyle = pal.primary;
      ctx.lineWidth = 1.5;
      ctx.beginPath();
      ctx.rect(x + 1, y, barWidth - 2, h);
      ctx.fill();
      ctx.stroke();
    }

    // Fitted Normal Distribution overlay curve
    if (stats.sd > 0) {
      ctx.strokeStyle = pal.accent;
      ctx.lineWidth = 2.5;
      ctx.beginPath();

      const points = 100;
      for (let i = 0; i <= points; i++) {
        const xVal = min + (i / points) * span;
        // Normal density
        const z = (xVal - stats.mean) / stats.sd;
        const density = (1 / (stats.sd * Math.sqrt(2 * Math.PI))) * Math.exp(-0.5 * z * z);
        // Expected frequency count in bin
        const expectedCount = density * binWidth * stats.n;
        const normY = expectedCount / yMax;
        const px = b.x + (i / points) * b.width;
        const py = b.y + b.height - Math.min(b.height, normY * b.height);

        if (i === 0) ctx.moveTo(px, py);
        else ctx.lineTo(px, py);
      }
      ctx.stroke();
    }
  },

  /**
   * Renders Scatter plot with OLS linear regression line and confidence bands
   * @param {ChartEngine} engine
   * @param {Array<{x: number, y: number}>} pairs
   * @param {object} reg - from Correlation.linearRegression
   * @param {string} title
   */
  renderScatterRegression(engine, pairs, reg, title = 'Linear Regression & Correlation', quad = null) {
    engine.lastRenderFn = () => this.renderScatterRegression(engine, pairs, reg, title, quad);
    engine.clear();
    const b = engine.getPlotBounds();
    const ctx = engine.ctx;
    const pal = engine.palette;

    if (!pairs || pairs.length < 2) return;

    const xVals = pairs.map(p => p.x);
    const yVals = pairs.map(p => p.y);

    const minX = Math.min(...xVals);
    const maxX = Math.max(...xVals);
    const minY = Math.min(...yVals);
    const maxY = Math.max(...yVals);

    const spanX = maxX - minX || 1;
    const spanY = maxY - minY || 1;

    const x0 = minX - 0.08 * spanX;
    const x1 = maxX + 0.08 * spanX;
    const y0 = minY - 0.08 * spanY;
    const y1 = maxY + 0.08 * spanY;

    const toX = (val) => b.x + ((val - x0) / (x1 - x0)) * b.width;
    const toY = (val) => b.y + b.height - ((val - y0) / (y1 - y0)) * b.height;

    const xTicks = [];
    const yTicks = [];
    for (let i = 0; i <= 4; i++) {
      const vx = x0 + (i / 4) * (x1 - x0);
      const vy = y0 + (i / 4) * (y1 - y0);
      xTicks.push({ norm: i / 4, label: vx.toFixed(1) });
      yTicks.push({ norm: i / 4, label: vy.toFixed(1) });
    }

    engine.drawAxes({ xTicks, yTicks, title, xLabel: 'Independent Variable (X)', yLabel: 'Dependent Variable (Y)' });

    const isQuadActive = quad && quad.isSignificantlyQuadratic && (quad.shape === 'U-Shaped' || quad.shape === 'Inverted U-Shaped');

    // Draw 95% Confidence Band for regression line if available
    if (reg && reg.seResidual && reg.sxx > 0 && !isQuadActive) {
      ctx.fillStyle = `${pal.primary}18`; // light tint
      ctx.beginPath();

      const upperBand = [];
      const lowerBand = [];
      const steps = 60;
      const tCrit = 1.96;

      for (let i = 0; i <= steps; i++) {
        const vx = x0 + (i / steps) * (x1 - x0);
        const yPred = reg.intercept + reg.slope * vx;
        const seFit = reg.seResidual * Math.sqrt(1 / reg.n + Math.pow(vx - reg.xStats.mean, 2) / reg.sxx);
        upperBand.push({ px: toX(vx), py: toY(yPred + tCrit * seFit) });
        lowerBand.push({ px: toX(vx), py: toY(yPred - tCrit * seFit) });
      }

      ctx.moveTo(upperBand[0].px, upperBand[0].py);
      for (const pt of upperBand) ctx.lineTo(pt.px, pt.py);
      for (let i = lowerBand.length - 1; i >= 0; i--) ctx.lineTo(lowerBand[i].px, lowerBand[i].py);
      ctx.closePath();
      ctx.fill();
    }

    // Linear Regression Line
    if (reg && reg.slope !== undefined) {
      ctx.strokeStyle = isQuadActive ? `${pal.primary}77` : pal.primary;
      ctx.lineWidth = isQuadActive ? 1.8 : 2.5;
      if (isQuadActive) ctx.setLineDash([4, 4]);
      ctx.beginPath();
      ctx.moveTo(toX(x0), toY(reg.intercept + reg.slope * x0));
      ctx.lineTo(toX(x1), toY(reg.intercept + reg.slope * x1));
      ctx.stroke();
      ctx.setLineDash([]);
    }

    // Quadratic fit curve & Nadir/Zenith point if U-shaped
    if (isQuadActive) {
      ctx.strokeStyle = '#f59e0b'; // amber-500
      ctx.lineWidth = 2.5;
      ctx.beginPath();
      const steps = 100;
      for (let i = 0; i <= steps; i++) {
        const vx = x0 + (i / steps) * (x1 - x0);
        const vy = quad.b0 + quad.b1 * vx + quad.b2 * vx * vx;
        const px = toX(vx);
        const py = toY(vy);
        if (i === 0) ctx.moveTo(px, py);
        else ctx.lineTo(px, py);
      }
      ctx.stroke();

      // Plot Nadir or Zenith point
      if (quad.vertexX >= minX - spanX * 0.1 && quad.vertexX <= maxX + spanX * 0.1) {
        const vx = quad.vertexX;
        const vy = quad.vertexY;
        const px = toX(vx);
        const py = toY(vy);

        // Glowing aura
        ctx.fillStyle = 'rgba(245, 158, 11, 0.25)';
        ctx.beginPath();
        ctx.arc(px, py, 9, 0, 2 * Math.PI);
        ctx.fill();

        // Core vertex dot
        ctx.fillStyle = '#f59e0b';
        ctx.strokeStyle = '#ffffff';
        ctx.lineWidth = 1.5;
        ctx.beginPath();
        ctx.arc(px, py, 5, 0, 2 * Math.PI);
        ctx.fill();
        ctx.stroke();

        // Vertex text badge
        ctx.fillStyle = '#f59e0b';
        ctx.font = '600 11px Inter, system-ui, sans-serif';
        const label = `${quad.b2 > 0 ? 'Nadir' : 'Zenith'}: (${vx.toFixed(1)}, ${vy.toFixed(1)})`;
        ctx.fillText(label, px + 8, py - 6);
      }

      // Mini legend
      ctx.font = '10px Inter, system-ui, sans-serif';
      ctx.fillStyle = pal.primary;
      ctx.fillRect(b.x + b.width - 165, b.y + 12, 12, 3);
      ctx.fillStyle = pal.textMuted || '#94a3b8';
      ctx.fillText(`Linear OLS (r = ${(reg.r !== undefined ? reg.r : Math.sqrt(Math.max(0, reg.rSquared || 0))).toFixed(2)})`, b.x + b.width - 148, b.y + 16);

      ctx.fillStyle = '#f59e0b';
      ctx.fillRect(b.x + b.width - 165, b.y + 26, 12, 3);
      ctx.fillStyle = '#f59e0b';
      ctx.fillText(`Quadratic (${quad.shape}, R² = ${quad.rSquaredQuad.toFixed(2)})`, b.x + b.width - 148, b.y + 30);
    }

    // Draw Scatter points
    ctx.fillStyle = pal.secondary;
    ctx.strokeStyle = '#ffffff';
    ctx.lineWidth = 1;
    for (const p of pairs) {
      ctx.beginPath();
      ctx.arc(toX(p.x), toY(p.y), 4.5, 0, 2 * Math.PI);
      ctx.fill();
      ctx.stroke();
    }
  },

  /**
   * Renders ROC curve with AUC shading and diagonal reference line
   * @param {ChartEngine} engine
   * @param {object} roc - from Diagnostic.computeROC
   * @param {string} title
   */
  renderROC(engine, roc, title = 'Receiver Operating Characteristic (ROC)') {
    engine.lastRenderFn = () => this.renderROC(engine, roc, title);
    engine.clear();
    const b = engine.getPlotBounds();
    const ctx = engine.ctx;
    const pal = engine.palette;

    if (!roc || !roc.points || roc.points.length === 0) return;

    // 0 to 1 on both axes
    const xTicks = [
      { norm: 0, label: '0.0' },
      { norm: 0.25, label: '0.25' },
      { norm: 0.5, label: '0.50' },
      { norm: 0.75, label: '0.75' },
      { norm: 1.0, label: '1.0' }
    ];
    const yTicks = [
      { norm: 0, label: '0.0' },
      { norm: 0.25, label: '0.25' },
      { norm: 0.5, label: '0.50' },
      { norm: 0.75, label: '0.75' },
      { norm: 1.0, label: '1.0' }
    ];

    engine.drawAxes({
      xTicks,
      yTicks,
      title,
      xLabel: '1 - Specificity (False Positive Rate)',
      yLabel: 'Sensitivity (True Positive Rate)'
    });

    const toX = (fpr) => b.x + Math.min(1, Math.max(0, fpr)) * b.width;
    const toY = (tpr) => b.y + b.height - Math.min(1, Math.max(0, tpr)) * b.height;

    // Diagonal Chance Line (dashed)
    ctx.strokeStyle = pal.axis;
    ctx.lineWidth = 1.5;
    ctx.setLineDash([5, 5]);
    ctx.beginPath();
    ctx.moveTo(b.x, b.y + b.height);
    ctx.lineTo(b.x + b.width, b.y);
    ctx.stroke();
    ctx.setLineDash([]);

    // Area Under Curve Shading
    ctx.fillStyle = `${pal.primary}22`;
    ctx.beginPath();
    ctx.moveTo(b.x, b.y + b.height);
    for (const pt of roc.points) {
      ctx.lineTo(toX(pt.fpr), toY(pt.tpr));
    }
    ctx.lineTo(b.x + b.width, b.y + b.height);
    ctx.closePath();
    ctx.fill();

    // ROC Line
    ctx.strokeStyle = pal.primary;
    ctx.lineWidth = 3;
    ctx.beginPath();
    for (let i = 0; i < roc.points.length; i++) {
      const pt = roc.points[i];
      if (i === 0) ctx.moveTo(toX(pt.fpr), toY(pt.tpr));
      else ctx.lineTo(toX(pt.fpr), toY(pt.tpr));
    }
    ctx.stroke();

    // Optimal Cutoff marker (Youden's J)
    if (roc.optimalCutoff) {
      const bestFPR = 1 - roc.optimalCutoff.spec;
      const bestTPR = roc.optimalCutoff.sens;
      const optX = toX(bestFPR);
      const optY = toY(bestTPR);

      ctx.fillStyle = pal.danger;
      ctx.strokeStyle = '#ffffff';
      ctx.lineWidth = 2;
      ctx.beginPath();
      ctx.arc(optX, optY, 6, 0, 2 * Math.PI);
      ctx.fill();
      ctx.stroke();

      ctx.fillStyle = pal.textBold;
      ctx.font = `600 11px ${engine.options.fontFamily}`;
      ctx.fillText(
        `Optimal J = ${roc.optimalCutoff.youdenJ.toFixed(3)} (Cutoff: ${roc.optimalCutoff.threshold})`,
        optX + 10,
        optY - 8
      );
    }
  },

  /**
   * Renders Teaching Distribution: Empirical Histogram + Theoretical PDF Overlay Curve
   */
  renderTeachingDistribution(engine, data, distKey, params = {}, title = 'Generated Distribution & Theoretical PDF') {
    engine.lastRenderFn = () => this.renderTeachingDistribution(engine, data, distKey, params, title);
    engine.clear();
    const b = engine.getPlotBounds();
    const ctx = engine.ctx;
    const pal = engine.palette;

    if (!data || data.length < 5) {
      ctx.fillStyle = pal.textDim || '#94a3b8';
      ctx.font = `14px ${engine.options.fontFamily}`;
      ctx.textAlign = 'center';
      ctx.fillText('Click "Generate New Sample" to simulate distribution data.', b.x + b.width / 2, b.y + b.height / 2);
      return;
    }

    const n = data.length;
    let minVal = Infinity;
    let maxVal = -Infinity;
    let sum = 0;
    for (let i = 0; i < n; i++) {
      const v = data[i];
      if (v < minVal) minVal = v;
      if (v > maxVal) maxVal = v;
      sum += v;
    }
    const sampleMean = sum / n;

    const span = maxVal - minVal || 1;
    const plotMin = minVal - 0.05 * span;
    const plotMax = maxVal + 0.05 * span;
    const plotSpan = plotMax - plotMin;

    // Number of bins via Freedman-Diaconis or Sturges
    const numBins = Math.max(12, Math.min(35, Math.round(1 + 3.322 * Math.log10(n) * 1.5)));
    const binWidth = plotSpan / numBins;

    const bins = new Array(numBins).fill(0);
    for (let i = 0; i < n; i++) {
      const idx = Math.min(numBins - 1, Math.max(0, Math.floor((data[i] - plotMin) / binWidth)));
      bins[idx]++;
    }

    const maxCount = Math.max(...bins, 1);
    const yMax = maxCount * 1.25;

    // Build Axes Ticks
    const yTicks = [
      { norm: 0, label: '0' },
      { norm: 0.5, label: (yMax * 0.5).toFixed(0) },
      { norm: 1.0, label: yMax.toFixed(0) }
    ];

    const xTicks = [];
    for (let i = 0; i <= 5; i++) {
      const v = plotMin + (i / 5) * plotSpan;
      xTicks.push({ norm: i / 5, label: v.toFixed(1) });
    }

    engine.drawAxes({
      yTicks,
      xTicks,
      title,
      xLabel: 'Observation Value (X)',
      yLabel: 'Frequency Count'
    });

    const toX = (val) => b.x + ((val - plotMin) / plotSpan) * b.width;
    const toY = (count) => b.y + b.height - (count / yMax) * b.height;

    // 1. Draw Histogram Bars
    const barWidth = b.width / numBins;
    for (let i = 0; i < numBins; i++) {
      const count = bins[i];
      if (count === 0) continue;
      const h = (count / yMax) * b.height;
      const x = b.x + i * barWidth;
      const y = b.y + b.height - h;

      ctx.fillStyle = `${pal.primary}33`;
      ctx.strokeStyle = pal.primary;
      ctx.lineWidth = 1.5;
      ctx.beginPath();
      ctx.rect(x + 1, y, Math.max(1, barWidth - 2), h);
      ctx.fill();
      ctx.stroke();
    }

    // 2. Draw Theoretical PDF Overlay Curve
    if (Teaching.pdf && Teaching.pdf[distKey]) {
      const pdfFn = Teaching.pdf[distKey];
      const densityScale = n * binWidth; // Converts density f(x) to frequency count

      ctx.strokeStyle = pal.accent || '#38bdf8';
      ctx.lineWidth = 2.8;
      ctx.beginPath();

      const steps = 180;
      let started = false;
      for (let s = 0; s <= steps; s++) {
        const xVal = plotMin + (s / steps) * plotSpan;
        let pdfVal = 0;

        switch (distKey) {
          case 'normal':
            pdfVal = pdfFn.call(Teaching.pdf, xVal, params.mean, params.sd);
            break;
          case 'studentsT':
            pdfVal = pdfFn.call(Teaching.pdf, xVal, params.df, params.mean, params.scale);
            break;
          case 'uniform':
            pdfVal = pdfFn.call(Teaching.pdf, xVal, params.min, params.max);
            break;
          case 'exponential':
            pdfVal = pdfFn.call(Teaching.pdf, xVal, params.rate);
            break;
          case 'logNormal':
            pdfVal = pdfFn.call(Teaching.pdf, xVal, params.mu, params.sigma);
            break;
          case 'bimodal':
            pdfVal = pdfFn.call(Teaching.pdf, xVal, params.m1, params.s1, params.m2, params.s2, params.p);
            break;
          case 'poisson':
            pdfVal = pdfFn.call(Teaching.pdf, Math.round(xVal), params.lambda);
            break;
          case 'chiSquare':
            pdfVal = pdfFn.call(Teaching.pdf, xVal, params.df);
            break;
          default:
            pdfVal = 0;
        }

        const countVal = pdfVal * densityScale;
        const px = toX(xVal);
        const py = Math.max(b.y - 10, toY(countVal));

        if (!started) {
          ctx.moveTo(px, py);
          started = true;
        } else {
          ctx.lineTo(px, py);
        }
      }
      ctx.stroke();
    }

    // 3. Mark Sample Mean
    if (isFinite(sampleMean)) {
      const meanX = toX(sampleMean);
      if (meanX >= b.x && meanX <= b.x + b.width) {
        ctx.strokeStyle = '#f59e0b';
        ctx.lineWidth = 2;
        ctx.setLineDash([4, 4]);
        ctx.beginPath();
        ctx.moveTo(meanX, b.y);
        ctx.lineTo(meanX, b.y + b.height);
        ctx.stroke();
        ctx.setLineDash([]);

        ctx.fillStyle = '#f59e0b';
        ctx.font = `600 11px ${engine.options.fontFamily}`;
        ctx.textAlign = 'center';
        ctx.fillText(`M = ${sampleMean.toFixed(2)}`, meanX, b.y + 14);
      }
    }

    // Legend
    ctx.textAlign = 'right';
    ctx.font = `500 11px ${engine.options.fontFamily}`;
    ctx.fillStyle = pal.primary;
    ctx.fillText('■ Empirical Sample Histogram', b.x + b.width - 10, b.y + 15);
    ctx.fillStyle = pal.accent || '#38bdf8';
    ctx.fillText('— Theoretical PDF Overlay', b.x + b.width - 10, b.y + 32);
    ctx.fillStyle = '#f59e0b';
    ctx.fillText('┆ Sample Mean', b.x + b.width - 10, b.y + 49);
  },

  /**
   * Renders CLT Parent Population Distribution
   */
  renderCltParent(engine, parentInfo, lastSample = [], title = 'CLT Parent Population Distribution') {
    engine.lastRenderFn = () => this.renderCltParent(engine, parentInfo, lastSample, title);
    engine.clear();
    const b = engine.getPlotBounds();
    const ctx = engine.ctx;
    const pal = engine.palette;

    if (!parentInfo) return;

    const minX = parentInfo.min;
    const maxX = parentInfo.max;
    const span = maxX - minX;

    // Evaluate PDF to find peak for Y-scaling
    const steps = 150;
    let maxDensity = 0;
    const densities = [];
    for (let i = 0; i <= steps; i++) {
      const x = minX + (i / steps) * span;
      const d = parentInfo.pdf(x);
      densities.push({ x, d });
      if (d > maxDensity) maxDensity = d;
    }
    const yMax = (maxDensity || 0.5) * 1.3;

    const toX = (val) => b.x + ((val - minX) / span) * b.width;
    const toY = (d) => b.y + b.height - (d / yMax) * b.height;

    // Draw Axes
    const xTicks = [];
    for (let i = 0; i <= 4; i++) {
      const v = minX + (i / 4) * span;
      xTicks.push({ norm: i / 4, label: v.toFixed(1) });
    }
    const yTicks = [
      { norm: 0, label: '0' },
      { norm: 0.5, label: (yMax * 0.5).toFixed(2) },
      { norm: 1.0, label: yMax.toFixed(2) }
    ];

    engine.drawAxes({
      yTicks,
      xTicks,
      title: `${title} (${parentInfo.name})`,
      xLabel: 'Observation Value (X)',
      yLabel: 'Probability Density f(X)'
    });

    // Fill Distribution Area
    ctx.fillStyle = `${pal.secondary || '#94a3b8'}22`;
    ctx.beginPath();
    ctx.moveTo(b.x, b.y + b.height);
    for (const pt of densities) {
      ctx.lineTo(toX(pt.x), toY(pt.d));
    }
    ctx.lineTo(b.x + b.width, b.y + b.height);
    ctx.closePath();
    ctx.fill();

    // Stroke Distribution Curve
    ctx.strokeStyle = pal.secondary || '#94a3b8';
    ctx.lineWidth = 2.5;
    ctx.beginPath();
    for (let i = 0; i < densities.length; i++) {
      const pt = densities[i];
      if (i === 0) ctx.moveTo(toX(pt.x), toY(pt.d));
      else ctx.lineTo(toX(pt.x), toY(pt.d));
    }
    ctx.stroke();

    // Mark True Population Mean (μ)
    const muX = toX(parentInfo.mean);
    if (muX >= b.x && muX <= b.x + b.width) {
      ctx.strokeStyle = '#ef4444';
      ctx.lineWidth = 2;
      ctx.setLineDash([4, 4]);
      ctx.beginPath();
      ctx.moveTo(muX, b.y);
      ctx.lineTo(muX, b.y + b.height);
      ctx.stroke();
      ctx.setLineDash([]);

      ctx.fillStyle = '#ef4444';
      ctx.font = `600 11px ${engine.options.fontFamily}`;
      ctx.textAlign = 'center';
      ctx.fillText(`True μ = ${parentInfo.mean.toFixed(2)}`, muX, b.y + 14);
    }

    // Draw Last Sample Draw points if present
    if (lastSample && lastSample.length > 0) {
      let sampleSum = 0;
      for (const val of lastSample) {
        sampleSum += val;
        const ptX = toX(val);
        const ptD = parentInfo.pdf(val);
        const ptY = toY(ptD);

        ctx.fillStyle = '#38bdf8';
        ctx.strokeStyle = '#ffffff';
        ctx.lineWidth = 1;
        ctx.beginPath();
        ctx.arc(ptX, ptY, 4, 0, 2 * Math.PI);
        ctx.fill();
        ctx.stroke();
      }

      // Mark the mean of this single sample
      const sampleMean = sampleSum / lastSample.length;
      const sMeanX = toX(sampleMean);
      ctx.strokeStyle = '#38bdf8';
      ctx.lineWidth = 2.2;
      ctx.beginPath();
      ctx.moveTo(sMeanX, b.y + b.height - 25);
      ctx.lineTo(sMeanX, b.y + b.height);
      ctx.stroke();

      ctx.fillStyle = '#38bdf8';
      ctx.font = `600 10px ${engine.options.fontFamily}`;
      ctx.textAlign = 'center';
      ctx.fillText(`Sample x̄ = ${sampleMean.toFixed(2)}`, sMeanX, b.y + b.height - 28);
    }

    // Legend
    ctx.textAlign = 'right';
    ctx.font = `500 11px ${engine.options.fontFamily}`;
    ctx.fillStyle = pal.secondary || '#94a3b8';
    ctx.fillText('— Parent Density f(X)', b.x + b.width - 10, b.y + 15);
    ctx.fillStyle = '#ef4444';
    ctx.fillText('┆ True Population Mean (μ)', b.x + b.width - 10, b.y + 32);
    if (lastSample && lastSample.length > 0) {
      ctx.fillStyle = '#38bdf8';
      ctx.fillText(`● Current Draw (n = ${lastSample.length})`, b.x + b.width - 10, b.y + 49);
    }
  },

  /**
   * Renders CLT Sampling Distribution of Sample Means with Gaussian Overlay
   */
  renderCltSampling(engine, cltData, title = 'Sampling Distribution of the Mean (x̄)') {
    engine.lastRenderFn = () => this.renderCltSampling(engine, cltData, title);
    engine.clear();
    const b = engine.getPlotBounds();
    const ctx = engine.ctx;
    const pal = engine.palette;

    if (!cltData || cltData.samplesDrawn === 0) {
      ctx.fillStyle = pal.textDim || '#94a3b8';
      ctx.font = `14px ${engine.options.fontFamily}`;
      ctx.textAlign = 'center';
      ctx.fillText('No samples drawn yet. Click "▶ Draw 1 Sample" or "⚡ Draw 100 Samples" above.', b.x + b.width / 2, b.y + b.height / 2);
      return;
    }

    const means = cltData.values;
    const k = means.length;
    const n = cltData.sampleSize;
    const trueMu = cltData.theoreticalMean;
    const trueSE = cltData.theoreticalSE;

    let minM = Infinity;
    let maxM = -Infinity;
    for (let i = 0; i < k; i++) {
      if (means[i] < minM) minM = means[i];
      if (means[i] > maxM) maxM = means[i];
    }

    // Bounds centered around trueMu with at least 3.5 SE on each side
    const seSpan = 3.5 * trueSE;
    const plotMin = Math.min(minM - 0.2 * trueSE, trueMu - seSpan);
    const plotMax = Math.max(maxM + 0.2 * trueSE, trueMu + seSpan);
    const plotSpan = plotMax - plotMin || 1;

    // Binning
    const numBins = Math.max(15, Math.min(40, Math.round(1 + 3.322 * Math.log10(k) * 2)));
    const binWidth = plotSpan / numBins;

    const bins = new Array(numBins).fill(0);
    for (let i = 0; i < k; i++) {
      const idx = Math.min(numBins - 1, Math.max(0, Math.floor((means[i] - plotMin) / binWidth)));
      bins[idx]++;
    }

    const maxCount = Math.max(...bins, 1);
    const yMax = maxCount * 1.25;

    const toX = (val) => b.x + ((val - plotMin) / plotSpan) * b.width;
    const toY = (count) => b.y + b.height - (count / yMax) * b.height;

    // Draw Axes
    const xTicks = [];
    for (let i = 0; i <= 5; i++) {
      const v = plotMin + (i / 5) * plotSpan;
      xTicks.push({ norm: i / 5, label: v.toFixed(2) });
    }
    const yTicks = [
      { norm: 0, label: '0' },
      { norm: 0.5, label: (yMax * 0.5).toFixed(0) },
      { norm: 1.0, label: yMax.toFixed(0) }
    ];

    engine.drawAxes({
      yTicks,
      xTicks,
      title: `${title} (k = ${k.toLocaleString()} samples, n = ${n})`,
      xLabel: 'Sample Mean Value (x̄)',
      yLabel: 'Frequency of Means'
    });

    // 1. Draw Histogram Bars
    const barWidth = b.width / numBins;
    for (let i = 0; i < numBins; i++) {
      const count = bins[i];
      if (count === 0) continue;
      const h = (count / yMax) * b.height;
      const x = b.x + i * barWidth;
      const y = b.y + b.height - h;

      ctx.fillStyle = `${pal.primary}44`;
      ctx.strokeStyle = pal.primary;
      ctx.lineWidth = 1.5;
      ctx.beginPath();
      ctx.rect(x + 1, y, Math.max(1, barWidth - 2), h);
      ctx.fill();
      ctx.stroke();
    }

    // 2. Theoretical CLT Normal Curve Overlay: N(μ, σ/√n)
    if (trueSE > 0) {
      const densityScale = k * binWidth;
      ctx.strokeStyle = '#22c55e'; // Bright Emerald Green
      ctx.lineWidth = 3;
      ctx.beginPath();

      const steps = 150;
      let started = false;
      for (let s = 0; s <= steps; s++) {
        const xVal = plotMin + (s / steps) * plotSpan;
        const z = (xVal - trueMu) / trueSE;
        const normDensity = (1.0 / (trueSE * Math.sqrt(2 * Math.PI))) * Math.exp(-0.5 * z * z);
        const countVal = normDensity * densityScale;
        const px = toX(xVal);
        const py = toY(countVal);

        if (!started) {
          ctx.moveTo(px, py);
          started = true;
        } else {
          ctx.lineTo(px, py);
        }
      }
      ctx.stroke();
    }

    // 3. Mark Theoretical Mean (μ)
    const muX = toX(trueMu);
    if (muX >= b.x && muX <= b.x + b.width) {
      ctx.strokeStyle = '#ef4444';
      ctx.lineWidth = 2;
      ctx.setLineDash([4, 4]);
      ctx.beginPath();
      ctx.moveTo(muX, b.y);
      ctx.lineTo(muX, b.y + b.height);
      ctx.stroke();
      ctx.setLineDash([]);
    }

    // 4. Mark Observed Mean of Means (x̄̄)
    if (cltData.observedMean !== null) {
      const obsX = toX(cltData.observedMean);
      if (obsX >= b.x && obsX <= b.x + b.width) {
        ctx.strokeStyle = '#f59e0b';
        ctx.lineWidth = 2;
        ctx.beginPath();
        ctx.moveTo(obsX, b.y);
        ctx.lineTo(obsX, b.y + b.height);
        ctx.stroke();
      }
    }

    // Legend
    ctx.textAlign = 'right';
    ctx.font = `500 11px ${engine.options.fontFamily}`;
    ctx.fillStyle = pal.primary;
    ctx.fillText(`■ Simulated Means (k = ${k.toLocaleString()})`, b.x + b.width - 10, b.y + 15);
    ctx.fillStyle = '#22c55e';
    ctx.fillText(`— CLT Normal Fit N(μ, σ/√n)`, b.x + b.width - 10, b.y + 32);
    ctx.fillStyle = '#ef4444';
    ctx.fillText(`┆ True Mean μ = ${trueMu.toFixed(2)}`, b.x + b.width - 10, b.y + 49);
    if (cltData.observedMean !== null) {
      ctx.fillStyle = '#f59e0b';
      ctx.fillText(`— Observed x̄̄ = ${cltData.observedMean.toFixed(2)} (SE: ${cltData.observedSE.toFixed(3)})`, b.x + b.width - 10, b.y + 66);
    }
  },

  /**
   * Renders Student's t-Distribution Convergence to Standard Normal N(0, 1)
   */
  renderTConvergence(engine, metrics, options = {}) {
    engine.lastRenderFn = () => this.renderTConvergence(engine, metrics, options);
    engine.clear();
    const b = engine.getPlotBounds();
    const ctx = engine.ctx;
    const pal = engine.palette;

    if (!metrics) return;

    const df = metrics.df;
    const showTailArea = options.showTailArea !== false;
    const title = options.title || `Student's t(ν = ${df}) Convergence to Standard Normal N(0, 1)`;

    const minX = -4.5;
    const maxX = 4.5;
    const spanX = maxX - minX;
    const yMax = 0.44; // Peak of N(0, 1) is 0.39894

    const toX = (val) => b.x + ((val - minX) / spanX) * b.width;
    const toY = (d) => b.y + b.height - (Math.max(0, d) / yMax) * b.height;

    // Generate grid & axes ticks
    const xTicks = [];
    for (let x = -4; x <= 4; x += 1) {
      xTicks.push({ norm: (x - minX) / spanX, label: `${x > 0 ? '+' : ''}${x}` });
    }
    const yTicks = [
      { norm: 0, label: '0.00' },
      { norm: 0.1 / yMax, label: '0.10' },
      { norm: 0.2 / yMax, label: '0.20' },
      { norm: 0.3 / yMax, label: '0.30' },
      { norm: 0.4 / yMax, label: '0.40' }
    ];

    engine.drawAxes({
      xTicks,
      yTicks,
      title,
      xLabel: 'Standardized Value (t / z)',
      yLabel: 'Probability Density f(x)'
    });

    const steps = 240;
    const pointsNorm = [];
    const pointsT = [];

    for (let i = 0; i <= steps; i++) {
      const x = minX + (i / steps) * spanX;
      // Normal density: (1 / sqrt(2pi)) * exp(-0.5 * x^2)
      const normD = (1.0 / Math.sqrt(2 * Math.PI)) * Math.exp(-0.5 * x * x);
      // Student t density
      const tD = Teaching.pdf.studentsT(x, df, 0, 1);
      pointsNorm.push({ x, y: normD });
      pointsT.push({ x, y: tD });
    }

    // 1. Shaded Tail Area (|x| >= 1.960) under Student's t curve
    if (showTailArea) {
      const zCrit = 1.95996;
      ctx.fillStyle = 'rgba(239, 68, 68, 0.22)'; // Amber/Red glow for fat tail risk

      // Left Tail [-4.5, -1.960]
      ctx.beginPath();
      ctx.moveTo(toX(minX), toY(0));
      for (const pt of pointsT) {
        if (pt.x <= -zCrit) {
          ctx.lineTo(toX(pt.x), toY(pt.y));
        }
      }
      const tAtLeftCrit = Teaching.pdf.studentsT(-zCrit, df, 0, 1);
      ctx.lineTo(toX(-zCrit), toY(tAtLeftCrit));
      ctx.lineTo(toX(-zCrit), toY(0));
      ctx.closePath();
      ctx.fill();

      // Right Tail [1.960, 4.5]
      ctx.beginPath();
      ctx.moveTo(toX(zCrit), toY(0));
      const tAtRightCrit = Teaching.pdf.studentsT(zCrit, df, 0, 1);
      ctx.lineTo(toX(zCrit), toY(tAtRightCrit));
      for (const pt of pointsT) {
        if (pt.x >= zCrit) {
          ctx.lineTo(toX(pt.x), toY(pt.y));
        }
      }
      ctx.lineTo(toX(maxX), toY(0));
      ctx.closePath();
      ctx.fill();
    }

    // 2. Render Standard Normal Reference Curve N(0, 1) [Emerald/Cyan dashed line]
    ctx.strokeStyle = '#10b981';
    ctx.lineWidth = 2.2;
    ctx.setLineDash([5, 4]);
    ctx.beginPath();
    for (let i = 0; i < pointsNorm.length; i++) {
      const pt = pointsNorm[i];
      if (i === 0) ctx.moveTo(toX(pt.x), toY(pt.y));
      else ctx.lineTo(toX(pt.x), toY(pt.y));
    }
    ctx.stroke();
    ctx.setLineDash([]);

    // 3. Render Student's t(ν) Curve [Vibrant Violet solid line]
    ctx.strokeStyle = '#a855f7';
    ctx.lineWidth = 2.8;
    ctx.beginPath();
    for (let i = 0; i < pointsT.length; i++) {
      const pt = pointsT[i];
      if (i === 0) ctx.moveTo(toX(pt.x), toY(pt.y));
      else ctx.lineTo(toX(pt.x), toY(pt.y));
    }
    ctx.stroke();

    // 4. Mark Critical Values Lines
    const zCrit = 1.95996;
    const tCrit = metrics.tCrit;

    // Draw Gaussian +/- 1.96 lines
    [-zCrit, zCrit].forEach(zVal => {
      const xPix = toX(zVal);
      ctx.strokeStyle = '#64748b';
      ctx.lineWidth = 1.5;
      ctx.setLineDash([3, 3]);
      ctx.beginPath();
      ctx.moveTo(xPix, toY(0));
      ctx.lineTo(xPix, toY(0.18));
      ctx.stroke();
      ctx.setLineDash([]);
    });

    // Draw Student's t +/- tCrit lines if within bounds
    if (tCrit <= 4.4) {
      [-tCrit, tCrit].forEach(tVal => {
        const xPix = toX(tVal);
        ctx.strokeStyle = '#c084fc';
        ctx.lineWidth = 1.8;
        ctx.setLineDash([4, 3]);
        ctx.beginPath();
        ctx.moveTo(xPix, toY(0));
        ctx.lineTo(xPix, toY(0.24));
        ctx.stroke();
        ctx.setLineDash([]);
      });
    }

    // 5. Annotations & Peak Height Indicator
    const normPeakY = toY(metrics.normPeak);
    const tPeakY = toY(metrics.tPeak);
    const midX = toX(0);

    // Peak difference indicator line at x=0
    if (Math.abs(metrics.peakDiffPct) > 1.5) {
      ctx.strokeStyle = '#f59e0b';
      ctx.lineWidth = 1.5;
      ctx.beginPath();
      ctx.moveTo(midX, tPeakY);
      ctx.lineTo(midX, normPeakY);
      ctx.stroke();

      // Peak label
      ctx.fillStyle = '#f59e0b';
      ctx.font = `600 10px ${engine.options.fontFamily}`;
      ctx.textAlign = 'left';
      ctx.fillText(`Δ Peak: ${metrics.peakDiffPct.toFixed(1)}%`, midX + 8, (normPeakY + tPeakY) / 2 + 4);
    }

    // 6. Legend
    ctx.textAlign = 'right';
    ctx.font = `500 11px ${engine.options.fontFamily}`;

    // Normal line
    ctx.fillStyle = '#10b981';
    ctx.fillText('— — Standard Normal N(0, 1) [Peak: 0.3989]', b.x + b.width - 10, b.y + 15);

    // Student t line
    ctx.fillStyle = '#a855f7';
    ctx.fillText(`—— Student's t (ν = ${df}) [Peak: ${metrics.tPeak.toFixed(4)}]`, b.x + b.width - 10, b.y + 32);

    // Critical values
    ctx.fillStyle = '#94a3b8';
    ctx.fillText(`┆ 95% Cutoffs: z = ±1.960 vs t = ±${tCrit.toFixed(3)} (${metrics.critDiffPct >= 0 ? '+' : ''}${metrics.critDiffPct.toFixed(1)}%)`, b.x + b.width - 10, b.y + 49);

    // Tail risk
    if (showTailArea) {
      ctx.fillStyle = '#ef4444';
      ctx.fillText(`░░ Fat Tail Risk: P(|T| > 1.96) = ${(metrics.tailProb * 100).toFixed(1)}% vs 5.0%`, b.x + b.width - 10, b.y + 66);
    }
  },

  /**
   * Two-Sample Overlap, Dispersion (SD vs. SEM), and Alpha Significance Plot
   * @param {ChartEngine} engine
   * @param {object} metrics Output from Teaching.significanceOverlap.getMetrics()
   * @param {object} options
   */
  renderTwoSampleOverlap(engine, metrics, options = {}) {
    engine.lastRenderFn = () => this.renderTwoSampleOverlap(engine, metrics, options);
    engine.clear();
    const b = engine.getPlotBounds ? engine.getPlotBounds() : engine.getBounds();
    const ctx = engine.ctx;
    const pal = engine.palette;

    const {
      mean1 = 10,
      mean2 = 12,
      delta = 2,
      sd = 2.5,
      sd1 = sd || 2.5,
      sd2 = sd || 2.5,
      n = 16,
      n1 = n || 16,
      n2 = n || 16,
      sem = 0.625,
      sem1 = sem || 0.625,
      sem2 = sem || 0.625,
      alpha = 0.05,
      viewMode = 'means',
      seDiff = 0.88,
      df = 30,
      tCrit = 2.04,
      zCrit = 1.96,
      deltaCrit = 1.8,
      tStat = 0,
      pValue = 0.05,
      isSignificant = false,
      cohensD = 0.8,
      patientOVL = 0.5,
      meansOVL = 0.05,
      moe = 1.0,
      ci1 = [mean1 - 1, mean1 + 1],
      ci2 = [mean2 - 1, mean2 + 1]
    } = metrics || {};

    // View mode branching
    if (viewMode === 'null') {
      // -------------------------------------------------------------
      // NULL HYPOTHESIS VIEW: Sampling Distribution of Difference H0
      // -------------------------------------------------------------
      const safeSEDiff = Math.max(0.001, Number.isFinite(seDiff) ? seDiff : 0.88);
      const safeDelta = Number.isFinite(delta) ? delta : 2.0;
      const xSpan = Math.max(0.1, Math.max(4.2 * safeSEDiff, safeDelta + 2.5 * safeSEDiff));
      const minX = -xSpan;
      const maxX = xSpan;
      const peakY = 1.0 / (safeSEDiff * Math.sqrt(2 * Math.PI));
      const maxY = peakY * 1.25;

      const toX = (val) => {
        const v = Number.isFinite(val) ? val : minX;
        return b.x + ((v - minX) / (maxX - minX)) * b.width;
      };
      const toY = (val) => {
        const v = Number.isFinite(val) ? val : 0;
        return b.y + b.height - (v / maxY) * b.height;
      };

      // Draw Grid & Axes
      ctx.strokeStyle = pal.grid;
      ctx.lineWidth = 1;
      const xSteps = 6;
      for (let i = 0; i <= xSteps; i++) {
        const xVal = minX + (i / xSteps) * (maxX - minX);
        const xPix = toX(xVal);
        ctx.beginPath();
        ctx.moveTo(xPix, b.y);
        ctx.lineTo(xPix, b.y + b.height);
        ctx.stroke();

        ctx.fillStyle = pal.textMuted;
        ctx.font = `500 10px ${engine.options.fontFamily}`;
        ctx.textAlign = 'center';
        ctx.fillText(xVal.toFixed(2), xPix, b.y + b.height + 15);
      }

      // Generate points for H0 distribution N(0, seDiff^2)
      const numPts = 250;
      const pts = [];
      for (let i = 0; i <= numPts; i++) {
        const xVal = minX + (i / numPts) * (maxX - minX);
        const z = xVal / seDiff;
        const yVal = (1.0 / (seDiff * Math.sqrt(2 * Math.PI))) * Math.exp(-0.5 * z * z);
        pts.push({ x: xVal, y: yVal });
      }

      // Rejection Zones Shading (x <= -deltaCrit and x >= deltaCrit)
      ctx.fillStyle = 'rgba(239, 68, 68, 0.25)'; // Red Rejection Area
      // Left tail
      ctx.beginPath();
      ctx.moveTo(toX(minX), toY(0));
      for (const pt of pts) {
        if (pt.x <= -deltaCrit) ctx.lineTo(toX(pt.x), toY(pt.y));
      }
      ctx.lineTo(toX(-deltaCrit), toY(0));
      ctx.closePath();
      ctx.fill();

      // Right tail
      ctx.beginPath();
      ctx.moveTo(toX(deltaCrit), toY(0));
      for (const pt of pts) {
        if (pt.x >= deltaCrit) ctx.lineTo(toX(pt.x), toY(pt.y));
      }
      ctx.lineTo(toX(maxX), toY(0));
      ctx.closePath();
      ctx.fill();

      // Retention Zone Shading (Green tinted center)
      ctx.fillStyle = 'rgba(16, 185, 129, 0.08)';
      ctx.beginPath();
      ctx.moveTo(toX(-deltaCrit), toY(0));
      for (const pt of pts) {
        if (pt.x >= -deltaCrit && pt.x <= deltaCrit) ctx.lineTo(toX(pt.x), toY(pt.y));
      }
      ctx.lineTo(toX(deltaCrit), toY(0));
      ctx.closePath();
      ctx.fill();

      // Draw H0 Curve
      ctx.strokeStyle = '#38bdf8';
      ctx.lineWidth = 2.4;
      ctx.beginPath();
      pts.forEach((pt, idx) => {
        if (idx === 0) ctx.moveTo(toX(pt.x), toY(pt.y));
        else ctx.lineTo(toX(pt.x), toY(pt.y));
      });
      ctx.stroke();

      // Critical Boundary Lines
      [-deltaCrit, deltaCrit].forEach((cVal, idx) => {
        const xPix = toX(cVal);
        ctx.strokeStyle = '#ef4444';
        ctx.lineWidth = 1.8;
        ctx.setLineDash([4, 4]);
        ctx.beginPath();
        ctx.moveTo(xPix, toY(0));
        ctx.lineTo(xPix, toY(peakY * 0.75));
        ctx.stroke();
        ctx.setLineDash([]);

        ctx.fillStyle = '#ef4444';
        ctx.font = `600 10px ${engine.options.fontFamily}`;
        ctx.textAlign = idx === 0 ? 'right' : 'left';
        ctx.fillText(idx === 0 ? `-Δcrit (${cVal.toFixed(2)})` : `+Δcrit (${cVal.toFixed(2)})`, xPix + (idx === 0 ? -5 : 5), toY(peakY * 0.75));
      });

      // Observed Difference Indicator
      const obsXPix = toX(delta);
      const obsColor = isSignificant ? '#10b981' : '#f59e0b';
      ctx.strokeStyle = obsColor;
      ctx.lineWidth = 2.8;
      ctx.beginPath();
      ctx.moveTo(obsXPix, toY(0));
      ctx.lineTo(obsXPix, toY(peakY * 0.95));
      ctx.stroke();

      // Arrow head at top
      ctx.fillStyle = obsColor;
      ctx.beginPath();
      ctx.moveTo(obsXPix, toY(peakY * 0.98));
      ctx.lineTo(obsXPix - 6, toY(peakY * 0.90));
      ctx.lineTo(obsXPix + 6, toY(peakY * 0.90));
      ctx.closePath();
      ctx.fill();

      // Observed label
      ctx.fillStyle = obsColor;
      ctx.font = `700 11px ${engine.options.fontFamily}`;
      ctx.textAlign = 'center';
      ctx.fillText(`Observed Δ = ${delta.toFixed(2)} (t = ${tStat.toFixed(2)})`, obsXPix, toY(peakY * 1.06));

      // Status Pill at Top Right
      const badgeText = isSignificant ? `✓ SIGNIFICANT (p = ${pValue.toFixed(4)} < α)` : `✗ NOT SIGNIFICANT (p = ${pValue.toFixed(4)} ≥ α)`;
      ctx.fillStyle = isSignificant ? 'rgba(16, 185, 129, 0.2)' : 'rgba(239, 68, 68, 0.2)';
      ctx.strokeStyle = isSignificant ? '#10b981' : '#ef4444';
      ctx.lineWidth = 1;
      const badgeWidth = ctx.measureText(badgeText).width + 24;
      ctx.fillRect(b.x + b.width - badgeWidth - 10, b.y + 10, badgeWidth, 24);
      ctx.strokeRect(b.x + b.width - badgeWidth - 10, b.y + 10, badgeWidth, 24);

      ctx.fillStyle = isSignificant ? '#10b981' : '#ef4444';
      ctx.font = `700 11px ${engine.options.fontFamily}`;
      ctx.textAlign = 'right';
      ctx.fillText(badgeText, b.x + b.width - 22, b.y + 26);

      // Legend
      ctx.textAlign = 'left';
      ctx.font = `500 11px ${engine.options.fontFamily}`;
      ctx.fillStyle = '#38bdf8';
      ctx.fillText(`— Null Distribution H₀: Δ ~ N(0, SE²diff), SE = ${seDiff.toFixed(3)}`, b.x + 12, b.y + 20);
      ctx.fillStyle = '#ef4444';
      ctx.fillText(`░ Rejection Region (α = ${alpha.toFixed(3)}, t_crit = ±${tCrit.toFixed(3)})`, b.x + 12, b.y + 36);
      return;
    }

    // -------------------------------------------------------------
    // MEANS, PATIENTS, OR DUAL VIEW
    // -------------------------------------------------------------
    const activeSigma1 = Math.max(0.001, (viewMode === 'patients') ? sd1 : sem1);
    const activeSigma2 = Math.max(0.001, (viewMode === 'patients') ? sd2 : sem2);
    const maxSpread = Math.max(sd1, sd2, 4.0);
    const minX = mean1 - Math.max(3.8 * maxSpread, 4.0);
    const maxX = Math.max(mean2 + Math.max(3.8 * maxSpread, 4.0), mean1 + 7.5);
    const xSpan = Math.max(0.1, maxX - minX);

    const peak1 = 1.0 / (activeSigma1 * Math.sqrt(2 * Math.PI));
    const peak2 = 1.0 / (activeSigma2 * Math.sqrt(2 * Math.PI));
    const maxDensity = Math.max(0.001, peak1, peak2);
    const maxY = maxDensity * 1.30;

    const toX = (val) => {
      const v = Number.isFinite(val) ? val : minX;
      return b.x + ((v - minX) / xSpan) * b.width;
    };
    const toY = (val) => {
      const v = Number.isFinite(val) ? val : 0;
      return b.y + b.height - (v / maxY) * b.height;
    };

    // Draw Grid & X-axis
    ctx.strokeStyle = pal.grid;
    ctx.lineWidth = 1;
    const xSteps = 7;
    for (let i = 0; i <= xSteps; i++) {
      const xVal = minX + (i / xSteps) * xSpan;
      const xPix = toX(xVal);
      ctx.beginPath();
      ctx.moveTo(xPix, b.y);
      ctx.lineTo(xPix, b.y + b.height);
      ctx.stroke();

      ctx.fillStyle = pal.textMuted;
      ctx.font = `500 10px ${engine.options.fontFamily}`;
      ctx.textAlign = 'center';
      ctx.fillText(xVal.toFixed(1), xPix, b.y + b.height + 15);
    }

    // Generate points for both distributions
    const numPts = 320;
    const pts1 = [];
    const pts2 = [];
    const ptsOverlap = [];

    const normPDF = (x, mu, s) => (1.0 / (s * Math.sqrt(2 * Math.PI))) * Math.exp(-0.5 * Math.pow((x - mu) / s, 2));

    for (let i = 0; i <= numPts; i++) {
      const x = minX + (i / numPts) * xSpan;
      const y1 = normPDF(x, mean1, activeSigma1);
      const y2 = normPDF(x, mean2, activeSigma2);
      const yOverlap = Math.min(y1, y2);

      pts1.push({ x, y: y1 });
      pts2.push({ x, y: y2 });
      ptsOverlap.push({ x, y: yOverlap });
    }

    // 1. Shaded Overlap Area
    ctx.fillStyle = isSignificant ? 'rgba(16, 185, 129, 0.22)' : 'rgba(239, 68, 68, 0.24)';
    ctx.beginPath();
    ctx.moveTo(toX(minX), toY(0));
    ptsOverlap.forEach(pt => ctx.lineTo(toX(pt.x), toY(pt.y)));
    ctx.lineTo(toX(maxX), toY(0));
    ctx.closePath();
    ctx.fill();

    // 2. Dual Mode Background Patient Density (translucent dashed curves)
    if (viewMode === 'dual') {
      const patientPeak1 = 1.0 / (sd1 * Math.sqrt(2 * Math.PI));
      const patientPeak2 = 1.0 / (sd2 * Math.sqrt(2 * Math.PI));
      const dualScale1 = (maxDensity * 0.45) / patientPeak1;
      const dualScale2 = (maxDensity * 0.45) / patientPeak2;

      // Group 1 Patient dashed
      ctx.strokeStyle = 'rgba(6, 182, 212, 0.45)';
      ctx.lineWidth = 1.6;
      ctx.setLineDash([4, 4]);
      ctx.beginPath();
      for (let i = 0; i <= numPts; i++) {
        const x = minX + (i / numPts) * xSpan;
        const y = normPDF(x, mean1, sd1) * dualScale1;
        if (i === 0) ctx.moveTo(toX(x), toY(y));
        else ctx.lineTo(toX(x), toY(y));
      }
      ctx.stroke();

      // Group 2 Patient dashed
      ctx.strokeStyle = 'rgba(168, 85, 247, 0.45)';
      ctx.beginPath();
      for (let i = 0; i <= numPts; i++) {
        const x = minX + (i / numPts) * xSpan;
        const y = normPDF(x, mean2, sd2) * dualScale2;
        if (i === 0) ctx.moveTo(toX(x), toY(y));
        else ctx.lineTo(toX(x), toY(y));
      }
      ctx.stroke();
      ctx.setLineDash([]);
    }

    // 3. Render Distribution 1 (Group 1 / Control, Cyan)
    ctx.strokeStyle = '#06b6d4';
    ctx.lineWidth = 2.6;
    ctx.beginPath();
    pts1.forEach((pt, idx) => {
      if (idx === 0) ctx.moveTo(toX(pt.x), toY(pt.y));
      else ctx.lineTo(toX(pt.x), toY(pt.y));
    });
    ctx.stroke();

    // 4. Render Distribution 2 (Group 2 / Treatment, Violet)
    ctx.strokeStyle = '#a855f7';
    ctx.lineWidth = 2.6;
    ctx.beginPath();
    pts2.forEach((pt, idx) => {
      if (idx === 0) ctx.moveTo(toX(pt.x), toY(pt.y));
      else ctx.lineTo(toX(pt.x), toY(pt.y));
    });
    ctx.stroke();

    // 5. Mean Centers & Drop Lines
    [
      { mu: mean1, color: '#06b6d4', label: `Group 1 (μ₁ = ${mean1.toFixed(1)}${viewMode === 'patients' ? `, SD₁=${sd1.toFixed(2)}` : `, SEM₁=${sem1.toFixed(2)}`})` },
      { mu: mean2, color: '#a855f7', label: `Group 2 (μ₂ = ${mean2.toFixed(1)}${viewMode === 'patients' ? `, SD₂=${sd2.toFixed(2)}` : `, SEM₂=${sem2.toFixed(2)}`})` }
    ].forEach((grp) => {
      const xPix = toX(grp.mu);
      ctx.strokeStyle = grp.color;
      ctx.lineWidth = 1.5;
      ctx.setLineDash([3, 3]);
      ctx.beginPath();
      ctx.moveTo(xPix, toY(0));
      ctx.lineTo(xPix, toY(maxDensity));
      ctx.stroke();
      ctx.setLineDash([]);

      ctx.fillStyle = grp.color;
      ctx.font = `600 10px ${engine.options.fontFamily}`;
      ctx.textAlign = 'center';
      ctx.fillText(grp.label, xPix, toY(maxDensity) - 8);
    });

    // 6. Critical Separation Boundary Line (at mu1 + deltaCrit)
    const critX = mean1 + deltaCrit;
    if (Number.isFinite(critX) && critX <= maxX && critX >= minX) {
      const critXPix = toX(critX);
      ctx.strokeStyle = '#f59e0b';
      ctx.lineWidth = 2.0;
      ctx.setLineDash([5, 3]);
      ctx.beginPath();
      ctx.moveTo(critXPix, toY(0));
      ctx.lineTo(critXPix, toY(maxDensity * 0.85));
      ctx.stroke();
      ctx.setLineDash([]);

      ctx.fillStyle = '#f59e0b';
      ctx.font = `600 10px ${engine.options.fontFamily}`;
      ctx.textAlign = 'center';
      ctx.fillText(`Significance Boundary (Δcrit = ${deltaCrit.toFixed(2)})`, critXPix, toY(maxDensity * 0.85) - 6);
    }

    // 7. Error Bars / Confidence Intervals at the Base
    if (viewMode === 'means' || viewMode === 'dual') {
      const barY1 = toY(maxDensity * 0.05);
      const barY2 = toY(maxDensity * 0.12);

      // CI 1
      if (Array.isArray(ci1) && Number.isFinite(ci1[0]) && Number.isFinite(ci1[1])) {
        ctx.strokeStyle = '#06b6d4';
        ctx.lineWidth = 2.0;
        ctx.beginPath();
        ctx.moveTo(toX(ci1[0]), barY1);
        ctx.lineTo(toX(ci1[1]), barY1);
        ctx.stroke();
        // Caps
        [ci1[0], ci1[1]].forEach(cx => {
          ctx.beginPath();
          ctx.moveTo(toX(cx), barY1 - 4);
          ctx.lineTo(toX(cx), barY1 + 4);
          ctx.stroke();
        });
      }

      // CI 2
      if (Array.isArray(ci2) && Number.isFinite(ci2[0]) && Number.isFinite(ci2[1])) {
        ctx.strokeStyle = '#a855f7';
        ctx.lineWidth = 2.0;
        ctx.beginPath();
        ctx.moveTo(toX(ci2[0]), barY2);
        ctx.lineTo(toX(ci2[1]), barY2);
        ctx.stroke();
        // Caps
        [ci2[0], ci2[1]].forEach(cx => {
          ctx.beginPath();
          ctx.moveTo(toX(cx), barY2 - 4);
          ctx.lineTo(toX(cx), barY2 + 4);
          ctx.stroke();
        });
      }

      ctx.fillStyle = '#64748b';
      ctx.font = `500 9px ${engine.options.fontFamily}`;
      ctx.textAlign = 'left';
      ctx.fillText(`(1-α)% CIs: ±t_crit·SEM`, b.x + 10, barY2 - 8);
    }

    // 8. Overlap Badge in the Middle
    const midXPix = toX((mean1 + mean2) / 2);
    const activeOVL = viewMode === 'patients' ? patientOVL : meansOVL;
    ctx.fillStyle = isSignificant ? '#10b981' : '#ef4444';
    ctx.font = `700 11px ${engine.options.fontFamily}`;
    ctx.textAlign = 'center';
    ctx.fillText(`Overlap: ${(activeOVL * 100).toFixed(1)}%`, midXPix, toY(maxDensity * 0.35));

    // 9. Status Pill at Top Right
    const badgeText = isSignificant ? `✓ SIGNIFICANT (p = ${pValue.toFixed(4)} < α)` : `✗ NOT SIGNIFICANT (p = ${pValue.toFixed(4)} ≥ α)`;
    ctx.fillStyle = isSignificant ? 'rgba(16, 185, 129, 0.2)' : 'rgba(239, 68, 68, 0.2)';
    ctx.strokeStyle = isSignificant ? '#10b981' : '#ef4444';
    ctx.lineWidth = 1;
    const badgeWidth = ctx.measureText(badgeText).width + 24;
    ctx.fillRect(b.x + b.width - badgeWidth - 10, b.y + 10, badgeWidth, 24);
    ctx.strokeRect(b.x + b.width - badgeWidth - 10, b.y + 10, badgeWidth, 24);

    ctx.fillStyle = isSignificant ? '#10b981' : '#ef4444';
    ctx.font = `700 11px ${engine.options.fontFamily}`;
    ctx.textAlign = 'right';
    ctx.fillText(badgeText, b.x + b.width - 22, b.y + 26);

    // 10. Dynamic Legend
    ctx.textAlign = 'left';
    ctx.font = `500 11px ${engine.options.fontFamily}`;

    if (viewMode === 'patients') {
      ctx.fillStyle = '#06b6d4';
      ctx.fillText(`— Group 1 Patient Population N(μ₁, SD₁²), SD₁ = ${sd1.toFixed(2)}`, b.x + 12, b.y + 20);
      ctx.fillStyle = '#a855f7';
      ctx.fillText(`— Group 2 Patient Population N(μ₂, SD₂²), SD₂ = ${sd2.toFixed(2)}`, b.x + 12, b.y + 36);
      ctx.fillStyle = '#94a3b8';
      ctx.fillText(`░ Patient Overlap = ${(patientOVL * 100).toFixed(1)}% (Cohen's d = ${cohensD.toFixed(2)})`, b.x + 12, b.y + 52);
    } else if (viewMode === 'dual') {
      ctx.fillStyle = '#06b6d4';
      ctx.fillText(`— Group 1 Means (Solid, SEM₁ = ${sem1.toFixed(3)}) & Patients (Dashed, SD₁ = ${sd1.toFixed(2)})`, b.x + 12, b.y + 20);
      ctx.fillStyle = '#a855f7';
      ctx.fillText(`— Group 2 Means (Solid, SEM₂ = ${sem2.toFixed(3)}) & Patients (Dashed, SD₂ = ${sd2.toFixed(2)})`, b.x + 12, b.y + 36);
      ctx.fillStyle = '#94a3b8';
      ctx.fillText(`░ Means Overlap = ${(meansOVL * 100).toFixed(1)}% vs Patient Overlap = ${(patientOVL * 100).toFixed(1)}%`, b.x + 12, b.y + 52);
    } else {
      ctx.fillStyle = '#06b6d4';
      ctx.fillText(`— Group 1 Sampling Distribution (SEM₁ = ${sem1.toFixed(3)}, n₁ = ${n1})`, b.x + 12, b.y + 20);
      ctx.fillStyle = '#a855f7';
      ctx.fillText(`— Group 2 Sampling Distribution (SEM₂ = ${sem2.toFixed(3)}, n₂ = ${n2})`, b.x + 12, b.y + 36);
      ctx.fillStyle = '#f59e0b';
      ctx.fillText(`┆ Boundary: Δcrit = ${deltaCrit.toFixed(2)} at α = ${alpha.toFixed(3)} (Welch df = ${df.toFixed(1)})`, b.x + 12, b.y + 52);
    }
  },

  /**
   * Render Statistical Power Simulation (Dual Distribution, Power Curve, or 2x2 Matrix)
   * @param {ChartEngine} engine
   * @param {object} metrics Output from Teaching.powerSimulation.getMetrics()
   * @param {object} options
   */
  renderPowerSimulation(engine, metrics, options = {}) {
    engine.lastRender = () => this.renderPowerSimulation(engine, metrics, options);
    engine.lastRenderFn = engine.lastRender;
    if (engine.canvas && typeof engine.initHiDPI === 'function') {
      const rect = engine.canvas.getBoundingClientRect();
      if (rect.width > 50 && (engine.width <= 100 || Math.abs(engine.width - rect.width) > 30)) {
        engine.initHiDPI();
      }
    }
    engine.clear();
    const b = engine.getPlotBounds ? engine.getPlotBounds() : engine.getBounds();
    const ctx = engine.ctx;
    const pal = engine.palette || {};
    const font = engine.options?.fontFamily || '-apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif';
    const textMuted = pal.textMuted || pal.text || '#94a3b8';

    const {
      sd = 4.0,
      sem = 0.50,
      n = 64,
      totalN = 128,
      power = 0.80,
      beta = 0.20,
      delta = 2.0,
      alpha = 0.05,
      zCrit = 1.96,
      seDiff = 0.707,
      lambda = 2.83,
      xCrit = 1.386,
      cohensD = 0.50,
      viewMode = 'distributions',
      matrix = {},
      curvePoints = [],
      powerRating = 'ADEQUATE'
    } = metrics || {};

    // -------------------------------------------------------------
    // VIEW MODE 1: POWER VS SAMPLE SIZE CURVE
    // -------------------------------------------------------------
    if (viewMode === 'curve') {
      const minN = 4;
      const maxN = 250;
      const toX = (val) => {
        const clamped = Math.max(minN, Math.min(maxN, Number.isFinite(val) ? val : minN));
        return b.x + ((clamped - minN) / (maxN - minN)) * b.width;
      };
      const toY = (val) => {
        const clamped = Math.max(0, Math.min(1.0, Number.isFinite(val) ? val : 0));
        return b.y + b.height - clamped * b.height;
      };

      // Draw Grid & Axes
      ctx.strokeStyle = pal.grid || 'rgba(255, 255, 255, 0.08)';
      ctx.lineWidth = 1;

      // X grid
      const nTicks = [4, 25, 50, 75, 100, 150, 200, 250];
      nTicks.forEach(nVal => {
        const xPix = toX(nVal);
        ctx.beginPath();
        ctx.moveTo(xPix, b.y);
        ctx.lineTo(xPix, b.y + b.height);
        ctx.stroke();

        ctx.fillStyle = textMuted;
        ctx.font = `500 10px ${font}`;
        ctx.textAlign = 'center';
        ctx.fillText(`n=${nVal}`, xPix, b.y + b.height + 15);
      });

      // Y grid
      const pTicks = [0.2, 0.4, 0.6, 0.8, 0.9, 1.0];
      pTicks.forEach(pVal => {
        const yPix = toY(pVal);
        ctx.beginPath();
        ctx.moveTo(b.x, yPix);
        ctx.lineTo(b.x + b.width, yPix);
        ctx.stroke();

        ctx.fillStyle = textMuted;
        ctx.font = `500 10px ${font}`;
        ctx.textAlign = 'right';
        ctx.fillText(`${(pVal * 100).toFixed(0)}%`, b.x - 8, yPix + 4);
      });

      // 80% Benchmark Line (Amber/Green dashed)
      const y80 = toY(0.80);
      ctx.strokeStyle = '#10b981';
      ctx.lineWidth = 1.5;
      ctx.setLineDash([4, 4]);
      ctx.beginPath();
      ctx.moveTo(b.x, y80);
      ctx.lineTo(b.x + b.width, y80);
      ctx.stroke();
      ctx.setLineDash([]);
      ctx.fillStyle = '#10b981';
      ctx.font = `600 10px ${font}`;
      ctx.textAlign = 'right';
      ctx.fillText('80% Regulatory Standard', b.x + b.width - 10, y80 - 6);

      // 90% Benchmark Line (Cyan dashed)
      const y90 = toY(0.90);
      ctx.strokeStyle = '#06b6d4';
      ctx.lineWidth = 1.5;
      ctx.setLineDash([4, 4]);
      ctx.beginPath();
      ctx.moveTo(b.x, y90);
      ctx.lineTo(b.x + b.width, y90);
      ctx.stroke();
      ctx.setLineDash([]);
      ctx.fillStyle = '#06b6d4';
      ctx.fillText('90% High Rigor', b.x + b.width - 10, y90 - 6);

      // Shaded area under curve
      if (curvePoints && curvePoints.length > 1) {
        ctx.fillStyle = 'rgba(16, 185, 129, 0.12)';
        ctx.beginPath();
        ctx.moveTo(toX(curvePoints[0].n), toY(0));
        curvePoints.forEach(pt => ctx.lineTo(toX(pt.n), toY(pt.power)));
        ctx.lineTo(toX(curvePoints[curvePoints.length - 1].n), toY(0));
        ctx.closePath();
        ctx.fill();

        // Draw Power Curve line
        ctx.strokeStyle = '#10b981';
        ctx.lineWidth = 3.0;
        ctx.beginPath();
        curvePoints.forEach((pt, idx) => {
          if (idx === 0) ctx.moveTo(toX(pt.n), toY(pt.power));
          else ctx.lineTo(toX(pt.n), toY(pt.power));
        });
        ctx.stroke();
      }

      // Current Point Marker
      const curXPix = toX(n);
      const curYPix = toY(power);

      // Drop lines to axes
      ctx.strokeStyle = '#a855f7';
      ctx.lineWidth = 1.6;
      ctx.setLineDash([3, 3]);
      ctx.beginPath();
      ctx.moveTo(curXPix, b.y + b.height);
      ctx.lineTo(curXPix, curYPix);
      ctx.lineTo(b.x, curYPix);
      ctx.stroke();
      ctx.setLineDash([]);

      // Glowing dot
      ctx.fillStyle = 'rgba(168, 85, 247, 0.35)';
      ctx.beginPath();
      ctx.arc(curXPix, curYPix, 9, 0, Math.PI * 2);
      ctx.fill();

      ctx.fillStyle = '#a855f7';
      ctx.strokeStyle = '#ffffff';
      ctx.lineWidth = 2.0;
      ctx.beginPath();
      ctx.arc(curXPix, curYPix, 5, 0, Math.PI * 2);
      ctx.fill();
      ctx.stroke();

      // Tooltip Callout Box
      const calloutText = `n = ${n} | Power = ${(power * 100).toFixed(1)}% (β = ${(beta * 100).toFixed(1)}%)`;
      ctx.font = `700 11px ${font}`;
      const textWidth = ctx.measureText(calloutText).width;
      const boxW = textWidth + 18;
      const boxH = 26;
      const boxX = Math.min(b.x + b.width - boxW - 8, Math.max(b.x + 8, curXPix - boxW / 2));
      const boxY = Math.max(b.y + 8, curYPix - 38);

      ctx.fillStyle = '#1e293b';
      ctx.strokeStyle = '#a855f7';
      ctx.lineWidth = 1.5;
      ctx.beginPath();
      if (ctx.roundRect) {
        ctx.roundRect(boxX, boxY, boxW, boxH, [6]);
      } else {
        ctx.rect(boxX, boxY, boxW, boxH);
      }
      ctx.fill();
      ctx.stroke();

      ctx.fillStyle = '#f8fafc';
      ctx.textAlign = 'center';
      ctx.fillText(calloutText, boxX + boxW / 2, boxY + 17);

      // Legend
      ctx.textAlign = 'left';
      ctx.font = `500 11px ${font}`;
      ctx.fillStyle = '#10b981';
      ctx.fillText(`— Statistical Power Curve: P(n) at Δ = ${delta.toFixed(2)}, SD = ${sd.toFixed(2)}, α = ${alpha.toFixed(3)}`, b.x + 12, b.y + 20);
      ctx.fillStyle = '#a855f7';
      ctx.fillText(`● Current Operating Point (n = ${n} per group, N = ${totalN}, SEM = ${sem.toFixed(3)})`, b.x + 12, b.y + 36);
      return;
    }

    // -------------------------------------------------------------
    // VIEW MODE 2: 2x2 DECISION ERROR MATRIX
    // -------------------------------------------------------------
    if (viewMode === 'matrix') {
      const pad = 12;
      const cardW = (b.width - pad * 3) / 2;
      const cardH = (b.height - pad * 3) / 2;

      const cells = [
        {
          col: 0, row: 0,
          title: 'SPECIFICITY (1 − α)',
          val: `${((1 - alpha) * 100).toFixed(1)}%`,
          sub: 'True Negative Rate',
          desc: 'H₀ is TRUE (no effect), decision is RETAIN H₀. Correct clinical conclusion: drug is correctly recognized as having no effect.',
          color: '#38bdf8',
          bg: 'rgba(56, 189, 248, 0.08)',
          border: 'rgba(56, 189, 248, 0.35)',
          pct: 1 - alpha
        },
        {
          col: 1, row: 0,
          title: 'TYPE I ERROR (α)',
          val: `${(alpha * 100).toFixed(1)}%`,
          sub: 'False Positive Rate (Significance Level)',
          desc: 'H₀ is TRUE (no effect), but decision is REJECT H₀. False alarm: ineffective drug erroneously declared effective.',
          color: '#ef4444',
          bg: 'rgba(239, 68, 68, 0.08)',
          border: 'rgba(239, 68, 68, 0.35)',
          pct: alpha
        },
        {
          col: 0, row: 1,
          title: 'TYPE II ERROR (β)',
          val: `${(beta * 100).toFixed(1)}%`,
          sub: 'False Negative Rate',
          desc: 'H₁ is TRUE (real effect Δ), but decision is RETAIN H₀. Missed discovery: effective therapy discarded due to lack of power!',
          color: '#f59e0b',
          bg: 'rgba(245, 158, 11, 0.08)',
          border: 'rgba(245, 158, 11, 0.35)',
          pct: beta
        },
        {
          col: 1, row: 1,
          title: 'STATISTICAL POWER (1 − β)',
          val: `${(power * 100).toFixed(1)}%`,
          sub: 'True Positive Rate (Sensitivity)',
          desc: 'H₁ is TRUE (real effect Δ), and decision is REJECT H₀. Successful trial: effective therapy correctly discovered and verified!',
          color: '#10b981',
          bg: 'rgba(16, 185, 129, 0.08)',
          border: 'rgba(16, 185, 129, 0.35)',
          pct: power
        }
      ];

      cells.forEach(c => {
        const cx = b.x + pad + c.col * (cardW + pad);
        const cy = b.y + pad + c.row * (cardH + pad);

        ctx.fillStyle = c.bg;
        ctx.strokeStyle = c.border;
        ctx.lineWidth = 1.5;
        ctx.beginPath();
        if (ctx.roundRect) {
          ctx.roundRect(cx, cy, cardW, cardH, [8]);
        } else {
          ctx.rect(cx, cy, cardW, cardH);
        }
        ctx.fill();
        ctx.stroke();

        // Card Title
        ctx.fillStyle = c.color;
        ctx.font = `700 11px ${font}`;
        ctx.textAlign = 'left';
        ctx.fillText(c.title, cx + 12, cy + 20);

        // Subtitle
        ctx.fillStyle = textMuted;
        ctx.font = `500 9px ${font}`;
        ctx.fillText(c.sub, cx + 12, cy + 34);

        // Large Percentage Value
        ctx.fillStyle = c.color;
        ctx.font = `800 24px ${font}`;
        ctx.fillText(c.val, cx + 12, cy + 64);

        // Progress bar indicator
        const barX = cx + 12;
        const barY = cy + 72;
        const barW = cardW - 24;
        const barH = 6;
        ctx.fillStyle = 'rgba(255, 255, 255, 0.1)';
        ctx.fillRect(barX, barY, barW, barH);
        ctx.fillStyle = c.color;
        ctx.fillRect(barX, barY, barW * Math.max(0, Math.min(1.0, c.pct)), barH);

        // Description text wrapped
        ctx.fillStyle = '#cbd5e1';
        ctx.font = `400 9.5px ${font}`;
        const words = c.desc.split(' ');
        let line = '';
        let lineY = cy + 93;
        for (let w = 0; w < words.length; w++) {
          const testLine = line + words[w] + ' ';
          if (ctx.measureText(testLine).width > cardW - 24 && w > 0) {
            ctx.fillText(line, cx + 12, lineY);
            line = words[w] + ' ';
            lineY += 12;
          } else {
            line = testLine;
          }
        }
        ctx.fillText(line, cx + 12, lineY);
      });
      return;
    }

    // -------------------------------------------------------------
    // VIEW MODE 3: DUAL DISTRIBUTION (H0 VS H1) WITH SHADED REGIONS
    // -------------------------------------------------------------
    const safeSEDiff = Math.max(0.001, Number.isFinite(seDiff) ? seDiff : 0.707);
    const safeDelta = Number.isFinite(delta) ? delta : 2.0;

    const minX = -3.5 * safeSEDiff;
    const maxX = safeDelta + 3.8 * safeSEDiff;
    const xSpan = Math.max(0.1, maxX - minX);

    const peakY = 1.0 / (safeSEDiff * Math.sqrt(2 * Math.PI));
    const maxY = peakY * 1.32;

    const toX = (val) => {
      const v = Number.isFinite(val) ? val : minX;
      return b.x + ((v - minX) / xSpan) * b.width;
    };
    const toY = (val) => {
      const v = Number.isFinite(val) ? val : 0;
      return b.y + b.height - (v / maxY) * b.height;
    };

    // Grid lines & X-axis
    ctx.strokeStyle = pal.grid || 'rgba(255, 255, 255, 0.08)';
    ctx.lineWidth = 1;
    const xSteps = 7;
    for (let i = 0; i <= xSteps; i++) {
      const xVal = minX + (i / xSteps) * xSpan;
      const xPix = toX(xVal);
      ctx.beginPath();
      ctx.moveTo(xPix, b.y);
      ctx.lineTo(xPix, b.y + b.height);
      ctx.stroke();

      ctx.fillStyle = textMuted;
      ctx.font = `500 10px ${font}`;
      ctx.textAlign = 'center';
      ctx.fillText(xVal.toFixed(2), xPix, b.y + b.height + 15);
    }

    // Generate Points for H0: N(0, seDiff^2) and H1: N(delta, seDiff^2)
    const numPts = 320;
    const ptsH0 = [];
    const ptsH1 = [];
    const normPDF = (x, mu, s) => (1.0 / (s * Math.sqrt(2 * Math.PI))) * Math.exp(-0.5 * Math.pow((x - mu) / s, 2));

    for (let i = 0; i <= numPts; i++) {
      const x = minX + (i / numPts) * xSpan;
      ptsH0.push({ x, y: normPDF(x, 0, safeSEDiff) });
      ptsH1.push({ x, y: normPDF(x, safeDelta, safeSEDiff) });
    }

    const safeXCrit = Number.isFinite(xCrit) ? xCrit : (zCrit * safeSEDiff);

    // 1. Shading: Statistical Power (1 - beta) under H1 (where x >= safeXCrit, Emerald Green)
    ctx.fillStyle = 'rgba(16, 185, 129, 0.38)';
    ctx.beginPath();
    ctx.moveTo(toX(safeXCrit), toY(0));
    ctx.lineTo(toX(safeXCrit), toY(normPDF(safeXCrit, safeDelta, safeSEDiff)));
    for (const pt of ptsH1) {
      if (pt.x >= safeXCrit) ctx.lineTo(toX(pt.x), toY(pt.y));
    }
    ctx.lineTo(toX(maxX), toY(0));
    ctx.closePath();
    ctx.fill();

    // 2. Shading: Type II Error (beta) under H1 (where x < safeXCrit, Amber)
    ctx.fillStyle = 'rgba(245, 158, 11, 0.35)';
    ctx.beginPath();
    ctx.moveTo(toX(minX), toY(0));
    for (const pt of ptsH1) {
      if (pt.x <= safeXCrit) ctx.lineTo(toX(pt.x), toY(pt.y));
    }
    ctx.lineTo(toX(safeXCrit), toY(normPDF(safeXCrit, safeDelta, safeSEDiff)));
    ctx.lineTo(toX(safeXCrit), toY(0));
    ctx.closePath();
    ctx.fill();

    // 3. Shading: Type I Error (alpha/2) under H0 (where x >= safeXCrit, Red rejection tail)
    ctx.fillStyle = 'rgba(239, 68, 68, 0.35)';
    ctx.beginPath();
    ctx.moveTo(toX(safeXCrit), toY(0));
    ctx.lineTo(toX(safeXCrit), toY(normPDF(safeXCrit, 0, safeSEDiff)));
    for (const pt of ptsH0) {
      if (pt.x >= safeXCrit) ctx.lineTo(toX(pt.x), toY(pt.y));
    }
    ctx.lineTo(toX(maxX), toY(0));
    ctx.closePath();
    ctx.fill();

    // 4. Shading: Retention Zone (1 - alpha/2) under H0 (Soft blue tint)
    ctx.fillStyle = 'rgba(56, 189, 248, 0.08)';
    ctx.beginPath();
    ctx.moveTo(toX(minX), toY(0));
    for (const pt of ptsH0) {
      if (pt.x <= safeXCrit) ctx.lineTo(toX(pt.x), toY(pt.y));
    }
    ctx.lineTo(toX(safeXCrit), toY(normPDF(safeXCrit, 0, safeSEDiff)));
    ctx.lineTo(toX(safeXCrit), toY(0));
    ctx.closePath();
    ctx.fill();

    // 5. Draw H0 Curve (Cyan)
    ctx.strokeStyle = '#06b6d4';
    ctx.lineWidth = 2.6;
    ctx.beginPath();
    ptsH0.forEach((pt, idx) => {
      if (idx === 0) ctx.moveTo(toX(pt.x), toY(pt.y));
      else ctx.lineTo(toX(pt.x), toY(pt.y));
    });
    ctx.stroke();

    // 6. Draw H1 Curve (Violet/Purple)
    ctx.strokeStyle = '#a855f7';
    ctx.lineWidth = 2.6;
    ctx.beginPath();
    ptsH1.forEach((pt, idx) => {
      if (idx === 0) ctx.moveTo(toX(pt.x), toY(pt.y));
      else ctx.lineTo(toX(pt.x), toY(pt.y));
    });
    ctx.stroke();

    // 7. Center Means Vertical Drop Lines (mu=0 and mu=delta)
    [
      { mu: 0, color: '#06b6d4', label: 'Null H₀: Δ = 0' },
      { mu: safeDelta, color: '#a855f7', label: `Alternative H₁: Δ = ${safeDelta.toFixed(2)}` }
    ].forEach(grp => {
      const xPix = toX(grp.mu);
      ctx.strokeStyle = grp.color;
      ctx.lineWidth = 1.6;
      ctx.setLineDash([3, 3]);
      ctx.beginPath();
      ctx.moveTo(xPix, toY(0));
      ctx.lineTo(xPix, toY(peakY));
      ctx.stroke();
      ctx.setLineDash([]);

      ctx.fillStyle = grp.color;
      ctx.font = `600 10px ${font}`;
      ctx.textAlign = 'center';
      ctx.fillText(grp.label, xPix, toY(peakY) - 8);
    });

    // 8. Effect Size Bracket Δ between mu0 and mu1
    const yBracket = toY(peakY * 0.45);
    ctx.strokeStyle = '#a855f7';
    ctx.lineWidth = 1.6;
    ctx.beginPath();
    ctx.moveTo(toX(0), yBracket);
    ctx.lineTo(toX(safeDelta), yBracket);
    ctx.stroke();
    // Bracket tick ends
    [0, safeDelta].forEach(muVal => {
      ctx.beginPath();
      ctx.moveTo(toX(muVal), yBracket - 4);
      ctx.lineTo(toX(muVal), yBracket + 4);
      ctx.stroke();
    });
    ctx.fillStyle = '#a855f7';
    ctx.font = `700 10px ${font}`;
    ctx.textAlign = 'center';
    ctx.fillText(`True Effect Δ = ${safeDelta.toFixed(2)} (d = ${cohensD.toFixed(2)})`, toX(safeDelta / 2), yBracket - 6);

    // 9. Critical Threshold Line xcrit (Red dashed)
    if (Number.isFinite(safeXCrit) && safeXCrit >= minX && safeXCrit <= maxX) {
      const critXPix = toX(safeXCrit);
      ctx.strokeStyle = '#ef4444';
      ctx.lineWidth = 2.2;
      ctx.setLineDash([5, 3]);
      ctx.beginPath();
      ctx.moveTo(critXPix, toY(0));
      ctx.lineTo(critXPix, toY(peakY * 1.05));
      ctx.stroke();
      ctx.setLineDash([]);

      ctx.fillStyle = '#ef4444';
      ctx.font = `700 10px ${font}`;
      ctx.textAlign = 'center';
      ctx.fillText(`Significance Cutoff xcrit = ${safeXCrit.toFixed(2)}`, critXPix, toY(peakY * 1.05) - 6);
    }

    // 10. Shading Zone Labels
    // Power label in green area
    const pwrX = toX(Math.max(safeXCrit + 0.3 * safeSEDiff, safeDelta));
    if (pwrX < b.x + b.width - 60) {
      ctx.fillStyle = '#10b981';
      ctx.font = `700 11px ${font}`;
      ctx.textAlign = 'center';
      ctx.fillText(`Power (1 − β): ${(power * 100).toFixed(1)}%`, pwrX, toY(peakY * 0.28));
    }

    // Beta label in amber area
    const betaX = toX(Math.min(safeXCrit - 0.2 * safeSEDiff, safeDelta - 0.3 * safeSEDiff));
    if (betaX > b.x + 50) {
      ctx.fillStyle = '#f59e0b';
      ctx.font = `700 11px ${font}`;
      ctx.textAlign = 'center';
      ctx.fillText(`Beta (β): ${(beta * 100).toFixed(1)}%`, betaX, toY(peakY * 0.16));
    }

    // 11. Status Badge at Top-Right
    const isAdequate = power >= 0.80;
    const badgeText = isAdequate
      ? `✓ ADEQUATE POWER: ${(power * 100).toFixed(1)}% (β = ${(beta * 100).toFixed(1)}%)`
      : `⚠ UNDERPOWERED: ${(power * 100).toFixed(1)}% (β = ${(beta * 100).toFixed(1)}%)`;
    ctx.fillStyle = isAdequate ? 'rgba(16, 185, 129, 0.2)' : 'rgba(239, 68, 68, 0.2)';
    ctx.strokeStyle = isAdequate ? '#10b981' : '#ef4444';
    ctx.lineWidth = 1;
    ctx.font = `700 11px ${font}`;
    const badgeWidth = ctx.measureText(badgeText).width + 24;
    ctx.fillRect(b.x + b.width - badgeWidth - 10, b.y + 10, badgeWidth, 24);
    ctx.strokeRect(b.x + b.width - badgeWidth - 10, b.y + 10, badgeWidth, 24);

    ctx.fillStyle = isAdequate ? '#10b981' : '#ef4444';
    ctx.textAlign = 'right';
    ctx.fillText(badgeText, b.x + b.width - 22, b.y + 26);

    // 12. Dynamic Legend
    ctx.textAlign = 'left';
    ctx.font = `500 11px ${font}`;
    ctx.fillStyle = '#06b6d4';
    ctx.fillText(`— Null Distribution H₀: Δ ~ N(0, SE²diff), SE = ${safeSEDiff.toFixed(3)} (SEM = ${sem.toFixed(3)})`, b.x + 12, b.y + 20);
    ctx.fillStyle = '#a855f7';
    ctx.fillText(`— Alternative Distribution H₁: Δ ~ N(${safeDelta.toFixed(2)}, SE²diff), n = ${n} per group (N = ${totalN})`, b.x + 12, b.y + 36);
    ctx.fillStyle = '#10b981';
    ctx.fillText(`░ Green Shaded Area: Power (1 − β) = ${(power * 100).toFixed(1)}% | ░ Amber Area: β = ${(beta * 100).toFixed(1)}%`, b.x + 12, b.y + 52);
  },

  /**
   * Section 6: Bayesian Statistics & Logic (3Blue1Brown Model)
   * Visualizing Bayes' Theorem via the 1x1 Unit Square, Natural Frequencies, and Odds Updating.
   *
   * @param {ChartEngine} engine
   * @param {object} metrics Output from Teaching.bayesianSimulation.getMetrics()
   * @param {object} options
   */
  renderBayesianSimulation(engine, metrics, options = {}) {
    engine.lastRender = () => this.renderBayesianSimulation(engine, metrics, options);
    engine.lastRenderFn = engine.lastRender;
    if (engine.canvas && typeof engine.initHiDPI === 'function') {
      const rect = engine.canvas.getBoundingClientRect();
      if (rect.width > 50 && (engine.width <= 100 || Math.abs(engine.width - rect.width) > 30)) {
        engine.initHiDPI();
      }
    }
    engine.clear();
    const b = engine.getPlotBounds ? engine.getPlotBounds() : engine.getBounds();
    const ctx = engine.ctx;
    const pal = engine.palette || {};
    const font = engine.options?.fontFamily || '-apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif';
    const textMuted = pal.textMuted || pal.text || '#94a3b8';

    const {
      prior = 0.0476,
      notPrior = 0.9524,
      likelihood = 0.40,
      falsePositive = 0.10,
      sampleSize = 210,
      areaHAndE = 0.01905,
      areaNotHAndE = 0.09524,
      pEvidence = 0.11429,
      posterior = 0.1667,
      posteriorNotH = 0.8333,
      priorOdds = 0.05,
      bayesFactor = 4.0,
      posteriorOdds = 0.20,
      beliefShift = 0.119,
      countH = 10,
      countNotH = 200,
      countHAndE = 4,
      countNotHAndE = 20,
      countTotalE = 24,
      trajectory = [],
      evidenceRating = 'Substantial / Moderate',
      viewMode = 'square',
      preset = 'steve'
    } = metrics || {};

    // -------------------------------------------------------------
    // VIEW MODE 1: 3BLUE1BROWN 1x1 UNIT SQUARE (GEOMETRY OF BAYES)
    // -------------------------------------------------------------
    if (viewMode === 'square') {
      const padTop = 32;
      const padBottom = 40;
      const labelLeftW = 90;
      const labelRightW = 88;

      const availableH = b.height - padTop - padBottom;
      const showRightCard = b.width >= 720;
      const rightCardW = showRightCard ? Math.min(380, Math.max(260, Math.floor(b.width * 0.36))) : 0;
      const gap = showRightCard ? 28 : 0;

      const availableW = b.width - labelLeftW - labelRightW - rightCardW - gap;
      const sqSize = Math.floor(Math.max(220, Math.min(availableH, availableW)));

      const sqX = b.x + labelLeftW;
      const sqY = b.y + padTop;

      // Draw Main 1x1 Possibility Space Square
      ctx.fillStyle = 'rgba(15, 23, 42, 0.75)';
      ctx.strokeStyle = 'rgba(255, 255, 255, 0.30)';
      ctx.lineWidth = 2.0;
      ctx.beginPath();
      if (ctx.roundRect) ctx.roundRect(sqX, sqY, sqSize, sqSize, [6]);
      else ctx.rect(sqX, sqY, sqSize, sqSize);
      ctx.fill();
      ctx.stroke();

      // Horizontal division coordinate for Prior P(H) vs P(~H)
      const splitX = sqX + sqSize * prior;
      const colWidthH = splitX - sqX;
      const colWidthNotH = sqX + sqSize - splitX;

      // Vertical heights for Evidence Shading
      const hShadeH = sqSize * likelihood;
      const hShadeNotH = sqSize * falsePositive;

      const yTopShadeH = sqY + sqSize - hShadeH;
      const yTopShadeNotH = sqY + sqSize - hShadeNotH;

      // 1. Shaded Region: P(H and E) = P(H) * P(E|H) (Emerald Green)
      if (colWidthH > 0.5 && hShadeH > 0.5) {
        let fillGreen = 'rgba(16, 185, 129, 0.55)';
        try {
          if (typeof ctx.createLinearGradient === 'function') {
            const gradGreen = ctx.createLinearGradient(sqX, yTopShadeH, sqX, sqY + sqSize);
            if (gradGreen && typeof gradGreen.addColorStop === 'function') {
              gradGreen.addColorStop(0, 'rgba(16, 185, 129, 0.70)');
              gradGreen.addColorStop(1, 'rgba(16, 185, 129, 0.40)');
              fillGreen = gradGreen;
            }
          }
        } catch { /* fallback */ }
        ctx.fillStyle = fillGreen;
        ctx.fillRect(sqX, yTopShadeH, colWidthH, hShadeH);

        ctx.strokeStyle = '#10b981';
        ctx.lineWidth = 2.0;
        ctx.strokeRect(sqX, yTopShadeH, colWidthH, hShadeH);
      }

      // 2. Shaded Region: P(~H and E) = P(~H) * P(E|~H) (Amber/Orange)
      if (colWidthNotH > 0.5 && hShadeNotH > 0.5) {
        let fillAmber = 'rgba(245, 158, 11, 0.55)';
        try {
          if (typeof ctx.createLinearGradient === 'function') {
            const gradAmber = ctx.createLinearGradient(splitX, yTopShadeNotH, splitX, sqY + sqSize);
            if (gradAmber && typeof gradAmber.addColorStop === 'function') {
              gradAmber.addColorStop(0, 'rgba(245, 158, 11, 0.70)');
              gradAmber.addColorStop(1, 'rgba(245, 158, 11, 0.40)');
              fillAmber = gradAmber;
            }
          }
        } catch { /* fallback */ }
        ctx.fillStyle = fillAmber;
        ctx.fillRect(splitX, yTopShadeNotH, colWidthNotH, hShadeNotH);

        ctx.strokeStyle = '#f59e0b';
        ctx.lineWidth = 2.0;
        ctx.strokeRect(splitX, yTopShadeNotH, colWidthNotH, hShadeNotH);
      }

      // 3. Vertical Dividing Line between H and ~H
      ctx.strokeStyle = 'rgba(255, 255, 255, 0.5)';
      ctx.lineWidth = 2.0;
      ctx.setLineDash([5, 4]);
      ctx.beginPath();
      ctx.moveTo(splitX, sqY);
      ctx.lineTo(splitX, sqY + sqSize);
      ctx.stroke();
      ctx.setLineDash([]);

      // Top Dimension Labels: Prior Hypotheses
      ctx.font = `700 12.5px ${font}`;
      ctx.textAlign = 'center';

      const rawHX = sqX + colWidthH / 2;
      const rawNotHX = splitX + colWidthNotH / 2;
      const minHeaderDist = 80;
      let headerHX = rawHX;
      let headerNotHX = rawNotHX;
      if (headerNotHX - headerHX < minHeaderDist) {
        headerHX = Math.max(sqX - 6, headerNotHX - minHeaderDist);
      }

      // Left Column Header: P(H)
      ctx.fillStyle = '#10b981';
      ctx.fillText(`P(H) = ${(prior * 100).toFixed(1)}%`, headerHX, sqY - 16);
      ctx.font = `500 10.5px ${font}`;
      ctx.fillStyle = textMuted;
      ctx.fillText(preset === 'steve' ? 'Librarians' : 'Hypothesis H', headerHX, sqY - 3);

      // Right Column Header: P(~H)
      ctx.font = `700 12.5px ${font}`;
      ctx.fillStyle = '#f59e0b';
      ctx.fillText(`P(¬H) = ${(notPrior * 100).toFixed(1)}%`, headerNotHX, sqY - 16);
      ctx.font = `500 10.5px ${font}`;
      ctx.fillStyle = textMuted;
      ctx.fillText(preset === 'steve' ? 'Farmers (20× more)' : 'Alternative ¬H', headerNotHX, sqY - 3);

      // Height Labels for Likelihoods
      // P(E|H) on left edge
      ctx.textAlign = 'right';
      ctx.fillStyle = '#10b981';
      ctx.font = `700 12px ${font}`;
      const yMidH = yTopShadeH + hShadeH / 2 + 4;
      ctx.fillText(`P(E|H) = ${(likelihood * 100).toFixed(0)}%`, sqX - 10, yMidH);
      ctx.font = `500 10px ${font}`;
      ctx.fillStyle = textMuted;
      ctx.fillText(preset === 'steve' ? 'True Pos' : 'Likelihood', sqX - 10, yMidH + 13);

      // P(E|~H) on right edge
      ctx.textAlign = 'left';
      ctx.fillStyle = '#f59e0b';
      ctx.font = `700 12px ${font}`;
      const yMidNotH = yTopShadeNotH + hShadeNotH / 2 + 4;
      ctx.fillText(`P(E|¬H) = ${(falsePositive * 100).toFixed(0)}%`, sqX + sqSize + 10, yMidNotH);
      ctx.font = `500 10px ${font}`;
      ctx.fillStyle = textMuted;
      ctx.fillText(preset === 'steve' ? 'False Alarm' : 'False Pos', sqX + sqSize + 10, yMidNotH + 13);

      // Unshaded Region Labels: Ruled Out by Evidence
      ctx.fillStyle = 'rgba(148, 163, 184, 0.45)';
      ctx.font = `italic 11px ${font}`;
      ctx.textAlign = 'center';
      if (sqSize - hShadeNotH > 40) {
        ctx.fillText('Ruled out by Evidence ¬E (Dimmed space)', splitX + colWidthNotH / 2, sqY + (sqSize - hShadeNotH) / 2);
      }

      // Bottom Area Value Labels
      let areaHX = rawHX;
      let areaNotHX = rawNotHX;
      if (areaNotHX - areaHX < 85) {
        areaHX = Math.max(sqX - 6, areaNotHX - 85);
      }

      ctx.font = `700 12px ${font}`;
      ctx.fillStyle = '#10b981';
      ctx.textAlign = 'center';
      ctx.fillText(`Area: ${(areaHAndE * 100).toFixed(1)}%`, areaHX, sqY + sqSize + 17);
      ctx.font = `500 10px ${font}`;
      ctx.fillStyle = textMuted;
      ctx.fillText('P(H ∩ E)', areaHX, sqY + sqSize + 30);

      ctx.font = `700 12px ${font}`;
      ctx.fillStyle = '#f59e0b';
      ctx.textAlign = 'center';
      ctx.fillText(`Area: ${(areaNotHAndE * 100).toFixed(1)}%`, areaNotHX, sqY + sqSize + 17);
      ctx.font = `500 10px ${font}`;
      ctx.fillStyle = textMuted;
      ctx.fillText('P(¬H ∩ E)', areaNotHX, sqY + sqSize + 30);

      // -------------------------------------------------------------
      // RIGHT SIDE: RESTRICTED SPACE PROPORTION CALLOUT (BAYES RULE)
      // -------------------------------------------------------------
      if (showRightCard && rightCardW >= 220) {
        const rightX = sqX + sqSize + labelRightW + gap;
        const rightY = sqY - 8;
        const cardH = sqSize + 44;

        // Card Box
        ctx.fillStyle = 'rgba(30, 41, 59, 0.75)';
        ctx.strokeStyle = 'rgba(56, 189, 248, 0.35)';
        ctx.lineWidth = 1.5;
        ctx.beginPath();
        if (ctx.roundRect) ctx.roundRect(rightX, rightY, rightCardW, cardH, [10]);
        else ctx.rect(rightX, rightY, rightCardW, cardH);
        ctx.fill();
        ctx.stroke();

        // Card Header
        ctx.textAlign = 'left';
        ctx.font = `700 12.5px ${font}`;
        ctx.fillStyle = '#38bdf8';
        ctx.fillText('RESTRICTED POSSIBILITY SPACE P(E)', rightX + 16, rightY + 24);

        ctx.font = `400 10px ${font}`;
        ctx.fillStyle = textMuted;
        ctx.fillText('The evidence discards all unshaded space.', rightX + 16, rightY + 40);

        // Total Evidence Area P(E)
        ctx.font = `700 12px ${font}`;
        ctx.fillStyle = '#cbd5e1';
        ctx.fillText(`Total Evidence Area = ${(pEvidence * 100).toFixed(2)}%`, rightX + 16, rightY + 66);

        // Visual Proportion Bar
        const barX = rightX + 16;
        const barY = rightY + 78;
        const barW = rightCardW - 32;
        const barH = 26;

        const pwrW = barW * Math.max(0, Math.min(1.0, posterior));
        const remW = barW - pwrW;

        // Green segment: Posterior P(H|E)
        ctx.fillStyle = '#10b981';
        ctx.fillRect(barX, barY, pwrW, barH);
        // Amber segment: P(~H|E)
        ctx.fillStyle = '#f59e0b';
        ctx.fillRect(barX + pwrW, barY, remW, barH);

        // Bar border
        ctx.strokeStyle = '#ffffff';
        ctx.lineWidth = 1.0;
        ctx.strokeRect(barX, barY, barW, barH);

        // Labels inside bar
        ctx.font = `700 11px ${font}`;
        if (pwrW > 35) {
          ctx.fillStyle = '#ffffff';
          ctx.textAlign = 'center';
          ctx.fillText(`${(posterior * 100).toFixed(1)}%`, barX + pwrW / 2, barY + 17);
        }
        if (remW > 45) {
          ctx.fillStyle = '#0f172a';
          ctx.textAlign = 'center';
          ctx.fillText(`${(posteriorNotH * 100).toFixed(1)}%`, barX + pwrW + remW / 2, barY + 17);
        }

        // Legend under bar
        ctx.textAlign = 'left';
        ctx.font = `600 10.5px ${font}`;
        ctx.fillStyle = '#10b981';
        ctx.fillText(`■ P(H|E): ${(posterior * 100).toFixed(1)}%`, barX, barY + 44);
        ctx.fillStyle = '#f59e0b';
        ctx.fillText(`■ P(¬H|E): ${(posteriorNotH * 100).toFixed(1)}%`, barX + barW / 2, barY + 44);

        // 3Blue1Brown Equation Breakdown
        const eqY = barY + 54;
        const eqH = 96;
        ctx.fillStyle = 'rgba(255, 255, 255, 0.06)';
        ctx.fillRect(barX, eqY, barW, eqH);
        ctx.strokeStyle = 'rgba(255, 255, 255, 0.15)';
        ctx.strokeRect(barX, eqY, barW, eqH);

        ctx.fillStyle = '#94a3b8';
        ctx.font = `600 10px ${font}`;
        ctx.fillText('BAYES\' PROPORTION RULE (3BLUE1BROWN):', barX + 10, eqY + 16);

        ctx.font = `700 11.5px monospace`;
        ctx.fillStyle = '#00d2ff';
        ctx.fillText('P(H|E) = Green Area / Total Shaded Area', barX + 10, eqY + 33);

        ctx.fillStyle = '#f8fafc';
        ctx.font = `600 10.5px monospace`;
        ctx.fillText(`• By Area %:   ${(areaHAndE * 100).toFixed(2)}% / ${(pEvidence * 100).toFixed(2)}% = ${(posterior * 100).toFixed(1)}%`, barX + 10, eqY + 52);

        // Calculation by Number / Natural Counts
        const rawCountHAndE = sampleSize * areaHAndE;
        const rawCountNotHAndE = sampleSize * areaNotHAndE;
        const rawCountTotalE = rawCountHAndE + rawCountNotHAndE;
        const fmtCount = (v) => Math.abs(v - Math.round(v)) < 0.05 ? Math.round(v).toString() : v.toFixed(1);
        const strH = fmtCount(rawCountHAndE);
        const strNotH = fmtCount(rawCountNotHAndE);
        const strTot = fmtCount(rawCountTotalE);

        ctx.fillStyle = '#38bdf8';
        ctx.font = `600 10.5px monospace`;
        ctx.fillText(`• By Number:   ${strH} / (${strH} + ${strNotH}) = ${strH} / ${strTot} = ${(posterior * 100).toFixed(1)}%`, barX + 10, eqY + 70);

        ctx.fillStyle = textMuted;
        ctx.font = `400 9px ${font}`;
        ctx.fillText(`(Representative cohort of N = ${sampleSize} individuals)`, barX + 10, eqY + 86);

        // Bottom takeaway callout
        const shiftY = eqY + eqH + 18;
        ctx.fillStyle = beliefShift >= 0 ? '#10b981' : '#ef4444';
        ctx.font = `700 12px ${font}`;
        ctx.fillText(
          `Belief Shift: ${(beliefShift >= 0 ? '+' : '')}${(beliefShift * 100).toFixed(1)}%`,
          barX,
          shiftY
        );

        ctx.font = `600 11px ${font}`;
        ctx.fillStyle = '#cbd5e1';
        ctx.fillText(
          `Prior Odds ${priorOdds < 0.1 ? '1:' + (1/priorOdds).toFixed(1) : priorOdds.toFixed(2)} → Posterior Odds ${posteriorOdds < 0.1 ? '1:' + (1/posteriorOdds).toFixed(1) : posteriorOdds.toFixed(2)}`,
          barX,
          shiftY + 18
        );
      }
      return;
    }

    // -------------------------------------------------------------
    // VIEW MODE 2: REPRESENTATIVE SAMPLE & NATURAL FREQUENCIES
    // -------------------------------------------------------------
    if (viewMode === 'sample') {
      const topY = b.y + 20;

      // Header Banner
      ctx.textAlign = 'left';
      ctx.font = `700 13px ${font}`;
      ctx.fillStyle = '#38bdf8';
      ctx.fillText(`Thinking With Counts: Representative Cohort of N = ${sampleSize} Individuals`, b.x + 16, topY);

      ctx.font = `400 10px ${font}`;
      ctx.fillStyle = textMuted;
      ctx.fillText('Daniel Kahneman & Amos Tversky discovered that cognitive errors drop from 85% to 0% when framed as natural counts.', b.x + 16, topY + 16);

      // Two Cohort Cards: Hypothesis H vs Alternative ~H
      const cardY = topY + 30;
      const cardW = (b.width - 48) / 2;
      const cardH = 145;

      // Left Card: Hypothesis H (e.g. Librarians)
      const card1X = b.x + 16;
      ctx.fillStyle = 'rgba(16, 185, 129, 0.08)';
      ctx.strokeStyle = '#10b981';
      ctx.lineWidth = 1.5;
      ctx.beginPath();
      if (ctx.roundRect) ctx.roundRect(card1X, cardY, cardW, cardH, [8]);
      else ctx.rect(card1X, cardY, cardW, cardH);
      ctx.fill();
      ctx.stroke();

      ctx.fillStyle = '#10b981';
      ctx.font = `700 12px ${font}`;
      ctx.fillText(preset === 'steve' ? 'LIBRARIANS (Hypothesis H)' : 'HYPOTHESIS H POPULATION', card1X + 16, cardY + 24);

      ctx.font = `800 24px ${font}`;
      ctx.fillText(`${countH} people`, card1X + 16, cardY + 54);

      ctx.font = `500 10px ${font}`;
      ctx.fillStyle = textMuted;
      ctx.fillText(`Total in population: ${(prior * 100).toFixed(1)}% of N = ${sampleSize}`, card1X + 16, cardY + 70);

      // Matching Evidence Subsection
      ctx.fillStyle = 'rgba(16, 185, 129, 0.2)';
      ctx.fillRect(card1X + 16, cardY + 80, cardW - 32, 45);
      ctx.fillStyle = '#10b981';
      ctx.font = `700 14px ${font}`;
      ctx.fillText(`✓ ${countHAndE} Fit Description / Test Positive`, card1X + 26, cardY + 102);
      ctx.font = `500 9.5px ${font}`;
      ctx.fillStyle = '#cbd5e1';
      ctx.fillText(`${(likelihood * 100).toFixed(0)}% likelihood of evidence given H (${countH} × ${(likelihood * 100).toFixed(0)}% = ${countHAndE})`, card1X + 26, cardY + 117);

      // Right Card: Alternative ~H (e.g. Farmers)
      const card2X = card1X + cardW + 16;
      ctx.fillStyle = 'rgba(245, 158, 11, 0.08)';
      ctx.strokeStyle = '#f59e0b';
      ctx.lineWidth = 1.5;
      ctx.beginPath();
      if (ctx.roundRect) ctx.roundRect(card2X, cardY, cardW, cardH, [8]);
      else ctx.rect(card2X, cardY, cardW, cardH);
      ctx.fill();
      ctx.stroke();

      ctx.fillStyle = '#f59e0b';
      ctx.font = `700 12px ${font}`;
      ctx.fillText(preset === 'steve' ? 'FARMERS (Alternative ¬H)' : 'ALTERNATIVE ¬H POPULATION', card2X + 16, cardY + 24);

      ctx.font = `800 24px ${font}`;
      ctx.fillText(`${countNotH} people`, card2X + 16, cardY + 54);

      ctx.font = `500 10px ${font}`;
      ctx.fillStyle = textMuted;
      ctx.fillText(`Total in population: ${(notPrior * 100).toFixed(1)}% of N = ${sampleSize}`, card2X + 16, cardY + 70);

      // Matching Evidence Subsection
      ctx.fillStyle = 'rgba(245, 158, 11, 0.2)';
      ctx.fillRect(card2X + 16, cardY + 80, cardW - 32, 45);
      ctx.fillStyle = '#f59e0b';
      ctx.font = `700 14px ${font}`;
      ctx.fillText(`⚠ ${countNotHAndE} Fit Description / False Alarms`, card2X + 26, cardY + 102);
      ctx.font = `500 9.5px ${font}`;
      ctx.fillStyle = '#cbd5e1';
      ctx.fillText(`${(falsePositive * 100).toFixed(0)}% false alarm rate given ¬H (${countNotH} × ${(falsePositive * 100).toFixed(0)}% = ${countNotHAndE})`, card2X + 26, cardY + 117);

      // Combined Outcome Summary Box
      const summY = cardY + cardH + 16;
      const summW = b.width - 32;
      const summH = 75;

      ctx.fillStyle = 'rgba(30, 41, 59, 0.85)';
      ctx.strokeStyle = '#a855f7';
      ctx.lineWidth = 1.6;
      ctx.beginPath();
      if (ctx.roundRect) ctx.roundRect(b.x + 16, summY, summW, summH, [8]);
      else ctx.rect(b.x + 16, summY, summW, summH);
      ctx.fill();
      ctx.stroke();

      ctx.textAlign = 'left';
      ctx.fillStyle = '#a855f7';
      ctx.font = `700 12px ${font}`;
      ctx.fillText('NATURAL FREQUENCY RATIO CALCULATION:', b.x + 32, summY + 22);

      ctx.fillStyle = '#f8fafc';
      ctx.font = `500 11px ${font}`;
      ctx.fillText(
        `Total individuals who fit the description = ${countHAndE} (from H) + ${countNotHAndE} (from ¬H) = ${countTotalE} people out of ${sampleSize}.`,
        b.x + 32,
        summY + 40
      );

      ctx.fillStyle = '#00d2ff';
      ctx.font = `800 13px monospace`;
      ctx.fillText(
        `Posterior Probability P(H|E) = ${countHAndE} / ${countTotalE} = ${(posterior * 100).toFixed(1)}% (Alternative is ${(posteriorNotH * 100).toFixed(1)}%)`,
        b.x + 32,
        summY + 60
      );
      return;
    }

    // -------------------------------------------------------------
    // VIEW MODE 3: ODDS FORM & BAYES FACTOR SCALE
    // -------------------------------------------------------------
    if (viewMode === 'odds') {
      const topY = b.y + 20;

      ctx.textAlign = 'left';
      ctx.font = `700 13px ${font}`;
      ctx.fillStyle = '#38bdf8';
      ctx.fillText('Bayes\' Theorem in Odds Form: Prior Odds × Bayes Factor = Posterior Odds', b.x + 16, topY);

      ctx.font = `400 10px ${font}`;
      ctx.fillStyle = textMuted;
      ctx.fillText('Thinking in odds turns Bayesian multiplication into straightforward scaling.', b.x + 16, topY + 16);

      // Three Step Cards: Prior Odds -> Bayes Factor -> Posterior Odds
      const cardY = topY + 32;
      const cardW = (b.width - 64) / 3;
      const cardH = 160;

      // Card 1: Prior Odds
      const c1X = b.x + 16;
      ctx.fillStyle = 'rgba(56, 189, 248, 0.08)';
      ctx.strokeStyle = '#38bdf8';
      ctx.lineWidth = 1.5;
      ctx.beginPath();
      if (ctx.roundRect) ctx.roundRect(c1X, cardY, cardW, cardH, [8]);
      else ctx.rect(c1X, cardY, cardW, cardH);
      ctx.fill();
      ctx.stroke();

      ctx.fillStyle = '#38bdf8';
      ctx.font = `700 11px ${font}`;
      ctx.fillText('1. PRIOR ODDS', c1X + 16, cardY + 22);

      ctx.font = `800 22px ${font}`;
      ctx.fillText(priorOdds < 0.1 ? `1 : ${(1 / priorOdds).toFixed(1)}` : priorOdds.toFixed(2), c1X + 16, cardY + 54);

      ctx.font = `500 10px ${font}`;
      ctx.fillStyle = textMuted;
      ctx.fillText(`P(H) = ${(prior * 100).toFixed(1)}%`, c1X + 16, cardY + 74);
      ctx.fillText(`P(¬H) = ${(notPrior * 100).toFixed(1)}%`, c1X + 16, cardY + 90);
      ctx.fillText('O(H) = P(H) / P(¬H)', c1X + 16, cardY + 114);

      // Card 2: Bayes Factor (Likelihood Ratio)
      const c2X = c1X + cardW + 16;
      ctx.fillStyle = 'rgba(168, 85, 247, 0.08)';
      ctx.strokeStyle = '#a855f7';
      ctx.lineWidth = 1.5;
      ctx.beginPath();
      if (ctx.roundRect) ctx.roundRect(c2X, cardY, cardW, cardH, [8]);
      else ctx.rect(c2X, cardY, cardW, cardH);
      ctx.fill();
      ctx.stroke();

      ctx.fillStyle = '#a855f7';
      ctx.font = `700 11px ${font}`;
      ctx.fillText('2. BAYES FACTOR (LR)', c2X + 16, cardY + 22);

      ctx.font = `800 22px ${font}`;
      ctx.fillText(`${bayesFactor.toFixed(2)}×`, c2X + 16, cardY + 54);

      ctx.font = `500 10px ${font}`;
      ctx.fillStyle = textMuted;
      ctx.fillText(`P(E|H) = ${(likelihood * 100).toFixed(0)}%`, c2X + 16, cardY + 74);
      ctx.fillText(`P(E|¬H) = ${(falsePositive * 100).toFixed(0)}%`, c2X + 16, cardY + 90);
      ctx.fillText('BF = P(E|H) / P(E|¬H)', c2X + 16, cardY + 114);
      ctx.fillStyle = '#a855f7';
      ctx.font = `600 9px ${font}`;
      ctx.fillText(evidenceRating, c2X + 16, cardY + 138);

      // Card 3: Posterior Odds
      const c3X = c2X + cardW + 16;
      ctx.fillStyle = 'rgba(16, 185, 129, 0.08)';
      ctx.strokeStyle = '#10b981';
      ctx.lineWidth = 1.5;
      ctx.beginPath();
      if (ctx.roundRect) ctx.roundRect(c3X, cardY, cardW, cardH, [8]);
      else ctx.rect(c3X, cardY, cardW, cardH);
      ctx.fill();
      ctx.stroke();

      ctx.fillStyle = '#10b981';
      ctx.font = `700 11px ${font}`;
      ctx.fillText('3. POSTERIOR ODDS', c3X + 16, cardY + 22);

      ctx.font = `800 22px ${font}`;
      ctx.fillText(posteriorOdds < 0.1 ? `1 : ${(1 / posteriorOdds).toFixed(1)}` : posteriorOdds.toFixed(2), c3X + 16, cardY + 54);

      ctx.font = `500 10px ${font}`;
      ctx.fillStyle = textMuted;
      ctx.fillText(`O(H|E) = O(H) × BF`, c3X + 16, cardY + 74);
      ctx.fillText(`= ${priorOdds.toFixed(3)} × ${bayesFactor.toFixed(2)}`, c3X + 16, cardY + 90);

      ctx.fillStyle = '#10b981';
      ctx.font = `700 12px monospace`;
      ctx.fillText(`P(H|E) = ${(posterior * 100).toFixed(1)}%`, c3X + 16, cardY + 120);

      // Connectors / Multipliers between cards
      ctx.fillStyle = '#f8fafc';
      ctx.font = `800 18px ${font}`;
      ctx.textAlign = 'center';
      ctx.fillText('×', c1X + cardW + 8, cardY + 50);
      ctx.fillText('=', c2X + cardW + 8, cardY + 50);
      return;
    }

    // -------------------------------------------------------------
    // VIEW MODE 4: SEQUENTIAL EVIDENCE UPDATING TRAJECTORY
    // -------------------------------------------------------------
    if (viewMode === 'sequential') {
      const topY = b.y + 15;

      ctx.textAlign = 'left';
      ctx.font = `700 13px ${font}`;
      ctx.fillStyle = '#38bdf8';
      ctx.fillText('Sequential Belief Trajectory: Compounding Evidence Updates Belief Toward Certainty', b.x + 16, topY);

      ctx.font = `400 10px ${font}`;
      ctx.fillStyle = textMuted;
      ctx.fillText('Today\'s posterior becomes tomorrow\'s prior when observing repeated independent tests.', b.x + 16, topY + 16);

      // Trajectory Chart Plot Area
      const chartX = b.x + 45;
      const chartY = topY + 38;
      const chartW = b.width - 70;
      const chartH = b.height - 85;

      // Draw Grid
      ctx.strokeStyle = pal.grid || 'rgba(255, 255, 255, 0.08)';
      ctx.lineWidth = 1.0;
      const yTicks = [0, 0.25, 0.50, 0.75, 1.0];
      yTicks.forEach(yt => {
        const yPix = chartY + chartH - yt * chartH;
        ctx.beginPath();
        ctx.moveTo(chartX, yPix);
        ctx.lineTo(chartX + chartW, yPix);
        ctx.stroke();

        ctx.fillStyle = textMuted;
        ctx.font = `500 10px ${font}`;
        ctx.textAlign = 'right';
        ctx.fillText(`${(yt * 100).toFixed(0)}%`, chartX - 8, yPix + 4);
      });

      // 50% Threshold line
      const y50 = chartY + chartH - 0.5 * chartH;
      ctx.strokeStyle = 'rgba(255, 255, 255, 0.25)';
      ctx.lineWidth = 1.2;
      ctx.setLineDash([3, 3]);
      ctx.beginPath();
      ctx.moveTo(chartX, y50);
      ctx.lineTo(chartX + chartW, y50);
      ctx.stroke();
      ctx.setLineDash([]);
      ctx.fillStyle = textMuted;
      ctx.font = `500 9px ${font}`;
      ctx.textAlign = 'right';
      ctx.fillText('50% Ambiguity Line', chartX + chartW - 8, y50 - 4);

      // Plot Trajectory Steps
      const numSteps = trajectory.length - 1;
      const getX = (stepIdx) => chartX + (stepIdx / numSteps) * chartW;
      const getY = (pVal) => chartY + chartH - pVal * chartH;

      // Fill area under trajectory
      ctx.fillStyle = 'rgba(0, 210, 255, 0.12)';
      ctx.beginPath();
      ctx.moveTo(getX(0), chartY + chartH);
      trajectory.forEach((t, i) => ctx.lineTo(getX(i), getY(t.p)));
      ctx.lineTo(getX(numSteps), chartY + chartH);
      ctx.closePath();
      ctx.fill();

      // Draw Trajectory Line
      ctx.strokeStyle = '#00d2ff';
      ctx.lineWidth = 3.0;
      ctx.beginPath();
      trajectory.forEach((t, i) => {
        if (i === 0) ctx.moveTo(getX(i), getY(t.p));
        else ctx.lineTo(getX(i), getY(t.p));
      });
      ctx.stroke();

      // Draw Step Markers & Callouts
      trajectory.forEach((t, i) => {
        const xPix = getX(i);
        const yPix = getY(t.p);

        // Dot
        ctx.fillStyle = '#00d2ff';
        ctx.strokeStyle = '#ffffff';
        ctx.lineWidth = 2.0;
        ctx.beginPath();
        ctx.arc(xPix, yPix, 6, 0, Math.PI * 2);
        ctx.fill();
        ctx.stroke();

        // Label above dot
        ctx.fillStyle = '#f8fafc';
        ctx.font = `700 11px ${font}`;
        ctx.textAlign = 'center';
        ctx.fillText(`${(t.p * 100).toFixed(1)}%`, xPix, yPix - 12);

        // Step Label below X axis
        ctx.fillStyle = textMuted;
        ctx.font = `600 10px ${font}`;
        ctx.fillText(t.label, xPix, chartY + chartH + 18);
      });
    }
  },

  /**
   * Renders the Love Plot for Covariate Balance Assessment in PSM
   * Shows Absolute Standardized Mean Differences before and after matching
   * with vertical benchmark cutoffs at |SMD| = 0.10 and |SMD| = 0.05.
   */
  renderLovePlot(engine, balanceData, options = {}) {
    if (!engine) return;
    engine.lastRenderFn = () => this.renderLovePlot(engine, balanceData, options);
    engine.clear();

    const b = engine.getPlotBounds ? engine.getPlotBounds() : engine.getBounds();
    const ctx = engine.ctx;
    const font = '-apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif';
    const isDark = engine.options.theme !== 'light';

    const textMain = isDark ? '#f8fafc' : '#0f172a';
    const textMuted = isDark ? '#94a3b8' : '#64748b';
    const gridColor = isDark ? 'rgba(255, 255, 255, 0.08)' : 'rgba(0, 0, 0, 0.08)';

    const table = balanceData?.balanceTable || [];
    if (table.length === 0) {
      ctx.fillStyle = textMuted;
      ctx.font = `500 13px ${font}`;
      ctx.textAlign = 'center';
      ctx.fillText('No covariate balance data available for Love Plot', b.x + b.width / 2, b.y + b.height / 2);
      return;
    }

    const title = options.title || 'Love Plot: Covariate Balance (Pre- vs. Post-Matching)';
    const maxVal = Math.max(0.40, ...table.map(d => Math.max(d.absSmdPre || 0, d.absSmdPost || 0))) * 1.15;

    // Margins
    const labelWidth = Math.min(140, b.width * 0.28);
    const plotX = b.x + labelWidth;
    const plotW = b.width - labelWidth - 20;
    const plotY = b.y + 35;
    const plotH = b.height - 65;

    const xToPix = (smd) => plotX + (Math.max(0, smd) / maxVal) * plotW;
    const k = table.length;
    const yToPix = (idx) => plotY + ((idx + 0.6) / (k + 0.2)) * plotH;

    // Header Title
    ctx.fillStyle = textMain;
    ctx.font = `700 13px ${font}`;
    ctx.textAlign = 'left';
    ctx.fillText(title, b.x + 10, b.y + 18);

    // Legend
    ctx.font = `600 11px ${font}`;
    ctx.textAlign = 'right';
    const legX = b.x + b.width - 20;
    // Post match dot
    ctx.fillStyle = '#10b981';
    ctx.beginPath(); ctx.arc(legX - 165, b.y + 14, 5, 0, Math.PI * 2); ctx.fill();
    ctx.fillStyle = textMain;
    ctx.fillText('Matched (Post)', legX - 85, b.y + 18);
    // Pre match dot
    ctx.fillStyle = '#f43f5e';
    ctx.beginPath(); ctx.arc(legX - 70, b.y + 14, 5, 0, Math.PI * 2); ctx.fill();
    ctx.fillStyle = textMain;
    ctx.fillText('Unadjusted (Pre)', legX, b.y + 18);

    // X-Axis Grid & Ticks
    const xSteps = [0.0, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.80, 1.0].filter(v => v <= maxVal);
    ctx.lineWidth = 1.0;
    xSteps.forEach(val => {
      const px = xToPix(val);
      ctx.strokeStyle = gridColor;
      ctx.beginPath();
      ctx.moveTo(px, plotY);
      ctx.lineTo(px, plotY + plotH);
      ctx.stroke();

      ctx.fillStyle = textMuted;
      ctx.font = `500 10px ${font}`;
      ctx.textAlign = 'center';
      ctx.fillText(val.toFixed(2), px, plotY + plotH + 16);
    });

    // X Axis Label
    ctx.fillStyle = textMain;
    ctx.font = `600 11px ${font}`;
    ctx.textAlign = 'center';
    ctx.fillText('Absolute Standardized Mean Difference (|SMD|)', plotX + plotW / 2, plotY + plotH + 30);

    // Cutoff 0.10 Line (Standard Benchmark)
    const cut10X = xToPix(0.10);
    ctx.strokeStyle = '#00d2ff';
    ctx.lineWidth = 1.8;
    ctx.setLineDash([5, 4]);
    ctx.beginPath();
    ctx.moveTo(cut10X, plotY);
    ctx.lineTo(cut10X, plotY + plotH);
    ctx.stroke();

    // Cutoff 0.05 Line (Strict Benchmark)
    if (maxVal >= 0.05) {
      const cut05X = xToPix(0.05);
      ctx.strokeStyle = '#a855f7';
      ctx.lineWidth = 1.2;
      ctx.setLineDash([3, 3]);
      ctx.beginPath();
      ctx.moveTo(cut05X, plotY);
      ctx.lineTo(cut05X, plotY + plotH);
      ctx.stroke();
      ctx.setLineDash([]);
    } else {
      ctx.setLineDash([]);
    }

    // Benchmark Cutoff Annotations at top
    ctx.fillStyle = '#00d2ff';
    ctx.font = `700 9px ${font}`;
    ctx.textAlign = 'center';
    ctx.fillText('0.10 Threshold', cut10X, plotY - 4);

    // Draw Rows for each Covariate
    table.forEach((row, idx) => {
      const py = yToPix(idx);
      const preX = xToPix(row.absSmdPre);
      const postX = xToPix(row.absSmdPost);

      // Horizontal subtle track line
      ctx.strokeStyle = gridColor;
      ctx.lineWidth = 1;
      ctx.beginPath();
      ctx.moveTo(plotX, py);
      ctx.lineTo(plotX + plotW, py);
      ctx.stroke();

      // Connecting arrow/segment from pre to post
      ctx.strokeStyle = row.absSmdPost <= 0.10 ? 'rgba(16, 185, 129, 0.45)' : 'rgba(244, 63, 94, 0.45)';
      ctx.lineWidth = 2.0;
      ctx.beginPath();
      ctx.moveTo(preX, py);
      ctx.lineTo(postX, py);
      ctx.stroke();

      // Pre-matching Marker (Rose)
      ctx.fillStyle = '#f43f5e';
      ctx.strokeStyle = isDark ? '#0f172a' : '#ffffff';
      ctx.lineWidth = 1.5;
      ctx.beginPath();
      ctx.arc(preX, py, 5, 0, Math.PI * 2);
      ctx.fill();
      ctx.stroke();

      // Post-matching Marker (Emerald)
      ctx.fillStyle = '#10b981';
      ctx.beginPath();
      ctx.arc(postX, py, 6, 0, Math.PI * 2);
      ctx.fill();
      ctx.stroke();

      // Covariate Name Label
      ctx.fillStyle = textMain;
      ctx.font = `600 11px ${font}`;
      ctx.textAlign = 'right';
      ctx.fillText(row.covariate, plotX - 10, py + 4);

      // Delta improvement text
      if (row.percentReduction !== undefined) {
        ctx.fillStyle = row.absSmdPost <= 0.10 ? '#10b981' : '#f43f5e';
        ctx.font = `500 9px ${font}`;
        ctx.textAlign = 'left';
        ctx.fillText(`${(row.absSmdPost).toFixed(3)}`, Math.max(preX, postX) + 10, py + 3);
      }
    });
  },

  /**
   * Renders Propensity Score Overlap & Common Support Plot
   */
  renderPsmOverlapPlot(engine, overlapData, options = {}) {
    if (!engine) return;
    engine.lastRenderFn = () => this.renderPsmOverlapPlot(engine, overlapData, options);
    engine.clear();

    const b = engine.getPlotBounds ? engine.getPlotBounds() : engine.getBounds();
    const ctx = engine.ctx;
    const font = '-apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif';
    const isDark = engine.options.theme !== 'light';

    const textMain = isDark ? '#f8fafc' : '#0f172a';
    const textMuted = isDark ? '#94a3b8' : '#64748b';
    const gridColor = isDark ? 'rgba(255, 255, 255, 0.08)' : 'rgba(0, 0, 0, 0.08)';

    const bins = overlapData?.bins || [];
    if (bins.length === 0) return;

    const title = options.title || 'Propensity Score Distribution & Common Support Overlap';
    const cs = overlapData?.commonSupport || { min: 0, max: 1 };

    const plotX = b.x + 35;
    const plotW = b.width - 60;
    const plotY = b.y + 40;
    const plotH = b.height - 75;

    // Header Title
    ctx.fillStyle = textMain;
    ctx.font = `700 13px ${font}`;
    ctx.textAlign = 'left';
    ctx.fillText(title, b.x + 10, b.y + 18);

    // Legend
    ctx.font = `600 11px ${font}`;
    ctx.textAlign = 'right';
    const legX = b.x + b.width - 20;

    // Treated legend
    ctx.fillStyle = '#00d2ff';
    ctx.fillRect(legX - 180, b.y + 10, 12, 8);
    ctx.fillStyle = textMain;
    ctx.fillText('Treated Group', legX - 95, b.y + 18);

    // Control legend
    ctx.fillStyle = '#f59e0b';
    ctx.fillRect(legX - 85, b.y + 10, 12, 8);
    ctx.fillStyle = textMain;
    ctx.fillText('Control Group', legX, b.y + 18);

    // Find max proportion for Y scaling
    let maxProp = 0.05;
    bins.forEach(bin => {
      maxProp = Math.max(maxProp, bin.preTreated, bin.preControl, bin.postTreated, bin.postControl);
    });
    maxProp *= 1.25;

    const xToPix = (p) => plotX + Math.max(0, Math.min(1, p)) * plotW;
    const yToPix = (prop) => plotY + plotH - (prop / maxProp) * plotH;

    // Common Support Shaded Area
    const csMinX = xToPix(cs.min);
    const csMaxX = xToPix(cs.max);
    ctx.fillStyle = isDark ? 'rgba(0, 210, 255, 0.07)' : 'rgba(0, 210, 255, 0.12)';
    ctx.fillRect(csMinX, plotY, csMaxX - csMinX, plotH);

    // Common Support Vertical Dashed Lines
    ctx.strokeStyle = '#00d2ff';
    ctx.lineWidth = 1.5;
    ctx.setLineDash([4, 4]);
    ctx.beginPath();
    ctx.moveTo(csMinX, plotY);
    ctx.lineTo(csMinX, plotY + plotH);
    ctx.moveTo(csMaxX, plotY);
    ctx.lineTo(csMaxX, plotY + plotH);
    ctx.stroke();
    ctx.setLineDash([]);

    // Common Support Label
    ctx.fillStyle = '#00d2ff';
    ctx.font = `600 10px ${font}`;
    ctx.textAlign = 'center';
    ctx.fillText(`Common Support [${cs.min.toFixed(2)}, ${cs.max.toFixed(2)}]`, (csMinX + csMaxX) / 2, plotY + 14);

    // X Axis Ticks & Grid
    [0.0, 0.2, 0.4, 0.6, 0.8, 1.0].forEach(p => {
      const px = xToPix(p);
      ctx.strokeStyle = gridColor;
      ctx.lineWidth = 1;
      ctx.beginPath();
      ctx.moveTo(px, plotY);
      ctx.lineTo(px, plotY + plotH);
      ctx.stroke();

      ctx.fillStyle = textMuted;
      ctx.font = `500 10px ${font}`;
      ctx.textAlign = 'center';
      ctx.fillText(p.toFixed(1), px, plotY + plotH + 16);
    });

    ctx.fillStyle = textMain;
    ctx.font = `600 11px ${font}`;
    ctx.textAlign = 'center';
    ctx.fillText('Estimated Propensity Score e(X)', plotX + plotW / 2, plotY + plotH + 30);

    // Draw Treated Distribution (Pre-match, Cyan Area)
    ctx.fillStyle = 'rgba(0, 210, 255, 0.25)';
    ctx.strokeStyle = '#00d2ff';
    ctx.lineWidth = 2.0;
    ctx.beginPath();
    ctx.moveTo(xToPix(bins[0].mid), plotY + plotH);
    bins.forEach(bin => {
      ctx.lineTo(xToPix(bin.mid), yToPix(bin.preTreated));
    });
    ctx.lineTo(xToPix(bins[bins.length - 1].mid), plotY + plotH);
    ctx.closePath();
    ctx.fill();
    ctx.stroke();

    // Draw Control Distribution (Pre-match, Amber Area)
    ctx.fillStyle = 'rgba(245, 158, 11, 0.25)';
    ctx.strokeStyle = '#f59e0b';
    ctx.lineWidth = 2.0;
    ctx.beginPath();
    ctx.moveTo(xToPix(bins[0].mid), plotY + plotH);
    bins.forEach(bin => {
      ctx.lineTo(xToPix(bin.mid), yToPix(bin.preControl));
    });
    ctx.lineTo(xToPix(bins[bins.length - 1].mid), plotY + plotH);
    ctx.closePath();
    ctx.fill();
    ctx.stroke();

    // Draw Matched Treated (Dashed Cyan Line)
    ctx.strokeStyle = '#00d2ff';
    ctx.lineWidth = 2.0;
    ctx.setLineDash([4, 3]);
    ctx.beginPath();
    bins.forEach((bin, i) => {
      const px = xToPix(bin.mid);
      const py = yToPix(bin.postTreated);
      if (i === 0) ctx.moveTo(px, py);
      else ctx.lineTo(px, py);
    });
    ctx.stroke();

    // Draw Matched Control (Dashed Emerald Line)
    ctx.strokeStyle = '#10b981';
    ctx.lineWidth = 2.2;
    ctx.setLineDash([2, 2]);
    ctx.beginPath();
    bins.forEach((bin, i) => {
      const px = xToPix(bin.mid);
      const py = yToPix(bin.postControl);
      if (i === 0) ctx.moveTo(px, py);
      else ctx.lineTo(px, py);
    });
    ctx.stroke();
    ctx.setLineDash([]);
  }
};


