/**
 * Statis-Gravity - Statistical Plots Suite
 * Renders Box-and-Whisker, Distribution Histograms, Scatter + Regression CI, and ROC Curves.
 */

import { ChartEngine } from './chart-engine.js';
import { Distributions } from '../stats/distributions.js';

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
  }
};
