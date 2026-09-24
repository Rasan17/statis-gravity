/**
 * Statis-Gravity - High-DPI Canvas Chart Engine
 * Lightweight, hardware-accelerated, zero-dependency data visualization engine.
 */

export class ChartEngine {
  constructor(canvasId, options = {}) {
    this.canvas = typeof canvasId === 'string' ? document.getElementById(canvasId) : canvasId;
    if (!this.canvas) return;
    this.ctx = this.canvas.getContext('2d');
    this.options = Object.assign({
      padding: { top: 40, right: 30, bottom: 50, left: 60 },
      theme: 'dark',
      fontFamily: '-apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif'
    }, options);

    this.colors = {
      dark: {
        bg: '#0f172a',
        grid: 'rgba(255, 255, 255, 0.08)',
        axis: 'rgba(255, 255, 255, 0.25)',
        text: '#94a3b8',
        textBold: '#f8fafc',
        primary: '#00d2ff',
        primaryGlow: 'rgba(0, 210, 255, 0.25)',
        secondary: '#a855f7',
        accent: '#10b981',
        warning: '#f59e0b',
        danger: '#ef4444'
      },
      light: {
        bg: '#ffffff',
        grid: 'rgba(0, 0, 0, 0.06)',
        axis: 'rgba(0, 0, 0, 0.2)',
        text: '#64748b',
        textBold: '#0f172a',
        primary: '#0284c7',
        primaryGlow: 'rgba(2, 132, 199, 0.2)',
        secondary: '#9333ea',
        accent: '#059669',
        warning: '#d97706',
        danger: '#dc2626'
      }
    };

    this.initHiDPI();
    this.setupResizeObserver();
  }

  get palette() {
    return this.colors[this.options.theme] || this.colors.dark;
  }

  initHiDPI() {
    if (!this.canvas) return;
    const rect = this.canvas.getBoundingClientRect();
    const dpr = window.devicePixelRatio || 1;
    this.width = rect.width || this.canvas.width || 600;
    this.height = rect.height || this.canvas.height || 360;

    this.canvas.width = this.width * dpr;
    this.canvas.height = this.height * dpr;
    this.ctx.resetTransform?.();
    this.ctx.scale(dpr, dpr);
  }

  setupResizeObserver() {
    if (typeof ResizeObserver === 'function' && this.canvas) {
      this.resizeObserver = new ResizeObserver(() => {
        this.initHiDPI();
        if (this.lastRenderFn) this.lastRenderFn();
      });
      this.resizeObserver.observe(this.canvas);
    }
  }

  clear() {
    this.ctx.clearRect(0, 0, this.width, this.height);
  }

  getPlotBounds() {
    const p = this.options.padding;
    return {
      x: p.left,
      y: p.top,
      width: Math.max(10, this.width - p.left - p.right),
      height: Math.max(10, this.height - p.top - p.bottom)
    };
  }

  getBounds() {
    return this.getPlotBounds();
  }

  drawAxes({ xLabel = '', yLabel = '', xTicks = [], yTicks = [], title = '' } = {}) {
    const b = this.getPlotBounds();
    const ctx = this.ctx;
    const pal = this.palette;

    // Title
    if (title) {
      ctx.fillStyle = pal.textBold;
      ctx.font = `600 14px ${this.options.fontFamily}`;
      ctx.textAlign = 'left';
      ctx.fillText(title, b.x, b.y - 16);
    }

    // Grid lines & Y ticks
    ctx.lineWidth = 1;
    for (const tick of yTicks) {
      const yPos = b.y + b.height - (tick.norm * b.height);
      ctx.strokeStyle = pal.grid;
      ctx.beginPath();
      ctx.moveTo(b.x, yPos);
      ctx.lineTo(b.x + b.width, yPos);
      ctx.stroke();

      ctx.fillStyle = pal.text;
      ctx.font = `500 11px ${this.options.fontFamily}`;
      ctx.textAlign = 'right';
      ctx.textBaseline = 'middle';
      ctx.fillText(tick.label, b.x - 10, yPos);
    }

    // X ticks
    for (const tick of xTicks) {
      const xPos = b.x + (tick.norm * b.width);
      ctx.strokeStyle = pal.grid;
      ctx.beginPath();
      ctx.moveTo(xPos, b.y);
      ctx.lineTo(xPos, b.y + b.height);
      ctx.stroke();

      ctx.fillStyle = pal.text;
      ctx.font = `500 11px ${this.options.fontFamily}`;
      ctx.textAlign = 'center';
      ctx.textBaseline = 'top';
      ctx.fillText(tick.label, xPos, b.y + b.height + 8);
    }

    // Axis frames
    ctx.strokeStyle = pal.axis;
    ctx.lineWidth = 1.5;
    ctx.beginPath();
    ctx.moveTo(b.x, b.y);
    ctx.lineTo(b.x, b.y + b.height);
    ctx.lineTo(b.x + b.width, b.y + b.height);
    ctx.stroke();

    // Axis Labels
    ctx.fillStyle = pal.text;
    ctx.font = `600 11px ${this.options.fontFamily}`;
    if (xLabel) {
      ctx.textAlign = 'center';
      ctx.fillText(xLabel, b.x + b.width / 2, b.y + b.height + 34);
    }
    if (yLabel) {
      ctx.save();
      ctx.translate(b.x - 42, b.y + b.height / 2);
      ctx.rotate(-Math.PI / 2);
      ctx.textAlign = 'center';
      ctx.fillText(yLabel, 0, 0);
      ctx.restore();
    }
  }

  saveAsImage(fileName = 'statis-gravity-plot.png') {
    if (!this.canvas) return;
    const link = document.createElement('a');
    link.download = fileName;
    link.href = this.canvas.toDataURL('image/png');
    link.click();
  }
}
