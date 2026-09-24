/**
 * Statis-Gravity - Main Application Controller
 * Manages tab switching, reactive statistical computations, dynamic charts, and report exports.
 */

import { Descriptive } from './stats/descriptive.js';
import { Hypothesis } from './stats/hypothesis.js';
import { Anova } from './stats/anova.js';
import { Categorical } from './stats/categorical.js';
import { Correlation } from './stats/correlation.js';
import { Diagnostic } from './stats/diagnostic.js';
import { PowerAnalysis } from './stats/power.js';
import { Teaching } from './stats/teaching.js';

import { ChartEngine } from './visualization/chart-engine.js';
import { Plots } from './visualization/plots.js';

import { DataParser } from './data/parser.js';
import { Exporter } from './data/exporter.js';

class StatisGravityApp {
  constructor() {
    this.currentTheme = localStorage.getItem('sg_theme') || 'dark';
    this.engines = {};
    this.init();
  }

  init() {
    this.applyTheme(this.currentTheme);
    this.setupTabs();
    this.setupEventListeners();
    this.initChartEngines();
    this.loadInitialSamples();
  }

  applyTheme(theme) {
    this.currentTheme = theme;
    document.documentElement.setAttribute('data-theme', theme);
    localStorage.setItem('sg_theme', theme);
    for (const key in this.engines) {
      if (this.engines[key]) {
        this.engines[key].options.theme = theme;
        if (this.engines[key].lastRenderFn) {
          this.engines[key].lastRenderFn();
        }
      }
    }
  }

  setupTabs() {
    const tabs = document.querySelectorAll('.tab-btn');
    const panes = document.querySelectorAll('.tab-pane');

    tabs.forEach(tab => {
      tab.addEventListener('click', () => {
        const targetId = tab.dataset.target;
        tabs.forEach(t => t.classList.remove('active'));
        panes.forEach(p => p.classList.remove('active'));

        tab.classList.add('active');
        const targetPane = document.getElementById(targetId);
        if (targetPane) {
          targetPane.classList.add('active');
          // Re-render all charts in newly visible tab
          const canvases = targetPane.querySelectorAll('canvas');
          canvases.forEach(canvas => {
            if (canvas && this.engines[canvas.id]?.lastRenderFn) {
              this.engines[canvas.id].initHiDPI();
              this.engines[canvas.id].lastRenderFn();
            }
          });
        }
      });
    });
  }

  initChartEngines() {
    const canvasIds = [
      'descCanvas',
      'descBoxCanvas',
      'descViolinCanvas',
      'hypoCanvas',
      'anovaCanvas',
      'corrCanvas',
      'rocCanvas',
      'teachingDistCanvas',
      'teachingCltParentCanvas',
      'teachingCltSamplingCanvas'
    ];

    canvasIds.forEach(id => {
      const el = document.getElementById(id);
      if (el) {
        this.engines[id] = new ChartEngine(el, { theme: this.currentTheme });
      }
    });
  }

  setupEventListeners() {
    // Theme toggle
    document.getElementById('themeToggleBtn')?.addEventListener('click', () => {
      this.applyTheme(this.currentTheme === 'dark' ? 'light' : 'dark');
    });

    // 1. Descriptive Stats Events
    document.getElementById('descComputeBtn')?.addEventListener('click', () => this.runDescriptive());
    document.getElementById('descSampleBtn')?.addEventListener('click', () => {
      const data = DataParser.samples.icpDynamics.groupA.data;
      document.getElementById('descInput').value = data.join(', ');
      this.runDescriptive();
    });
    document.getElementById('descOutlierSampleBtn')?.addEventListener('click', () => {
      const data = DataParser.samples.drainOutputSkewed.data;
      document.getElementById('descInput').value = data.join(', ');
      this.runDescriptive();
    });

    // 2. Hypothesis Testing Events
    document.getElementById('hypoComputeBtn')?.addEventListener('click', () => this.runHypothesis());
    document.getElementById('hypoSampleBtn')?.addEventListener('click', () => {
      document.getElementById('hypoGroupA').value = DataParser.samples.icpDynamics.groupA.data.join(', ');
      document.getElementById('hypoGroupB').value = DataParser.samples.icpDynamics.groupB.data.join(', ');
      this.runHypothesis();
    });
    document.getElementById('hypoErrorBarMode')?.addEventListener('change', () => {
      this.runHypothesis();
    });

    // 3. ANOVA Events
    document.getElementById('anovaComputeBtn')?.addEventListener('click', () => this.runAnova());
    document.getElementById('anovaSampleBtn')?.addEventListener('click', () => {
      const g = DataParser.samples.cranialAsymmetry.groups;
      document.getElementById('anovaG1').value = g[0].data.join(', ');
      document.getElementById('anovaG2').value = g[1].data.join(', ');
      document.getElementById('anovaG3').value = g[2].data.join(', ');
      this.runAnova();
    });
    document.getElementById('anovaErrorBarMode')?.addEventListener('change', () => {
      this.runAnova();
    });

    // 4. Categorical / 2x2 Risk & Contingency
    const matrixInputs = ['catA', 'catB', 'catC', 'catD'];
    matrixInputs.forEach(id => {
      document.getElementById(id)?.addEventListener('input', () => this.runCategorical());
    });
    document.getElementById('catAnalysisMode')?.addEventListener('change', (e) => {
      const mode = e.target.value;
      const curA = parseFloat(document.getElementById('catA')?.value) || 0;
      const curB = parseFloat(document.getElementById('catB')?.value) || 0;
      const curC = parseFloat(document.getElementById('catC')?.value) || 0;
      const curD = parseFloat(document.getElementById('catD')?.value) || 0;

      // Auto-swap default sample when user switches frameworks on unmodified sample data
      if (mode === 'diagnostic' && curA === 14 && curB === 36 && curC === 186 && curD === 164) {
        const s = DataParser.samples.diagnostic2x2;
        document.getElementById('catA').value = s.a;
        document.getElementById('catB').value = s.b;
        document.getElementById('catC').value = s.c;
        document.getElementById('catD').value = s.d;
      } else if (mode === 'study' && curA === 92 && curB === 8 && curC === 12 && curD === 188) {
        const s = DataParser.samples.shuntComplications || DataParser.samples.shunt2x2;
        document.getElementById('catA').value = s.a;
        document.getElementById('catB').value = s.b;
        document.getElementById('catC').value = s.c;
        document.getElementById('catD').value = s.d;
      }
      this.runCategorical();
    });
    document.getElementById('catSampleBtn')?.addEventListener('click', () => {
      const mode = document.getElementById('catAnalysisMode')?.value || 'diagnostic';
      if (mode === 'diagnostic') {
        const s = DataParser.samples.diagnostic2x2;
        document.getElementById('catA').value = s.a;
        document.getElementById('catB').value = s.b;
        document.getElementById('catC').value = s.c;
        document.getElementById('catD').value = s.d;
      } else {
        const s = DataParser.samples.shuntComplications || DataParser.samples.shunt2x2;
        document.getElementById('catA').value = s.a;
        document.getElementById('catB').value = s.b;
        document.getElementById('catC').value = s.c;
        document.getElementById('catD').value = s.d;
      }
      this.runCategorical();
    });

    // 5. Correlation & Regression Events
    document.getElementById('corrComputeBtn')?.addEventListener('click', () => this.runCorrelation());
    document.getElementById('corrSampleBtn')?.addEventListener('click', () => {
      document.getElementById('corrX').value = '10, 12, 14, 15, 18, 20, 22, 24, 25, 28, 30, 32';
      document.getElementById('corrY').value = '15, 19, 21, 22, 27, 31, 33, 37, 39, 42, 47, 49';
      this.runCorrelation();
    });
    document.getElementById('corrUShapeBtn')?.addEventListener('click', () => {
      document.getElementById('corrX').value = '42, 46, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110';
      document.getElementById('corrY').value = '82, 68, 54, 38, 25, 18, 16, 19, 26, 36, 49, 62, 73, 85, 96';
      this.runCorrelation();
    });
    document.getElementById('corrNCurveBtn')?.addEventListener('click', () => {
      document.getElementById('corrX').value = '1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15';
      document.getElementById('corrY').value = '15, 28, 45, 62, 70, 65, 50, 36, 25, 20, 26, 42, 60, 78, 95';
      this.runCorrelation();
    });

    // 6. Diagnostic & ROC Events
    document.getElementById('rocComputeBtn')?.addEventListener('click', () => this.runDiagnostic());
    document.getElementById('rocSampleBtn')?.addEventListener('click', () => {
      const samples = DataParser.samples.diagnosticBiomarker.data;
      const lines = samples.map(s => `${s.score}, ${s.status}`).join('\n');
      document.getElementById('rocInput').value = lines;
      this.runDiagnostic();
    });

    // 7. Power Analysis Events
    const pwrInputs = [
      'pwrM1', 'pwrM2', 'pwrSD',
      'pwrPairedM1', 'pwrPairedM2', 'pwrPairedSD',
      'pwrP1', 'pwrP2', 'pwrAlpha', 'pwrPower', 'pwrGivenN', 'pwrPropTestType'
    ];
    pwrInputs.forEach(id => {
      document.getElementById(id)?.addEventListener('input', () => this.runPower());
      document.getElementById(id)?.addEventListener('change', () => this.runPower());
    });

    document.getElementById('pwrStudyDesign')?.addEventListener('change', (e) => {
      const design = e.target.value;
      const indepSec = document.getElementById('pwrSectionIndep');
      const pairedSec = document.getElementById('pwrSectionPaired');
      const contSec = document.getElementById('pwrSectionContingency');
      const cardFisher = document.getElementById('pwrCardFisher');
      const cardRisk = document.getElementById('pwrCardClinicalRisk');
      const subtitle = document.getElementById('pwrSubtitle');
      const givenNLabel = document.getElementById('pwrGivenNLabel');

      if (indepSec) indepSec.style.display = design === 'independent' ? 'block' : 'none';
      if (pairedSec) pairedSec.style.display = design === 'paired' ? 'block' : 'none';
      if (contSec) contSec.style.display = design === 'contingency' ? 'block' : 'none';

      if (cardFisher) cardFisher.style.display = design === 'contingency' ? 'block' : 'none';
      if (cardRisk) cardRisk.style.display = design === 'contingency' ? 'block' : 'none';

      if (subtitle) {
        if (design === 'independent') subtitle.innerText = 'Two Independent Groups (Independent Samples t-Test)';
        else if (design === 'paired') subtitle.innerText = 'Paired Mean Study Design (Before vs After / Paired t-Test)';
        else subtitle.innerText = 'Clinical Study: Test vs Control (2x2 Contingency: Chi-Square / Fisher\'s Exact)';
      }

      if (givenNLabel) {
        givenNLabel.innerText = design === 'paired'
          ? 'Available Number of Pairs (N):'
          : 'Available Sample Size per Group (n):';
      }

      this.runPower();
    });

    document.getElementById('pwrCalcGoal')?.addEventListener('change', (e) => {
      const goal = e.target.value;
      const targetGroup = document.getElementById('pwrTargetPowerGroup');
      const givenGroup = document.getElementById('pwrGivenNGroup');
      if (targetGroup) targetGroup.style.display = goal === 'sample_size' ? 'block' : 'none';
      if (givenGroup) givenGroup.style.display = goal === 'power' ? 'block' : 'none';
      this.runPower();
    });

    document.getElementById('powerComputeBtn')?.addEventListener('click', () => this.runPower());

    document.getElementById('powerSampleBtn')?.addEventListener('click', () => {
      const design = document.getElementById('pwrStudyDesign')?.value || 'independent';
      if (design === 'independent') {
        document.getElementById('pwrM1').value = '10';
        document.getElementById('pwrM2').value = '15';
        document.getElementById('pwrSD').value = '10';
      } else if (design === 'paired') {
        document.getElementById('pwrPairedM1').value = '120';
        document.getElementById('pwrPairedM2').value = '112';
        document.getElementById('pwrPairedSD').value = '10';
      } else if (design === 'contingency') {
        document.getElementById('pwrP1').value = '0.08';
        document.getElementById('pwrP2').value = '0.20';
        document.getElementById('pwrPropTestType').value = 'fisher';
      }
      this.runPower();
    });

    // Export & Copy Buttons
    document.querySelectorAll('.btn-copy-report').forEach(btn => {
      btn.addEventListener('click', async (e) => {
        const targetId = e.currentTarget.dataset.reportId;
        const text = document.getElementById(targetId)?.innerText;
        if (text) {
          await Exporter.copyToClipboard(text);
          const orig = e.currentTarget.innerText;
          e.currentTarget.innerText = 'Copied!';
          setTimeout(() => { e.currentTarget.innerText = orig; }, 1500);
        }
      });
    });

    document.querySelectorAll('.btn-export-docx').forEach(btn => {
      btn.addEventListener('click', (e) => {
        const tab = e.currentTarget.dataset.tab;
        this.exportDocx(tab);
      });
    });

    document.querySelectorAll('.btn-save-plot').forEach(btn => {
      btn.addEventListener('click', (e) => {
        const canvasId = e.currentTarget.dataset.canvasId;
        if (this.engines[canvasId]) {
          this.engines[canvasId].saveAsImage(`${canvasId}.png`);
        }
      });
    });

    // 8. Teaching & Simulation Events
    document.getElementById('teachingDistSelect')?.addEventListener('change', () => {
      this.renderTeachingParams();
      this.runTeachingDistribution();
    });

    document.getElementById('teachingNRange')?.addEventListener('input', (e) => {
      const val = e.target.value;
      const el = document.getElementById('teachingNVal');
      if (el) el.innerText = val;
      this.runTeachingDistribution();
    });

    document.getElementById('teachingGenBtn')?.addEventListener('click', () => {
      this.runTeachingDistribution();
    });

    document.getElementById('teachingCopyBtn')?.addEventListener('click', (e) => {
      if (this.currentTeachingData) {
        navigator.clipboard.writeText(this.currentTeachingData.join(', '));
        const orig = e.currentTarget.innerText;
        e.currentTarget.innerText = 'Copied!';
        setTimeout(() => { e.currentTarget.innerText = orig; }, 1500);
      }
    });

    document.getElementById('teachingSendToDescBtn')?.addEventListener('click', () => {
      if (this.currentTeachingData) {
        const descInput = document.getElementById('descInput');
        if (descInput) {
          descInput.value = this.currentTeachingData.join(', ');
        }
        const descTabBtn = document.querySelector('.tab-btn[data-target="tab-descriptive"]');
        if (descTabBtn) descTabBtn.click();
        this.runDescriptive();
      }
    });

    // CLT Simulation Controls
    document.getElementById('cltPopSelect')?.addEventListener('change', (e) => {
      Teaching.clt.setPopulation(e.target.value);
      this.updateCltUI(Teaching.clt.getSummary());
    });

    document.getElementById('cltNRange')?.addEventListener('input', (e) => {
      const n = parseInt(e.target.value);
      const el = document.getElementById('cltNVal');
      if (el) el.innerText = n;
      Teaching.clt.setSampleSize(n);
      this.updateCltUI(Teaching.clt.getSummary());
    });

    document.getElementById('cltStepBtn')?.addEventListener('click', () => {
      const summary = Teaching.clt.drawSamples(1);
      this.updateCltUI(summary);
    });

    document.getElementById('cltDraw100Btn')?.addEventListener('click', () => {
      const summary = Teaching.clt.drawSamples(100);
      this.updateCltUI(summary);
    });

    document.getElementById('cltDraw1000Btn')?.addEventListener('click', () => {
      const summary = Teaching.clt.drawSamples(1000);
      this.updateCltUI(summary);
    });

    document.getElementById('cltResetBtn')?.addEventListener('click', () => {
      Teaching.clt.reset();
      this.updateCltUI(Teaching.clt.getSummary());
    });

    // Disclaimer Modal Dismissal
    document.getElementById('disclaimerDismissBtn')?.addEventListener('click', () => {
      document.getElementById('disclaimerModal')?.classList.add('hidden');
    });
    document.getElementById('disclaimerModal')?.addEventListener('click', (e) => {
      if (e.target.id === 'disclaimerModal') {
        e.currentTarget.classList.add('hidden');
      }
    });
    document.addEventListener('keydown', (e) => {
      if (e.key === 'Escape') {
        document.getElementById('disclaimerModal')?.classList.add('hidden');
      }
    });
  }

  loadInitialSamples() {
    // Populate all initial views so app has immediate visual impact
    document.getElementById('descSampleBtn')?.click();
    document.getElementById('hypoSampleBtn')?.click();
    document.getElementById('anovaSampleBtn')?.click();
    document.getElementById('catSampleBtn')?.click();
    document.getElementById('corrSampleBtn')?.click();
    document.getElementById('rocSampleBtn')?.click();
    this.runPower();
    this.initTeachingModule();
  }

  // --- Run Modules ---

  runDescriptive() {
    const raw = document.getElementById('descInput')?.value || '';
    const series = DataParser.parseSeries(raw);
    const stats = Descriptive.calculate(series);

    if (stats.error) {
      alert(stats.error);
      return;
    }

    document.getElementById('descN').innerText = stats.n;
    const rangeEl = document.getElementById('descRange');
    if (rangeEl) rangeEl.innerText = `Range: ${stats.min.toFixed(1)} to ${stats.max.toFixed(1)}`;

    document.getElementById('descMean').innerText = stats.mean.toFixed(2);
    const sdSubEl = document.getElementById('descSDSub');
    if (sdSubEl) sdSubEl.innerText = `SD: ±${stats.sd.toFixed(2)}`;

    document.getElementById('descMedian').innerText = stats.median.toFixed(2);

    // Mode
    const modeEl = document.getElementById('descMode');
    const modeSubEl = document.getElementById('descModeSub');
    if (modeEl) {
      if (stats.modes && stats.modes.length > 0) {
        if (stats.modes.length === 1) {
          modeEl.innerText = stats.modes[0].toFixed(2);
          if (modeSubEl) modeSubEl.innerText = `Unimodal (Count: ${stats.maxFreq || 2})`;
        } else if (stats.modes.length === 2) {
          modeEl.innerText = `${stats.modes[0].toFixed(1)}, ${stats.modes[1].toFixed(1)}`;
          if (modeSubEl) modeSubEl.innerText = `Bimodal (Count: ${stats.maxFreq || 2})`;
        } else {
          modeEl.innerText = `${stats.modes.slice(0, 2).map(m => m.toFixed(1)).join(', ')}...`;
          if (modeSubEl) modeSubEl.innerText = `Multimodal (${stats.modes.length} modes)`;
        }
      } else {
        modeEl.innerText = 'No Mode';
        if (modeSubEl) modeSubEl.innerText = 'All unique (freq = 1)';
      }
    }

    // Interquartile Range (IQR)
    const iqrEl = document.getElementById('descIQR');
    const iqrSubEl = document.getElementById('descIQRSub');
    if (iqrEl) iqrEl.innerText = stats.iqr.toFixed(2);
    if (iqrSubEl) iqrSubEl.innerText = `Q1: ${stats.q1.toFixed(2)} | Q3: ${stats.q3.toFixed(2)}`;

    document.getElementById('descSD').innerText = stats.sd.toFixed(2);
    const varEl = document.getElementById('descVar');
    if (varEl) varEl.innerText = `Var: ${stats.variance.toFixed(2)}`;

    document.getElementById('descCI95').innerText = `[${stats.ci95[0].toFixed(2)}, ${stats.ci95[1].toFixed(2)}]`;
    const semEl = document.getElementById('descSEM');
    if (semEl) semEl.innerText = `SEM: ±${stats.sem.toFixed(2)}`;

    // Skewness
    const skewValEl = document.getElementById('descSkewnessVal');
    const skewSubEl = document.getElementById('descSkewnessSub');
    if (skewValEl) skewValEl.innerText = stats.skewness.toFixed(2);
    if (skewSubEl) skewSubEl.innerText = stats.skewnessInterpretation;

    // Kurtosis
    const kurtValEl = document.getElementById('descKurtosisVal');
    const kurtSubEl = document.getElementById('descKurtosisSub');
    if (kurtValEl) kurtValEl.innerText = (stats.kurtosis > 0 ? '+' : '') + stats.kurtosis.toFixed(2);
    if (kurtSubEl) kurtSubEl.innerText = stats.kurtosisInterpretation;

    // Outliers (Tukey's fences)
    const outValEl = document.getElementById('descOutliersVal');
    const outSubEl = document.getElementById('descOutliersSub');
    if (outValEl) {
      outValEl.innerText = `${stats.outliers.length} Detected`;
      if (stats.outliers.length > 0) {
        const extremeCount = stats.outliers.filter(o => o.type === 'Extreme').length;
        const mildCount = stats.outliers.length - extremeCount;
        const parts = [];
        if (mildCount > 0) parts.push(`${mildCount} Mild`);
        if (extremeCount > 0) parts.push(`${extremeCount} Extreme`);
        if (outSubEl) outSubEl.innerText = parts.join(', ');
      } else {
        if (outSubEl) outSubEl.innerText = `Fences [${stats.lowerFence.toFixed(1)}, ${stats.upperFence.toFixed(1)}]`;
      }
    }

    // Normality Assessment (Jarque-Bera Test)
    document.getElementById('descNormality').innerText = stats.normality.isNormal ? 'Normal ✓' : 'Non-Normal ⚠️';
    const normSubEl = document.getElementById('descNormalitySub');
    if (normSubEl) normSubEl.innerText = `JB = ${stats.normality.statistic.toFixed(2)} | p = ${stats.normality.pValue.toFixed(3)}`;
    const legacySkew = document.getElementById('descSkewness');
    if (legacySkew) legacySkew.innerText = `Skew: ${stats.skewness.toFixed(2)}`;

    // Outlier & Distribution Diagnostics Banner
    const diagBanner = document.getElementById('descDiagBanner');
    const diagBadge = document.getElementById('descDiagBadge');
    const diagTitle = document.getElementById('descDiagTitle');
    const diagDetail = document.getElementById('descDiagDetail');

    if (diagBanner && diagBadge && diagTitle && diagDetail) {
      const hasOutliers = stats.outliers && stats.outliers.length > 0;
      const isSkewed = Math.abs(stats.skewness) > 0.5 || !stats.normality.isNormal;
      const hasExtreme = stats.outliers && stats.outliers.some(o => o.type === 'Extreme');

      if (hasOutliers || isSkewed) {
        diagBanner.style.display = 'block';
        if (hasExtreme || Math.abs(stats.skewness) > 1.0) {
          diagBadge.className = 'badge badge-danger';
          diagBadge.innerText = 'High Anomaly / Skew';
        } else {
          diagBadge.className = 'badge badge-warning';
          diagBadge.innerText = 'Moderate Anomaly';
        }

        const issues = [];
        if (hasOutliers) {
          const outList = stats.outliers.map(o => `${o.value.toFixed(1)} (${o.type} ${o.direction}, Z = ${o.zScore > 0 ? '+' : ''}${o.zScore.toFixed(2)})`).join(', ');
          issues.push(`<strong>${stats.outliers.length} outlier(s)</strong> detected outside Tukey fences [${stats.lowerFence.toFixed(1)}, ${stats.upperFence.toFixed(1)}]: ${outList}.`);
        }
        if (isSkewed) {
          issues.push(`Distribution displays <strong>${stats.skewnessInterpretation}</strong> (G₁ = ${stats.skewness.toFixed(2)}) and <strong>${stats.kurtosisInterpretation}</strong> (excess G₂ = ${stats.kurtosis.toFixed(2)}).`);
        }

        diagTitle.innerText = hasOutliers ? 'Distribution Anomaly: Outliers & Skewness Detected' : 'Distribution Notice: Asymmetric Skewness';
        diagDetail.innerHTML = issues.join(' ') +
          `<br><span style="display:inline-block; margin-top:0.35rem; color: var(--cyan-primary);">💡 <strong>Clinical Guidance:</strong> Reporting <strong>Median (${stats.median.toFixed(2)}) &amp; IQR (${stats.iqr.toFixed(2)})</strong> is strongly recommended over Mean ± SD. For inferential comparison, use non-parametric tests (Mann-Whitney U or Kruskal-Wallis) to prevent outlier distortion.</span>`;
      } else {
        diagBanner.style.display = 'block';
        diagBadge.className = 'badge badge-sig';
        diagBadge.innerText = 'Normal Distribution';
        diagTitle.innerText = 'Standard Gaussian Distribution Verified';
        diagDetail.innerHTML = `No outliers detected within Tukey's fences [${stats.lowerFence.toFixed(1)}, ${stats.upperFence.toFixed(1)}]. Skewness (${stats.skewness.toFixed(2)}) and Kurtosis (${stats.kurtosis.toFixed(2)}) conform to normal standards (Jarque-Bera p = ${stats.normality.pValue.toFixed(3)}). <strong>Reporting Mean ± SD is statistically sound.</strong>`;
      }
    }

    // Descriptive Report Statement
    const repEl = document.getElementById('descReportText');
    if (repEl) {
      repEl.innerText = Exporter.formatDescriptiveReport(stats);
    }

    // Render 1: Histogram & Fitted Normal Curve
    if (this.engines.descCanvas) {
      Plots.renderHistogram(this.engines.descCanvas, stats, 'Frequency Distribution & Normal Fit');
    }

    // Render 2: Box & Whiskers Plot
    if (this.engines.descBoxCanvas) {
      Plots.renderBoxPlot(this.engines.descBoxCanvas, [{ name: 'Sample Data', stats, color: this.engines.descBoxCanvas.palette.primary }], 'Box & Whiskers (Tukey Fences & Outliers)');
    }

    // Render 3: Violin Density Plot (KDE)
    if (this.engines.descViolinCanvas) {
      Plots.renderViolinPlot(this.engines.descViolinCanvas, [{ name: 'Sample Data', stats, color: this.engines.descViolinCanvas.palette.secondary }], 'Violin Density Plot (KDE & Quartiles)');
    }
    this.results = this.results || {}; this.results.descriptive = stats;
  }

  runHypothesis() {
    const rawA = document.getElementById('hypoGroupA')?.value || '';
    const rawB = document.getElementById('hypoGroupB')?.value || '';
    const nameA = document.getElementById('hypoNameA')?.value || 'Group A';
    const nameB = document.getElementById('hypoNameB')?.value || 'Group B';
    const testType = document.getElementById('hypoTestType')?.value || 'welch';

    const dataA = DataParser.parseSeries(rawA);
    const dataB = DataParser.parseSeries(rawB);

    let res;
    if (testType === 'student') {
      res = Hypothesis.independentTTest(dataA, dataB);
    } else if (testType === 'welch') {
      res = Hypothesis.welchTTest(dataA, dataB);
    } else if (testType === 'paired') {
      res = Hypothesis.pairedTTest(dataA, dataB);
    } else {
      res = Hypothesis.mannWhitneyUTest(dataA, dataB);
    }

    if (res.error) {
      alert(res.error);
      return;
    }

    res.groupA = Object.assign(res.groupA || Descriptive.calculate(dataA), { name: nameA });
    res.groupB = Object.assign(res.groupB || Descriptive.calculate(dataB), { name: nameB });

    document.getElementById('hypoStat').innerText = (res.statistic || res.zScore || 0).toFixed(2);
    document.getElementById('hypoPVal').innerText = Exporter.formatP(res.pValue);
    document.getElementById('hypoPValBadge').className = `badge ${res.isSignificant ? 'badge-sig' : 'badge-ns'}`;
    document.getElementById('hypoPValBadge').innerText = res.isSignificant ? 'Significant (p < .05)' : 'Not Significant (ns)';
    document.getElementById('hypoEffect').innerText = (res.cohensD !== undefined ? res.cohensD.toFixed(2) : (res.rankBiserial || 0).toFixed(2));
    document.getElementById('hypoDiff').innerText = res.meanDiff !== undefined ? res.meanDiff.toFixed(2) : 'N/A';

    // Clinical Report Box
    const reportText = Exporter.formatTTestReport(res);
    document.getElementById('hypoReportText').innerText = reportText;

    // Render Dispersion Plot (95% CI, SEM, SD, or IQR Box & Whiskers)
    const errorBarMode = document.getElementById('hypoErrorBarMode')?.value || 'ci95';
    const modeDescriptions = {
      ci95: 'Error Bars: 95% Confidence Interval (Mean ± 95% CI)',
      sem: 'Error Bars: Standard Error of Mean (Mean ± 1 SEM)',
      sd: 'Error Bars: Standard Deviation (Mean ± 1 SD)',
      iqr: 'Distribution: Box & Whiskers (Median, Q1-Q3 IQR, Tukey Fences)'
    };
    const subElem = document.getElementById('hypoChartSub');
    if (subElem) {
      subElem.innerText = modeDescriptions[errorBarMode] || modeDescriptions.ci95;
    }

    if (this.engines.hypoCanvas) {
      Plots.renderErrorBarPlot(this.engines.hypoCanvas, [
        { name: nameA, stats: res.groupA },
        { name: nameB, stats: res.groupB }
      ], {
        mode: errorBarMode,
        title: `${nameA} vs ${nameB}`
      });
    }
    this.results = this.results || {}; this.results.hypothesis = res;
  }

  runAnova() {
    const g1 = DataParser.parseSeries(document.getElementById('anovaG1')?.value || '');
    const g2 = DataParser.parseSeries(document.getElementById('anovaG2')?.value || '');
    const g3 = DataParser.parseSeries(document.getElementById('anovaG3')?.value || '');

    const groups = [
      { name: document.getElementById('anovaName1')?.value || 'Cohort 1', data: g1 },
      { name: document.getElementById('anovaName2')?.value || 'Cohort 2', data: g2 },
      { name: document.getElementById('anovaName3')?.value || 'Cohort 3', data: g3 }
    ];

    const res = Anova.oneWay(groups);
    if (res.error) {
      alert(res.error);
      return;
    }

    document.getElementById('anovaF').innerText = res.fStatistic.toFixed(2);
    document.getElementById('anovaP').innerText = Exporter.formatP(res.pValue);
    document.getElementById('anovaEta').innerText = res.etaSquared.toFixed(3);
    document.getElementById('anovaOmega').innerText = res.omegaSquared.toFixed(3);

    // Report
    document.getElementById('anovaReportText').innerText = Exporter.formatAnovaReport(res);

    // Render Dispersion Plot (95% CI, SEM, SD, or IQR Box & Whiskers)
    const errorBarMode = document.getElementById('anovaErrorBarMode')?.value || 'ci95';
    const modeDescriptions = {
      ci95: 'Error Bars: 95% Confidence Interval (Mean ± 95% CI)',
      sem: 'Error Bars: Standard Error of Mean (Mean ± 1 SEM)',
      sd: 'Error Bars: Standard Deviation (Mean ± 1 SD)',
      iqr: 'Distribution: Box & Whiskers (Median, Q1-Q3 IQR, Tukey Fences)'
    };
    const subElem = document.getElementById('anovaChartSub');
    if (subElem) {
      subElem.innerText = modeDescriptions[errorBarMode] || modeDescriptions.ci95;
    }

    if (this.engines.anovaCanvas) {
      Plots.renderErrorBarPlot(this.engines.anovaCanvas, res.groups, {
        mode: errorBarMode,
        title: 'Multi-Cohort Comparison'
      });
    }
    this.results = this.results || {}; this.results.anova = res;
  }

  runCategorical() {
    const mode = document.getElementById('catAnalysisMode')?.value || 'diagnostic';
    const a = parseFloat(document.getElementById('catA')?.value) || 0;
    const b = parseFloat(document.getElementById('catB')?.value) || 0;
    const c = parseFloat(document.getElementById('catC')?.value) || 0;
    const d = parseFloat(document.getElementById('catD')?.value) || 0;

    const cornerHeader = document.getElementById('catCornerHeader');
    const col1Header = document.getElementById('catCol1Header');
    const col2Header = document.getElementById('catCol2Header');
    const row1Header = document.getElementById('catRow1Header');
    const row2Header = document.getElementById('catRow2Header');
    const labelA = document.getElementById('catLabelA');
    const labelB = document.getElementById('catLabelB');
    const labelC = document.getElementById('catLabelC');
    const labelD = document.getElementById('catLabelD');
    const explanation = document.getElementById('catExplanation');
    const diagGrid = document.getElementById('catDiagnosticMetrics');
    const studyGrid = document.getElementById('catStudyMetrics');
    const reportTitle = document.getElementById('catReportHeaderTitle');

    if (mode === 'diagnostic') {
      if (cornerHeader) cornerHeader.innerText = 'Test \\ Ref';
      if (col1Header) col1Header.innerText = 'Gold Standard (+)';
      if (col2Header) col2Header.innerText = 'Gold Standard (-)';
      if (row1Header) row1Header.innerText = 'New Test (+)';
      if (row2Header) row2Header.innerText = 'New Test (-)';
      if (labelA) labelA.innerText = 'True Positive (TP)';
      if (labelB) labelB.innerText = 'False Positive (FP)';
      if (labelC) labelC.innerText = 'False Negative (FN)';
      if (labelD) labelD.innerText = 'True Negative (TN)';
      if (explanation) explanation.innerText = 'Evaluates index test accuracy against the reference standard: Sensitivity, Specificity, PPV, NPV, Overall Accuracy (Wilson Score 95% CIs), and Likelihood Ratios.';
      if (diagGrid) diagGrid.style.display = 'grid';
      if (studyGrid) studyGrid.style.display = 'none';
      if (reportTitle) reportTitle.innerText = 'Diagnostic Performance & Accuracy Report (STARD compliant)';

      const res = Diagnostic.evaluate2x2(a, b, c, d);
      if (res.error) {
        document.getElementById('catReportText').innerText = res.error;
        return;
      }

      const formatPct = (val) => `${(val * 100).toFixed(1)}%`;
      const formatCI = (ci) => `95% CI: [${(ci[0] * 100).toFixed(1)}%, ${(ci[1] * 100).toFixed(1)}%]`;

      const sensEl = document.getElementById('diagSens');
      if (sensEl) sensEl.innerText = formatPct(res.sensitivity);
      const sensCIEl = document.getElementById('diagSensCI');
      if (sensCIEl) sensCIEl.innerText = formatCI(res.sensitivityCI95);

      const specEl = document.getElementById('diagSpec');
      if (specEl) specEl.innerText = formatPct(res.specificity);
      const specCIEl = document.getElementById('diagSpecCI');
      if (specCIEl) specCIEl.innerText = formatCI(res.specificityCI95);

      const ppvEl = document.getElementById('diagPPV');
      if (ppvEl) ppvEl.innerText = formatPct(res.ppv);
      const ppvCIEl = document.getElementById('diagPPVCI');
      if (ppvCIEl) ppvCIEl.innerText = formatCI(res.ppvCI95);

      const npvEl = document.getElementById('diagNPV');
      if (npvEl) npvEl.innerText = formatPct(res.npv);
      const npvCIEl = document.getElementById('diagNPVCI');
      if (npvCIEl) npvCIEl.innerText = formatCI(res.npvCI95);

      const accEl = document.getElementById('diagAcc');
      if (accEl) accEl.innerText = formatPct(res.accuracy);
      const accCIEl = document.getElementById('diagAccCI');
      if (accCIEl) accCIEl.innerText = formatCI(res.accuracyCI95);

      const lrEl = document.getElementById('diagLR');
      if (lrEl) lrEl.innerText = `LR+ ${isFinite(res.plr) ? res.plr.toFixed(2) : '∞'} | LR- ${isFinite(res.nlr) ? res.nlr.toFixed(2) : '0'}`;

      const youdenEl = document.getElementById('diagYouden');
      if (youdenEl) youdenEl.innerText = `Youden's J: ${res.youdenJ.toFixed(3)} | Prev: ${formatPct(res.prevalence)}`;

      const report = `Diagnostic test evaluation against gold standard reference (Total N = ${res.total}): Overall Accuracy = ${formatPct(res.accuracy)} (${formatCI(res.accuracyCI95)}). ` +
        `Sensitivity (TPR) = ${formatPct(res.sensitivity)} (${formatCI(res.sensitivityCI95)}), ` +
        `Specificity (TNR) = ${formatPct(res.specificity)} (${formatCI(res.specificityCI95)}). ` +
        `Positive Predictive Value (PPV) = ${formatPct(res.ppv)} (${formatCI(res.ppvCI95)}), ` +
        `Negative Predictive Value (NPV) = ${formatPct(res.npv)} (${formatCI(res.npvCI95)}). ` +
        `Positive Likelihood Ratio (LR+) = ${isFinite(res.plr) ? res.plr.toFixed(2) : 'N/A'}, ` +
        `Negative Likelihood Ratio (LR-) = ${isFinite(res.nlr) ? res.nlr.toFixed(2) : 'N/A'}, ` +
        `Youden's J index = ${res.youdenJ.toFixed(3)}, Sample Disease Prevalence = ${formatPct(res.prevalence)}.`;

      document.getElementById('catReportText').innerText = report;
    } else {
      if (cornerHeader) cornerHeader.innerText = 'Cohort \\ Event';
      if (col1Header) col1Header.innerText = 'Intervention';
      if (col2Header) col2Header.innerText = 'Control';
      if (row1Header) row1Header.innerText = 'Event (+)';
      if (row2Header) row2Header.innerText = 'No Event (-)';
      if (labelA) labelA.innerText = 'Treated Event';
      if (labelB) labelB.innerText = 'Control Event';
      if (labelC) labelC.innerText = 'Treated No Event';
      if (labelD) labelD.innerText = 'Control No Event';
      if (explanation) explanation.innerText = 'Instant computation of Odds Ratio (Woolf 95% CI), Relative Risk, Absolute Risk Reduction, Number Needed to Treat (NNT), Pearson Chi-Square, and Fisher\'s exact test.';
      if (diagGrid) diagGrid.style.display = 'none';
      if (studyGrid) studyGrid.style.display = 'grid';
      if (reportTitle) reportTitle.innerText = 'Epidemiological & Clinical Risk Report';

      const res = Categorical.twoByTwo(a, b, c, d);
      if (res.error) {
        document.getElementById('catReportText').innerText = res.error;
        return;
      }

      const r = res.riskMetrics;
      document.getElementById('catOR').innerText = `${r.oddsRatio.toFixed(2)} [${r.orCI95[0].toFixed(2)}, ${r.orCI95[1].toFixed(2)}]`;
      document.getElementById('catRR').innerText = `${r.relativeRisk.toFixed(2)} [${r.rrCI95[0].toFixed(2)}, ${r.rrCI95[1].toFixed(2)}]`;
      document.getElementById('catChiSq').innerText = res.chiSquare.standard.toFixed(2);
      document.getElementById('catPVal').innerText = Exporter.formatP(res.chiSquare.pValueStandard);
      document.getElementById('catFisher').innerText = Exporter.formatP(res.fishersExact.pValue);
      document.getElementById('catNNT').innerText = isFinite(r.nnt) ? r.nnt.toFixed(1) : '∞';

      const arrSubEl = document.getElementById('catARRSub');
      if (arrSubEl) arrSubEl.innerText = `ARR: ${(r.arr * 100).toFixed(1)}%`;
      const rrrEl = document.getElementById('catRRR');
      if (rrrEl) rrrEl.innerText = isFinite(r.rrr) ? `${(r.rrr * 100).toFixed(1)}%` : 'N/A';

      document.getElementById('catReportText').innerText = Exporter.formatCategoricalReport(res);
    }
    this.results = this.results || {}; this.results.categorical = res;
  }

  runCorrelation() {
    const xRaw = DataParser.parseSeries(document.getElementById('corrX')?.value || '');
    const yRaw = DataParser.parseSeries(document.getElementById('corrY')?.value || '');

    const corr = Correlation.pearson(xRaw, yRaw);
    const spearman = Correlation.spearman(xRaw, yRaw);
    const reg = Correlation.linearRegression(xRaw, yRaw);

    if (corr.error || reg.error) {
      alert(corr.error || reg.error);
      return;
    }

    document.getElementById('corrR').innerText = corr.r.toFixed(3);
    document.getElementById('corrR2').innerText = corr.rSquared.toFixed(3);
    document.getElementById('corrRho').innerText = spearman.rho.toFixed(3);
    document.getElementById('corrP').innerText = Exporter.formatP(corr.pValue);
    document.getElementById('regEq').innerText = reg.equation;

    // Check Assumptions
    const assump = Correlation.checkAssumptions(xRaw, yRaw, reg);
    if (!assump.error) {
      const ax = document.getElementById('corrAssumpX');
      const ay = document.getElementById('corrAssumpY');
      const ao = document.getElementById('corrAssumpOutlier');
      const ah = document.getElementById('corrAssumpHomo');
      const as = document.getElementById('corrAssumpShape');

      if (ax) ax.innerHTML = assump.xNormal ? '<span style="color: var(--emerald-primary);">Normal ✓</span>' : '<span style="color: var(--amber-primary);">Skewed ⚠️</span>';
      if (ay) ay.innerHTML = assump.yNormal ? '<span style="color: var(--emerald-primary);">Normal ✓</span>' : '<span style="color: var(--amber-primary);">Skewed ⚠️</span>';
      if (ao) ao.innerHTML = assump.outliers.length === 0 ? '<span style="color: var(--emerald-primary);">None ✓</span>' : `<span style="color: var(--rose-primary);">${assump.outliers.length} Outlier(s) ⚠️</span>`;
      if (ah) ah.innerHTML = assump.isHomoscedastic ? '<span style="color: var(--emerald-primary);">Equal Var ✓</span>' : '<span style="color: var(--amber-primary);">Heteroscedastic ⚠️</span>';
      if (as) as.innerHTML = assump.isUShaped ? `<span style="color: var(--amber-primary); font-weight: 700;">${assump.shape} ⚠️</span>` : '<span style="color: var(--emerald-primary);">Linear ✓</span>';

      const banner = document.getElementById('corrRecBanner');
      const badge = document.getElementById('corrRecBadge');
      const title = document.getElementById('corrRecTitle');
      const detail = document.getElementById('corrRecDetail');

      if (badge) {
        badge.className = `badge ${assump.isParametricOk ? 'badge-sig' : 'badge-ns'}`;
        badge.innerText = assump.isParametricOk ? 'Parametric Valid' : (assump.isUShaped ? 'Non-Linear Indicated' : 'Non-Parametric Indicated');
      }
      if (title) title.innerText = assump.recommendedTest;
      if (detail) detail.innerText = assump.recommendationDetail;
      if (banner) {
        banner.style.background = assump.isParametricOk ? 'rgba(0, 210, 255, 0.08)' : 'rgba(245, 158, 11, 0.08)';
        banner.style.borderColor = assump.isParametricOk ? 'var(--cyan-primary)' : 'var(--amber-primary)';
      }

      let report = `Bivariate Statistical Evaluation (N = ${xRaw.length} pairs):\n`;
      report += `• Assumption Diagnostics: Bivariate normality is ${assump.bivariateNormal ? 'satisfied' : 'violated'}; ` +
        `curve morphology is ${assump.isUShaped ? `strongly non-linear (${assump.shape}, nadir/zenith at X = ${assump.quad.vertexX.toFixed(2)})` : (assump.isLinear ? 'linear' : 'monotonic non-linear')}; ` +
        `homoscedasticity is ${assump.isHomoscedastic ? 'preserved' : 'violated'}; ` +
        `${assump.outliers.length === 0 ? 'no extreme leverage outliers detected' : `${assump.outliers.length} influential outlier(s) detected`}.\n`;
      report += `• Method Decision: ${assump.recommendedTest} is the mathematically recommended procedure.\n`;
      if (assump.isUShaped) {
        report += `• Non-Linear Quadratic Model: y = ${assump.quad.equation} (R² = ${assump.quad.rSquaredQuad.toFixed(3)}, Incremental F = ${assump.quad.fStat.toFixed(2)}, ${Exporter.formatP(assump.quad.pQuad)}).\n` +
          `• Clinical Morphology: ${assump.shape} relationship with physiological ${assump.quad.b2 > 0 ? 'optimum / nadir' : 'peak / zenith'} at X = ${assump.quad.vertexX.toFixed(2)} (predicted Y = ${assump.quad.vertexY.toFixed(2)}). Standard linear Pearson (r = ${corr.r.toFixed(3)}) fails to capture this strong biological curve.`;
      } else if (assump.isParametricOk) {
        report += `• Results: Pearson's correlation r = ${corr.r.toFixed(3)} (${corr.isSignificant ? 'statistically significant' : 'ns'}, ${Exporter.formatP(corr.pValue)}), R² = ${corr.rSquared.toFixed(3)}. Linear regression equation: ${reg.equation} (residual SE = ${reg.seResidual.toFixed(3)}).`;
      } else {
        report += `• Results: Spearman's rank correlation ρ = ${spearman.rho.toFixed(3)} (${spearman.pValue < 0.05 ? 'statistically significant' : 'ns'}, ${Exporter.formatP(spearman.pValue)}). (Reference Pearson r = ${corr.r.toFixed(3)}, R² = ${corr.rSquared.toFixed(3)}).`;
      }
      const rEl = document.getElementById('corrReportText');
      if (rEl) rEl.innerText = report;
    }

    // Render Scatter + Regression / Quadratic
    const pairs = Correlation.cleanPairs(xRaw, yRaw);
    if (this.engines.corrCanvas) {
      const plotTitle = assump && assump.isUShaped ? `Curvilinear Regression (${assump.shape})` : 'Scatter Plot & Linear Regression';
      Plots.renderScatterRegression(this.engines.corrCanvas, pairs, reg, plotTitle, assump?.quad);
    }
    this.results = this.results || {}; this.results.correlation = { n: xs.length, correlation: p, regression: reg, morphology: assump };
  }

  runDiagnostic() {
    const text = document.getElementById('rocInput')?.value || '';
    const lines = text.trim().split(/\r?\n/);
    const samples = [];
    for (const line of lines) {
      const parts = line.split(/[,\t\s]+/).filter(Boolean);
      if (parts.length >= 2) {
        samples.push({ score: parseFloat(parts[0]), status: parseInt(parts[1], 10) });
      }
    }

    const roc = Diagnostic.computeROC(samples);
    if (roc.error) {
      alert(roc.error);
      return;
    }

    // 2x2 at optimal cutoff
    const best = roc.optimalCutoff;
    const tp = Math.round(best.sens * roc.nPos);
    const fp = Math.round((1 - best.spec) * roc.nNeg);
    const fn = roc.nPos - tp;
    const tn = roc.nNeg - fp;
    const eval2x2 = Diagnostic.evaluate2x2(tp, fp, fn, tn);

    document.getElementById('rocAUC').innerText = `${roc.auc.toFixed(3)} [${roc.aucCI95[0].toFixed(3)}, ${roc.aucCI95[1].toFixed(3)}]`;
    document.getElementById('rocCutoff').innerText = typeof best.threshold === 'number' ? best.threshold.toFixed(2) : best.threshold;
    document.getElementById('rocSens').innerText = `${(eval2x2.sensitivity * 100).toFixed(1)}%`;
    document.getElementById('rocSpec').innerText = `${(eval2x2.specificity * 100).toFixed(1)}%`;
    document.getElementById('rocPPV').innerText = `${(eval2x2.ppv * 100).toFixed(1)}%`;
    document.getElementById('rocNPV').innerText = `${(eval2x2.npv * 100).toFixed(1)}%`;
    document.getElementById('rocPLR').innerText = eval2x2.plr.toFixed(2);
    document.getElementById('rocNLR').innerText = eval2x2.nlr.toFixed(2);

    document.getElementById('rocReportText').innerText = Exporter.formatDiagnosticReport(eval2x2, roc);

    if (this.engines.rocCanvas) {
      Plots.renderROC(this.engines.rocCanvas, roc, `ROC Curve (AUC = ${roc.auc.toFixed(3)})`);
    }
    this.results = this.results || {}; this.results.diagnostic = { eval2x2: d, roc, best };
  }

  runPower() {
    const design = document.getElementById('pwrStudyDesign')?.value || 'independent';
    const goal = document.getElementById('pwrCalcGoal')?.value || 'sample_size';
    const alpha = parseFloat(document.getElementById('pwrAlpha')?.value) || 0.05;
    const power = parseFloat(document.getElementById('pwrPower')?.value) || 0.80;
    const givenN = goal === 'power' ? (parseFloat(document.getElementById('pwrGivenN')?.value) || 50) : null;

    const effectLabel = document.getElementById('pwrEffectLabel');
    const effectVal = document.getElementById('pwrD');
    const effectSub = document.getElementById('pwrEffectSub');
    const nGroupLabel = document.getElementById('pwrNGroupLabel');
    const nGroupVal = document.getElementById('pwrNGroup');
    const nGroupSub = document.getElementById('pwrNGroupSub');
    const nTotalLabel = document.getElementById('pwrNTotalLabel');
    const nTotalVal = document.getElementById('pwrNTotal');
    const nTotalSub = document.getElementById('pwrNTotalSub');
    const achievedPowerVal = document.getElementById('pwrAchievedPower');
    const powerSub = document.getElementById('pwrPowerSub');
    const cardFisher = document.getElementById('pwrCardFisher');
    const cardRisk = document.getElementById('pwrCardClinicalRisk');
    const fisherNVal = document.getElementById('pwrFisherN');
    const fisherSub = document.getElementById('pwrFisherSub');
    const riskVal = document.getElementById('pwrClinicalRisk');
    const riskSub = document.getElementById('pwrClinicalRiskSub');
    const reportBox = document.getElementById('pwrReportText');

    if (design === 'independent') {
      if (cardFisher) cardFisher.style.display = 'none';
      if (cardRisk) cardRisk.style.display = 'none';

      const m1 = parseFloat(document.getElementById('pwrM1')?.value) || 10;
      const m2 = parseFloat(document.getElementById('pwrM2')?.value) || 15;
      const sd = parseFloat(document.getElementById('pwrSD')?.value) || 10;

      const res = PowerAnalysis.twoIndependentMeans({ m1, m2, sd, alpha, power, nPerGroup: givenN });
      if (res.error) {
        if (reportBox) reportBox.innerText = res.error;
        return;
      }

      if (effectLabel) effectLabel.innerText = "Effect Size (Cohen's d)";
      if (effectVal) effectVal.innerText = res.cohensD.toFixed(2);
      if (effectSub) effectSub.innerText = `|Δ| = ${res.diff.toFixed(2)} | SD = ${res.sd.toFixed(2)}`;

      if (nGroupLabel) nGroupLabel.innerText = 'Sample Size per Group (n)';
      if (nGroupVal) nGroupVal.innerText = res.nPerGroup;
      if (nGroupSub) nGroupSub.innerText = 'Equal 1:1 allocation';

      if (nTotalLabel) nTotalLabel.innerText = 'Total Study Sample Size (N)';
      if (nTotalVal) nTotalVal.innerText = res.totalN;
      if (nTotalSub) nTotalSub.innerText = 'Across both arms (n₁ + n₂)';

      if (achievedPowerVal) achievedPowerVal.innerText = `${(res.achievedPower * 100).toFixed(1)}%`;
      if (powerSub) powerSub.innerText = `Type II error β = ${(1 - res.achievedPower).toFixed(3)} | α = ${res.alpha}`;

      let report = '';
      if (goal === 'sample_size') {
        report = `For a two-arm parallel randomized study (independent samples t-test) to detect a standardized effect size of Cohen's d = ${res.cohensD.toFixed(2)} (mean difference |Δ| = ${res.diff.toFixed(2)}, pooled SD = ${res.sd.toFixed(2)}) with ${Math.round(res.power * 100)}% statistical power at two-sided α = ${res.alpha}, a minimum of n = ${res.nPerGroup} patients per group (total N = ${res.totalN}) is required. Allowing for a standard 10% loss-to-follow-up buffer, a target enrollment of N = ${Math.ceil(res.totalN / 0.9)} patients (${Math.ceil(res.nPerGroup / 0.9)} per arm) is recommended.`;
      } else {
        report = `With an enrolled sample size of n = ${res.nPerGroup} patients per group (total N = ${res.totalN}) comparing two independent continuous means (Cohen's d = ${res.cohensD.toFixed(2)}, |Δ| = ${res.diff.toFixed(2)}, SD = ${res.sd.toFixed(2)}), the achieved statistical power to detect a true difference at two-sided α = ${res.alpha} is ${(res.achievedPower * 100).toFixed(1)}% (Type II error rate β = ${(1 - res.achievedPower).toFixed(3)}).`;
      }
      if (reportBox) reportBox.innerText = report;

    } else if (design === 'paired') {
      if (cardFisher) cardFisher.style.display = 'none';
      if (cardRisk) cardRisk.style.display = 'none';

      const m1 = parseFloat(document.getElementById('pwrPairedM1')?.value) || 120;
      const m2 = parseFloat(document.getElementById('pwrPairedM2')?.value) || 112;
      const sdDiff = parseFloat(document.getElementById('pwrPairedSD')?.value) || 10;

      const res = PowerAnalysis.pairedMeans({ m1, m2, sdDiff, alpha, power, nPairs: givenN });
      if (res.error) {
        if (reportBox) reportBox.innerText = res.error;
        return;
      }

      if (effectLabel) effectLabel.innerText = "Effect Size (Cohen's d_z)";
      if (effectVal) effectVal.innerText = res.dz.toFixed(2);
      if (effectSub) effectSub.innerText = `|Δ| = ${res.diff.toFixed(2)} | σ_d = ${res.sdDiff.toFixed(2)}`;

      if (nGroupLabel) nGroupLabel.innerText = 'Required Number of Pairs (N)';
      if (nGroupVal) nGroupVal.innerText = res.nPairs;
      if (nGroupSub) nGroupSub.innerText = 'Matched pairs / repeat measures';

      if (nTotalLabel) nTotalLabel.innerText = 'Total Paired Subjects (N)';
      if (nTotalVal) nTotalVal.innerText = res.totalN;
      if (nTotalSub) nTotalSub.innerText = 'Pre-to-post repeated subjects';

      if (achievedPowerVal) achievedPowerVal.innerText = `${(res.achievedPower * 100).toFixed(1)}%`;
      if (powerSub) powerSub.innerText = `Type II error β = ${(1 - res.achievedPower).toFixed(3)} | α = ${res.alpha}`;

      let report = '';
      if (goal === 'sample_size') {
        report = `In a paired / before-and-after repeated measures study design (paired t-test), evaluating a mean intra-individual shift of |Δ| = ${res.diff.toFixed(2)} (SD of paired differences σ_d = ${res.sdDiff.toFixed(2)}, standardized effect size Cohen's d_z = ${res.dz.toFixed(2)}) with ${Math.round(res.power * 100)}% statistical power at two-sided α = ${res.alpha}, a minimum of N = ${res.nPairs} paired subjects is required. Accounting for an estimated 10% withdrawal or uninterpretable paired follow-up, an enrollment of N = ${Math.ceil(res.nPairs / 0.9)} subjects is recommended.`;
      } else {
        report = `For a paired study design with N = ${res.nPairs} evaluable paired subjects evaluating a mean difference of |Δ| = ${res.diff.toFixed(2)} (SD of differences σ_d = ${res.sdDiff.toFixed(2)}, Cohen's d_z = ${res.dz.toFixed(2)}), the achieved statistical power at two-sided α = ${res.alpha} is ${(res.achievedPower * 100).toFixed(1)}% (Type II error rate β = ${(1 - res.achievedPower).toFixed(3)}).`;
      }
      if (reportBox) reportBox.innerText = report;

    } else if (design === 'contingency') {
      if (cardFisher) cardFisher.style.display = 'block';
      if (cardRisk) cardRisk.style.display = 'block';

      const p1 = parseFloat(document.getElementById('pwrP1')?.value) || 0.15;
      const p2 = parseFloat(document.getElementById('pwrP2')?.value) || 0.30;
      const propTestType = document.getElementById('pwrPropTestType')?.value || 'chisq';

      const res = PowerAnalysis.contingency2x2({ p1, p2, alpha, power, nPerGroup: givenN, testType: propTestType });
      if (res.error) {
        if (reportBox) reportBox.innerText = res.error;
        return;
      }

      if (effectLabel) effectLabel.innerText = 'Absolute Risk Reduction (ARR)';
      if (effectVal) effectVal.innerText = `${(res.arr * 100).toFixed(1)}%`;
      if (effectSub) effectSub.innerText = `Cohen's h = ${res.cohenH.toFixed(2)} | p₁=${(res.p1*100).toFixed(0)}% vs p₂=${(res.p2*100).toFixed(0)}%`;

      if (nGroupLabel) nGroupLabel.innerText = 'Sample Size per Group (n)';
      if (nGroupVal) nGroupVal.innerText = res.nPerGroup;
      if (nGroupSub) nGroupSub.innerText = `${res.testType === 'fisher' ? 'Continuity-corrected (Fleiss)' : 'Standard uncorrected'} 1:1`;

      if (nTotalLabel) nTotalLabel.innerText = 'Total Study Sample Size (N)';
      if (nTotalVal) nTotalVal.innerText = res.totalN;
      if (nTotalSub) nTotalSub.innerText = `Test + Control (2 × ${res.nPerGroup})`;

      if (achievedPowerVal) achievedPowerVal.innerText = `${(res.achievedPower * 100).toFixed(1)}%`;
      if (powerSub) powerSub.innerText = `Type II error β = ${(1 - res.achievedPower).toFixed(3)} | α = ${res.alpha}`;

      if (fisherNVal) fisherNVal.innerText = `n = ${res.continuityCorrectedN || res.nPerGroup}`;
      if (fisherSub) fisherSub.innerText = `Total N = ${2 * (res.continuityCorrectedN || res.nPerGroup)} (Fleiss / CC)`;

      if (riskVal) riskVal.innerText = `RR: ${res.rr.toFixed(2)} | OR: ${res.or.toFixed(2)}`;
      if (riskSub) riskSub.innerText = `NNT: ${res.nnt.toFixed(1)} | RRR: ${(res.rrr * 100).toFixed(1)}%`;

      let report = '';
      if (goal === 'sample_size') {
        report = `In a clinical trial comparing test (intervention) versus control cohorts with binary outcomes evaluated via a 2x2 contingency table (${res.testName}), expecting an event rate of ${(res.p1 * 100).toFixed(1)}% in the test arm versus ${(res.p2 * 100).toFixed(1)}% in the control arm (Absolute Risk Reduction ARR = ${(res.arr * 100).toFixed(1)}%, Relative Risk RR = ${res.rr.toFixed(2)}, Odds Ratio OR = ${res.or.toFixed(2)}, Number Needed to Treat NNT = ${res.nnt.toFixed(1)}):\n` +
          `• Uncorrected Pearson Chi-Square test requires n = ${res.uncorrectedN} patients per arm (total N = ${2 * res.uncorrectedN}) for ${Math.round(res.power * 100)}% power at two-sided α = ${res.alpha}.\n` +
          `• Fisher's Exact test / Continuity-Corrected Chi-Square (Fleiss & Casagrande-Pike-Smith formula) requires n = ${res.continuityCorrectedN} patients per arm (total N = ${2 * res.continuityCorrectedN}).\n` +
          `With an anticipated 15% loss-to-follow-up buffer, a target enrollment of N = ${Math.ceil((2 * res.nPerGroup) / 0.85)} patients (${Math.ceil(res.nPerGroup / 0.85)} per arm) is recommended.`;
      } else {
        report = `In a 2x2 contingency cohort study of n = ${res.nPerGroup} patients per group (total N = ${res.totalN}) comparing event rates of ${(res.p1 * 100).toFixed(1)}% (test) vs ${(res.p2 * 100).toFixed(1)}% (control) (ARR = ${(res.arr * 100).toFixed(1)}%, RR = ${res.rr.toFixed(2)}, OR = ${res.or.toFixed(2)}):\n` +
          `• Pearson Chi-Square test achieved power = ${(res.powerChisq * 100).toFixed(1)}% (β = ${(1 - res.powerChisq).toFixed(3)}).\n` +
          `• Fisher's Exact test / Continuity-Corrected achieved power = ${(res.powerFisher * 100).toFixed(1)}% (β = ${(1 - res.powerFisher).toFixed(3)}).`;
      }
      if (reportBox) reportBox.innerText = report;
    }
  }
}

// Instantiate on DOM load
window.addEventListener('DOMContentLoaded', () => {
  window.sgApp = new StatisGravityApp();
  window.app = window.sgApp;
});
  exportDocx(tabId) {
    if (!this.results || !this.results[tabId]) {
      switch (tabId) {
        case 'descriptive': this.runDescriptive(); break;
        case 'hypothesis': this.runHypothesis(); break;
        case 'anova': this.runAnova(); break;
        case 'categorical': this.runCategorical(); break;
        case 'correlation': this.runCorrelation(); break;
        case 'diagnostic': this.runROC(); break;
        case 'power': this.runPower(); break;
        case 'teaching': this.runTeaching(); break;
      }
    }
    const res = this.results ? this.results[tabId] : null;
    if (!res) {
      alert('Please run the analysis first before exporting the report.');
      return;
    }

    let exportData = {};
    if (tabId === 'descriptive') {
      exportData = {
        name: 'Continuous Variable',
        n: res.n,
        mean: res.mean,
        ci95: res.ci95,
        sd: res.sd,
        variance: res.variance,
        sem: res.sem,
        median: res.median,
        q1: res.q1,
        q3: res.q3,
        iqr: res.iqr,
        modes: res.modes,
        maxFreq: res.maxFreq,
        min: res.min,
        max: res.max,
        skewness: res.skewness,
        skewnessInterpretation: res.skewnessInterpretation,
        kurtosis: res.kurtosis,
        kurtosisInterpretation: res.kurtosisInterpretation,
        normality: res.normality,
        lowerFence: res.lowerFence,
        upperFence: res.upperFence,
        outliers: res.outliers,
        reportText: document.getElementById('descReportText')?.innerText,
        values: res.values
      };
    } else if (tabId === 'hypothesis') {
      exportData = {
        nameA: document.getElementById('hypoNameA')?.value || 'Cohort A',
        nameB: document.getElementById('hypoNameB')?.value || 'Cohort B',
        testName: res.testName,
        statistic: res.statistic,
        df: res.df,
        pValue: res.pValue,
        isSignificant: res.isSignificant,
        meanDiff: res.meanDiff,
        ci95: res.ci95,
        cohensD: res.cohensD,
        rankBiserial: res.rankBiserial,
        reportText: document.getElementById('hypoReportText')?.innerText,
        groupA: res.groupA,
        groupB: res.groupB
      };
    } else if (tabId === 'anova') {
      exportData = {
        k: res.k,
        totalN: res.totalN,
        grandMean: res.grandMean,
        ssBetween: res.ssBetween,
        dfBetween: res.dfBetween,
        msBetween: res.msBetween,
        ssWithin: res.ssWithin,
        dfWithin: res.dfWithin,
        msWithin: res.msWithin,
        ssTotal: res.ssTotal,
        fStatistic: res.fStatistic,
        pValue: res.pValue,
        etaSquared: res.etaSquared,
        omegaSquared: res.omegaSquared,
        reportText: document.getElementById('anovaReportText')?.innerText,
        groups: res.groups,
        pairwise: res.pairwise
      };
    } else if (tabId === 'categorical') {
      const mode = document.getElementById('catAnalysisMode')?.value || 'diagnostic';
      const a = parseFloat(document.getElementById('catA')?.value) || 0;
      const b = parseFloat(document.getElementById('catB')?.value) || 0;
      const c = parseFloat(document.getElementById('catC')?.value) || 0;
      const d = parseFloat(document.getElementById('catD')?.value) || 0;
      exportData = {
        mode,
        a, b, c, d,
        totalN: a + b + c + d,
        reportText: document.getElementById('catReportText')?.innerText,
        diag: res.diagnostic,
        chiSquare: res.chiSquare,
        fishersExact: res.fishersExact,
        risk: res.riskMetrics
      };
    } else if (tabId === 'correlation') {
      exportData = {
        xName: 'Independent Variable X',
        yName: 'Dependent Variable Y',
        n: res.n,
        reportText: document.getElementById('corrReportText')?.innerText,
        corr: res.correlation,
        reg: res.regression,
        morphology: res.morphology
      };
    } else if (tabId === 'diagnostic') {
      exportData = {
        name: 'Clinical Biomarker Evaluation',
        totalN: res.eval2x2.total,
        posCount: res.eval2x2.tp + res.eval2x2.fn,
        negCount: res.eval2x2.fp + res.eval2x2.tn,
        reportText: document.getElementById('rocReportText')?.innerText,
        auc: res.roc.auc,
        aucCI95: res.roc.aucCI95,
        seAuc: res.roc.seAuc,
        optimalCutoff: res.roc.optimalCutoff,
        sensitivity: res.eval2x2.sensitivity,
        specificity: res.eval2x2.specificity,
        plr: res.eval2x2.plr,
        nlr: res.eval2x2.nlr
      };
    } else if (tabId === 'power') {
      exportData = {
        designLabel: document.getElementById('pwrStudyDesign')?.selectedOptions?.[0]?.text || 'Study Design',
        goal: document.getElementById('pwrCalcGoal')?.value || 'sample_size',
        alpha: parseFloat(document.getElementById('pwrAlpha')?.value) || 0.05,
        effectSize: res.effectSize || 0.5,
        effectSizeLabel: res.effectSizeLabel || "Cohen's d",
        nPerGroup: res.nPerGroup || res.nPairs || 0,
        totalN: res.totalN || 0,
        targetPower: parseFloat(document.getElementById('pwrPower')?.value) || 0.80,
        achievedPower: res.achievedPower || 0.80,
        reportText: document.getElementById('pwrReportText')?.innerText,
        additionalMetrics: res.additionalMetrics || {}
      };
    } else if (tabId === 'teaching') {
      const dist = res.dist || {};
      const clt = res.clt || {};
      exportData = {
        distName: dist.name,
        distN: dist.n,
        clinicalExample: dist.clinicalExample,
        sampleMean: dist.mean,
        theoMean: dist.theoMean,
        sampleSD: dist.sd,
        theoSD: dist.theoSD,
        skewness: dist.skewness,
        skewnessLabel: dist.skewnessInterpretation,
        kurtosis: dist.kurtosis,
        kurtosisLabel: dist.kurtosisInterpretation,
        jbStat: dist.normality?.statistic,
        jbP: dist.normality?.pValue,
        isNormal: dist.normality?.isNormal,
        cltParentName: clt.population?.name,
        cltN: clt.sampleSize,
        cltK: clt.samplesDrawn,
        cltResults: {
          theoMean: clt.theoreticalMean,
          obsMean: clt.observedMean,
          theoSE: clt.theoreticalSE,
          obsSE: clt.observedSE,
          skewness: clt.skewness,
          normalityP: clt.normality?.pValue,
          isNormal: clt.normality?.isNormal
        },
        reportText: document.getElementById('teachingReportText')?.innerText
      };
    }

    Exporter.exportTabToDocx(tabId, exportData);
  }

  // --- Teaching & Simulation Module ---

  initTeachingModule() {
    this.renderTeachingParams();
    this.runTeachingDistribution();
    this.updateCltUI(Teaching.clt.getSummary());
  }

  renderTeachingParams() {
    const distKey = document.getElementById('teachingDistSelect')?.value || 'normal';
    const container = document.getElementById('teachingParamsContainer');
    if (!container) return;

    let html = '';
    switch (distKey) {
      case 'normal':
        html = `
          <div style="display: grid; grid-template-columns: 1fr 1fr; gap: 0.75rem;">
            <div>
              <label class="form-label" style="font-size: 0.78rem;">Mean (μ):</label>
              <input id="tp_norm_mean" type="number" step="0.5" value="12" class="form-control" style="padding: 0.35rem 0.5rem;">
            </div>
            <div>
              <label class="form-label" style="font-size: 0.78rem;">Std Dev (σ):</label>
              <input id="tp_norm_sd" type="number" min="0.1" step="0.5" value="3" class="form-control" style="padding: 0.35rem 0.5rem;">
            </div>
          </div>
        `;
        break;
      case 'studentsT':
        html = `
          <div style="display: grid; grid-template-columns: 1fr 1fr; gap: 0.75rem;">
            <div>
              <label class="form-label" style="font-size: 0.78rem;">Degrees of Freedom (ν):</label>
              <input id="tp_t_df" type="number" min="1" max="100" step="1" value="4" class="form-control" style="padding: 0.35rem 0.5rem;">
            </div>
            <div>
              <label class="form-label" style="font-size: 0.78rem;">Scale (s):</label>
              <input id="tp_t_scale" type="number" min="0.1" step="0.5" value="2" class="form-control" style="padding: 0.35rem 0.5rem;">
            </div>
          </div>
        `;
        break;
      case 'logNormal':
        html = `
          <div style="display: grid; grid-template-columns: 1fr 1fr; gap: 0.75rem;">
            <div>
              <label class="form-label" style="font-size: 0.78rem;">Log-Mean (μ):</label>
              <input id="tp_log_mu" type="number" step="0.1" value="1.5" class="form-control" style="padding: 0.35rem 0.5rem;">
            </div>
            <div>
              <label class="form-label" style="font-size: 0.78rem;">Log-SD (σ):</label>
              <input id="tp_log_sigma" type="number" min="0.1" step="0.1" value="0.6" class="form-control" style="padding: 0.35rem 0.5rem;">
            </div>
          </div>
        `;
        break;
      case 'exponential':
        html = `
          <div>
            <label class="form-label" style="font-size: 0.78rem;">Rate Parameter (λ):</label>
            <input id="tp_exp_rate" type="number" min="0.05" step="0.1" value="0.4" class="form-control" style="padding: 0.35rem 0.5rem;">
            <div style="font-size: 0.7rem; color: var(--text-muted); margin-top: 0.2rem;">Theoretical Mean = 1/λ = 2.50</div>
          </div>
        `;
        break;
      case 'bimodal':
        html = `
          <div style="display: grid; grid-template-columns: 1fr 1fr; gap: 0.5rem; margin-bottom: 0.5rem;">
            <div>
              <label class="form-label" style="font-size: 0.75rem;">Peak 1 Mean:</label>
              <input id="tp_bim_m1" type="number" value="10" class="form-control" style="padding: 0.3rem 0.4rem;">
            </div>
            <div>
              <label class="form-label" style="font-size: 0.75rem;">Peak 1 SD:</label>
              <input id="tp_bim_s1" type="number" min="0.1" value="2" class="form-control" style="padding: 0.3rem 0.4rem;">
            </div>
          </div>
          <div style="display: grid; grid-template-columns: 1fr 1fr; gap: 0.5rem;">
            <div>
              <label class="form-label" style="font-size: 0.75rem;">Peak 2 Mean:</label>
              <input id="tp_bim_m2" type="number" value="25" class="form-control" style="padding: 0.3rem 0.4rem;">
            </div>
            <div>
              <label class="form-label" style="font-size: 0.75rem;">Peak 2 SD:</label>
              <input id="tp_bim_s2" type="number" min="0.1" value="3" class="form-control" style="padding: 0.3rem 0.4rem;">
            </div>
          </div>
        `;
        break;
      case 'uniform':
        html = `
          <div style="display: grid; grid-template-columns: 1fr 1fr; gap: 0.75rem;">
            <div>
              <label class="form-label" style="font-size: 0.78rem;">Lower Bound (a):</label>
              <input id="tp_uni_min" type="number" value="0" class="form-control" style="padding: 0.35rem 0.5rem;">
            </div>
            <div>
              <label class="form-label" style="font-size: 0.78rem;">Upper Bound (b):</label>
              <input id="tp_uni_max" type="number" value="20" class="form-control" style="padding: 0.35rem 0.5rem;">
            </div>
          </div>
        `;
        break;
      case 'poisson':
        html = `
          <div>
            <label class="form-label" style="font-size: 0.78rem;">Event Rate Parameter (λ):</label>
            <input id="tp_pois_lambda" type="number" min="0.5" step="0.5" value="5" class="form-control" style="padding: 0.35rem 0.5rem;">
            <div style="font-size: 0.7rem; color: var(--text-muted); margin-top: 0.2rem;">Mean = Variance = λ = 5.0</div>
          </div>
        `;
        break;
      case 'chiSquare':
        html = `
          <div>
            <label class="form-label" style="font-size: 0.78rem;">Degrees of Freedom (k):</label>
            <input id="tp_chi_df" type="number" min="1" max="50" step="1" value="4" class="form-control" style="padding: 0.35rem 0.5rem;">
            <div style="font-size: 0.7rem; color: var(--text-muted); margin-top: 0.2rem;">Mean = k = 4.0, Variance = 2k = 8.0</div>
          </div>
        `;
        break;
    }

    container.innerHTML = html;

    // Attach reactive input listeners to all inputs in container
    container.querySelectorAll('input').forEach(input => {
      input.addEventListener('input', () => this.runTeachingDistribution());
    });
  }

  runTeachingDistribution() {
    const distKey = document.getElementById('teachingDistSelect')?.value || 'normal';
    const n = parseInt(document.getElementById('teachingNRange')?.value) || 500;
    const meta = Teaching.distMetadata[distKey] || {};

    let data = [];
    const params = {};
    let theoMean = 0;
    let theoSD = 0;
    let theoSkew = '0.000';
    let theoKurt = '0.000';

    switch (distKey) {
      case 'normal': {
        const mean = parseFloat(document.getElementById('tp_norm_mean')?.value) || 12;
        const sd = Math.max(0.1, parseFloat(document.getElementById('tp_norm_sd')?.value) || 3);
        params.mean = mean;
        params.sd = sd;
        theoMean = mean;
        theoSD = sd;
        theoSkew = '0.000 (Symmetric)';
        theoKurt = '0.000 (Mesokurtic)';
        data = Teaching.generators.normal(n, mean, sd);
        break;
      }
      case 'studentsT': {
        const df = Math.max(1, parseInt(document.getElementById('tp_t_df')?.value) || 4);
        const scale = Math.max(0.1, parseFloat(document.getElementById('tp_t_scale')?.value) || 2);
        params.df = df;
        params.scale = scale;
        params.mean = 0;
        theoMean = 0;
        theoSD = df > 2 ? scale * Math.sqrt(df / (df - 2)) : NaN;
        theoSkew = '0.000';
        theoKurt = df > 4 ? (6 / (df - 4)).toFixed(2) + ' (Heavy tails)' : '∞ (Fat tails)';
        data = Teaching.generators.studentsT(n, df, 0, scale);
        break;
      }
      case 'logNormal': {
        const mu = parseFloat(document.getElementById('tp_log_mu')?.value) || 1.5;
        const sigma = Math.max(0.1, parseFloat(document.getElementById('tp_log_sigma')?.value) || 0.6);
        params.mu = mu;
        params.sigma = sigma;
        theoMean = Math.exp(mu + 0.5 * sigma * sigma);
        theoSD = Math.sqrt((Math.exp(sigma * sigma) - 1) * Math.exp(2 * mu + sigma * sigma));
        theoSkew = ((Math.exp(sigma * sigma) + 2) * Math.sqrt(Math.exp(sigma * sigma) - 1)).toFixed(2) + ' (Right skew)';
        theoKurt = 'Positive (Leptokurtic)';
        data = Teaching.generators.logNormal(n, mu, sigma);
        break;
      }
      case 'exponential': {
        const rate = Math.max(0.01, parseFloat(document.getElementById('tp_exp_rate')?.value) || 0.4);
        params.rate = rate;
        theoMean = 1 / rate;
        theoSD = 1 / rate;
        theoSkew = '2.000 (Strong right skew)';
        theoKurt = '6.000 (Heavy tails)';
        data = Teaching.generators.exponential(n, rate);
        break;
      }
      case 'bimodal': {
        const m1 = parseFloat(document.getElementById('tp_bim_m1')?.value) || 10;
        const s1 = Math.max(0.1, parseFloat(document.getElementById('tp_bim_s1')?.value) || 2);
        const m2 = parseFloat(document.getElementById('tp_bim_m2')?.value) || 25;
        const s2 = Math.max(0.1, parseFloat(document.getElementById('tp_bim_s2')?.value) || 3);
        params.m1 = m1;
        params.s1 = s1;
        params.m2 = m2;
        params.s2 = s2;
        params.p = 0.5;
        theoMean = 0.5 * (m1 + m2);
        theoSD = Math.sqrt(0.5 * (s1 * s1 + s2 * s2) + 0.25 * Math.pow(m1 - m2, 2));
        theoSkew = '~0.000 (Symmetric Peaks)';
        theoKurt = 'Platykurtic (Bimodal Trough)';
        data = Teaching.generators.bimodal(n, m1, s1, m2, s2, 0.5);
        break;
      }
      case 'uniform': {
        const min = parseFloat(document.getElementById('tp_uni_min')?.value) || 0;
        const max = parseFloat(document.getElementById('tp_uni_max')?.value) || 20;
        params.min = min;
        params.max = max;
        theoMean = (min + max) / 2;
        theoSD = Math.sqrt(Math.pow(max - min, 2) / 12);
        theoSkew = '0.000 (Flat)';
        theoKurt = '-1.200 (Platykurtic)';
        data = Teaching.generators.uniform(n, min, max);
        break;
      }
      case 'poisson': {
        const lambda = Math.max(0.1, parseFloat(document.getElementById('tp_pois_lambda')?.value) || 5);
        params.lambda = lambda;
        theoMean = lambda;
        theoSD = Math.sqrt(lambda);
        theoSkew = (1 / Math.sqrt(lambda)).toFixed(2);
        theoKurt = (1 / lambda).toFixed(2);
        data = Teaching.generators.poisson(n, lambda);
        break;
      }
      case 'chiSquare': {
        const df = Math.max(1, parseInt(document.getElementById('tp_chi_df')?.value) || 4);
        params.df = df;
        theoMean = df;
        theoSD = Math.sqrt(2 * df);
        theoSkew = Math.sqrt(8 / df).toFixed(2);
        theoKurt = (12 / df).toFixed(2);
        data = Teaching.generators.chiSquare(n, df);
        break;
      }
    }

    this.currentTeachingData = data;

    // Calculate empirical descriptive statistics
    const stats = Descriptive.calculate(data);

    // Update Metrics Cards
    const meanEl = document.getElementById('teachingSampleMean');
    if (meanEl) meanEl.innerText = stats.mean.toFixed(2);
    const theoMeanEl = document.getElementById('teachingTheoMean');
    if (theoMeanEl) theoMeanEl.innerText = `Theo: ${isNaN(theoMean) ? 'N/A' : theoMean.toFixed(2)}`;

    const sdEl = document.getElementById('teachingSampleSD');
    if (sdEl) sdEl.innerText = stats.sd.toFixed(2);
    const theoSDEl = document.getElementById('teachingTheoSD');
    if (theoSDEl) theoSDEl.innerText = `Theo: ${isNaN(theoSD) ? 'N/A' : theoSD.toFixed(2)}`;

    const skewEl = document.getElementById('teachingSkewness');
    if (skewEl) skewEl.innerText = stats.skewness.toFixed(2);
    const skewSubEl = document.getElementById('teachingSkewnessSub');
    if (skewSubEl) skewSubEl.innerText = stats.skewnessInterpretation;

    const kurtEl = document.getElementById('teachingKurtosis');
    if (kurtEl) kurtEl.innerText = stats.kurtosis.toFixed(2);
    const kurtSubEl = document.getElementById('teachingKurtosisSub');
    if (kurtSubEl) kurtSubEl.innerText = stats.kurtosisInterpretation;

    const normEl = document.getElementById('teachingNormality');
    const normSubEl = document.getElementById('teachingNormalitySub');
    if (normEl) {
      normEl.innerText = stats.normality.isNormal ? 'Normal ✓' : 'Non-Normal ✗';
      normEl.style.color = stats.normality.isNormal ? '#22c55e' : '#ef4444';
    }
    if (normSubEl) {
      normSubEl.innerText = `JB p = ${stats.normality.pValue.toFixed(4)}`;
    }

    // Update Context Card
    const titleEl = document.getElementById('teachingContextTitle');
    const bodyEl = document.getElementById('teachingContextBody');
    if (titleEl) titleEl.innerText = `${meta.name || 'Distribution'} [${meta.symbol || ''}]`;
    if (bodyEl) {
      bodyEl.innerHTML = `
        <strong>Clinical Example:</strong> ${meta.clinicalExample || ''}<br>
        <strong>Theoretical Behavior:</strong> ${meta.properties || ''}<br>
        <strong>Methodological Guideline:</strong> ${meta.testNote || ''}
      `;
    }

    // Update Data Preview (first 50 values)
    const previewEl = document.getElementById('teachingDataPreview');
    if (previewEl) {
      previewEl.value = data.slice(0, 50).map(v => v.toFixed(2)).join(', ') + (data.length > 50 ? ` ... (+${data.length - 50} more)` : '');
    }

    // Render Canvas
    const canvas = document.getElementById('teachingDistCanvas');
    if (canvas && this.engines['teachingDistCanvas']) {
      Plots.renderTeachingDistribution(
        this.engines['teachingDistCanvas'],
        data,
        distKey,
        params,
        `${meta.name || 'Distribution'} (Empirical N = ${n} vs Theoretical PDF)`
      );
    }

    // Cache result
    this.results.teachingDist = {
      name: meta.name || distKey,
      clinicalExample: meta.clinicalExample,
      n,
      mean: stats.mean,
      theoMean,
      sd: stats.sd,
      theoSD,
      skewness: stats.skewness,
      skewnessInterpretation: stats.skewnessInterpretation,
      theoSkewness: theoSkew,
      kurtosis: stats.kurtosis,
      kurtosisInterpretation: stats.kurtosisInterpretation,
      theoKurtosis: theoKurt,
      normality: stats.normality
    };

    if (!this.results.teaching) this.results.teaching = {};
    this.results.teaching.dist = this.results.teachingDist;
  }

  updateCltUI(summary) {
    if (!summary) return;

    // Update Metrics
    const kEl = document.getElementById('cltSamplesCount');
    if (kEl) kEl.innerText = summary.samplesDrawn.toLocaleString();

    const nSubEl = document.getElementById('cltNSub');
    if (nSubEl) nSubEl.innerText = `Each of size n = ${summary.sampleSize}`;

    const muEl = document.getElementById('cltTrueMu');
    if (muEl) muEl.innerText = summary.theoreticalMean.toFixed(2);

    const sigmaSubEl = document.getElementById('cltTrueSigmaSub');
    if (sigmaSubEl) sigmaSubEl.innerText = `Pop SD: ${summary.population.sd.toFixed(2)}`;

    const obsMeanEl = document.getElementById('cltObsMean');
    const meanDiffEl = document.getElementById('cltMeanDiffSub');
    if (obsMeanEl) {
      obsMeanEl.innerText = summary.observedMean !== null ? summary.observedMean.toFixed(2) : '--';
    }
    if (meanDiffEl) {
      meanDiffEl.innerText = summary.observedMean !== null
        ? `|Diff|: ${Math.abs(summary.observedMean - summary.theoreticalMean).toFixed(3)}`
        : '|Diff|: --';
    }

    const theoSEEl = document.getElementById('cltTheoSE');
    if (theoSEEl) theoSEEl.innerText = summary.theoreticalSE.toFixed(3);

    const obsSEEl = document.getElementById('cltObsSE');
    const seDiffEl = document.getElementById('cltSEDiffSub');
    if (obsSEEl) {
      obsSEEl.innerText = summary.observedSE !== null ? summary.observedSE.toFixed(3) : '--';
    }
    if (seDiffEl) {
      seDiffEl.innerText = summary.observedSE !== null
        ? `Diff: ${(summary.observedSE - summary.theoreticalSE).toFixed(3)}`
        : 'Formula: σ / √n';
    }

    const normEl = document.getElementById('cltNormality');
    const normSubEl = document.getElementById('cltNormalitySub');
    if (normEl) {
      if (summary.normality) {
        normEl.innerText = summary.normality.isNormal ? 'Normal ✓' : 'Non-Normal';
        normEl.style.color = summary.normality.isNormal ? '#22c55e' : '#f59e0b';
      } else {
        normEl.innerText = '--';
        normEl.style.color = 'inherit';
      }
    }
    if (normSubEl) {
      normSubEl.innerText = summary.normality ? `JB p = ${summary.normality.pValue.toFixed(3)}` : 'JB Test: --';
    }

    // Render Canvas Charts
    if (this.engines['teachingCltParentCanvas']) {
      Plots.renderCltParent(
        this.engines['teachingCltParentCanvas'],
        summary.population,
        summary.lastSample,
        'Parent Population Distribution'
      );
    }

    if (this.engines['teachingCltSamplingCanvas']) {
      Plots.renderCltSampling(
        this.engines['teachingCltSamplingCanvas'],
        summary,
        'Sampling Distribution of Mean (x̄)'
      );
    }

    // Update Report Text
    const reportEl = document.getElementById('teachingReportText');
    if (reportEl) {
      if (summary.samplesDrawn === 0) {
        reportEl.innerText = `Simulation initialized with ${summary.population.name} and sample size n = ${summary.sampleSize}. Click "▶ Step (Draw 1 Sample)" or "⚡ Draw 100" to begin demonstrating convergence.`;
      } else {
        const normStatement = summary.normality?.isNormal
          ? 'conforms strictly to a Gaussian bell curve (Jarque-Bera p ≥ 0.05), mathematically proving the Central Limit Theorem.'
          : 'is rapidly converging toward Gaussian symmetry as iterations accumulate.';

        reportEl.innerText =
          `A Central Limit Theorem simulation was conducted using a ${summary.population.name} (True μ = ${summary.theoreticalMean.toFixed(2)}, σ = ${summary.population.sd.toFixed(2)}). ` +
          `A total of k = ${summary.samplesDrawn.toLocaleString()} independent samples (each of size n = ${summary.sampleSize}) were drawn. ` +
          `The grand mean of sample means was x̄̄ = ${summary.observedMean.toFixed(2)} (bias = ${Math.abs(summary.observedMean - summary.theoreticalMean).toFixed(3)}), ` +
          `with an empirical Standard Error of s_x̄ = ${summary.observedSE.toFixed(3)} compared to theoretical SE = σ/√n = ${summary.theoreticalSE.toFixed(3)}. ` +
          `The sampling distribution displays Skewness G₁ = ${summary.skewness.toFixed(3)} and ${normStatement} ` +
          `In clinical trials, this proves why parametric t-tests and ANOVA remain robust for sample sizes n ≥ 30 even when raw clinical metrics are non-normal.`;
      }
    }

    // Cache
    this.results.teachingClt = summary;
    if (!this.results.teaching) this.results.teaching = {};
    this.results.teaching.clt = summary;
  }

  runTeaching() {
    this.runTeachingDistribution();
    this.updateCltUI(Teaching.clt.getSummary());
  }

