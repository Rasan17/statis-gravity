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
      'rocCanvas'
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

    // 4. Categorical Events
    const matrixInputs = ['catA', 'catB', 'catC', 'catD'];
    matrixInputs.forEach(id => {
      document.getElementById(id)?.addEventListener('input', () => this.runCategorical());
    });
    document.getElementById('catSampleBtn')?.addEventListener('click', () => {
      const s = DataParser.samples.shuntComplications;
      document.getElementById('catA').value = s.a;
      document.getElementById('catB').value = s.b;
      document.getElementById('catC').value = s.c;
      document.getElementById('catD').value = s.d;
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
    document.getElementById('powerComputeBtn')?.addEventListener('click', () => this.runPower());

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

    document.querySelectorAll('.btn-save-plot').forEach(btn => {
      btn.addEventListener('click', (e) => {
        const canvasId = e.currentTarget.dataset.canvasId;
        if (this.engines[canvasId]) {
          this.engines[canvasId].saveAsImage(`${canvasId}.png`);
        }
      });
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

    // Box Plot
    if (this.engines.anovaCanvas) {
      Plots.renderBoxPlot(this.engines.anovaCanvas, res.groups, 'Multi-Cohort Comparison');
    }
  }

  runCategorical() {
    const a = parseFloat(document.getElementById('catA')?.value) || 0;
    const b = parseFloat(document.getElementById('catB')?.value) || 0;
    const c = parseFloat(document.getElementById('catC')?.value) || 0;
    const d = parseFloat(document.getElementById('catD')?.value) || 0;

    const res = Categorical.twoByTwo(a, b, c, d);
    if (res.error) return;

    const r = res.riskMetrics;
    document.getElementById('catOR').innerText = `${r.oddsRatio.toFixed(2)} [${r.orCI95[0].toFixed(2)}, ${r.orCI95[1].toFixed(2)}]`;
    document.getElementById('catRR').innerText = `${r.relativeRisk.toFixed(2)} [${r.rrCI95[0].toFixed(2)}, ${r.rrCI95[1].toFixed(2)}]`;
    document.getElementById('catChiSq').innerText = res.chiSquare.standard.toFixed(2);
    document.getElementById('catPVal').innerText = Exporter.formatP(res.chiSquare.pValueStandard);
    document.getElementById('catFisher').innerText = Exporter.formatP(res.fishersExact.pValue);
    document.getElementById('catNNT').innerText = isFinite(r.nnt) ? r.nnt.toFixed(1) : '∞';

    document.getElementById('catReportText').innerText = Exporter.formatCategoricalReport(res);
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
  }

  runPower() {
    const m1 = parseFloat(document.getElementById('pwrM1')?.value) || 10;
    const m2 = parseFloat(document.getElementById('pwrM2')?.value) || 15;
    const sd = parseFloat(document.getElementById('pwrSD')?.value) || 10;
    const alpha = parseFloat(document.getElementById('pwrAlpha')?.value) || 0.05;
    const power = parseFloat(document.getElementById('pwrPower')?.value) || 0.80;

    const res = PowerAnalysis.sampleSizeMeans(m1, m2, sd, alpha, power);
    if (res.error) {
      alert(res.error);
      return;
    }

    document.getElementById('pwrD').innerText = res.cohensD.toFixed(2);
    document.getElementById('pwrNGroup').innerText = res.nPerGroup;
    document.getElementById('pwrNTotal').innerText = res.totalN;

    document.getElementById('pwrReportText').innerText =
      `To detect an effect size of Cohen's d = ${res.cohensD.toFixed(2)} (|Δ| = ${res.diff.toFixed(2)}, SD = ${res.sd.toFixed(2)}) ` +
      `with ${Math.round(res.power * 100)}% power at two-sided α = ${res.alpha}, a minimum sample size of ` +
      `n = ${res.nPerGroup} patients per group (total N = ${res.totalN}) is required.`;
  }
}

// Instantiate on DOM load
window.addEventListener('DOMContentLoaded', () => {
  new StatisGravityApp();
});
