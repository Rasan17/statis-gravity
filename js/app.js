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
import { Randomiser } from './stats/randomiser.js';
import { Psm } from './stats/psm.js';

import { ChartEngine } from './visualization/chart-engine.js';
import { Plots } from './visualization/plots.js';

import { DataParser } from './data/parser.js';
import { Exporter } from './data/exporter.js';

class StatisGravityApp {
  constructor() {
    this.currentTheme = localStorage.getItem('sg_theme') || 'dark';
    this.engines = {};
    this.results = { teaching: {} };
    this.init();
  }

  init() {
    this.applyTheme(this.currentTheme);
    this.setupTabs();
    this.initAnovaGroups();
    this.setupEventListeners();
    this.initChartEngines();
    this.loadInitialSamples();
    this.initRandomiser();
    this.initPsm();
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
            const eng = this.engines[canvas.id];
            if (canvas && eng) {
              eng.initHiDPI();
              if (typeof eng.lastRender === 'function') eng.lastRender();
              else if (typeof eng.lastRenderFn === 'function') eng.lastRenderFn();
            }
          });
          if (targetId === 'tab-teaching') {
            if (typeof this.runTwoSampleOverlap === 'function') this.runTwoSampleOverlap();
            if (typeof this.runPowerSimulation === 'function') this.runPowerSimulation();
          }
          if (targetId === 'tab-teaching-bayesian') {
            if (typeof this.runBayesianSimulation === 'function') this.runBayesianSimulation();
          }
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
      'teachingCltSamplingCanvas',
      'teachingTCanvas',
      'teachingOverlapCanvas',
      'teachingPowerCanvas',
      'teachingBayesCanvas',
      'psmLovePlotCanvas',
      'psmOverlapCanvas'
    ];

    canvasIds.forEach(id => {
      const el = document.getElementById(id);
      if (el) {
        const opts = id === 'teachingBayesCanvas'
          ? { theme: this.currentTheme, padding: { top: 20, right: 20, bottom: 24, left: 20 } }
          : { theme: this.currentTheme };
        this.engines[id] = new ChartEngine(el, opts);
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
    document.getElementById('hypoTestType')?.addEventListener('change', () => this.runHypothesis());
    document.getElementById('hypoErrorBarMode')?.addEventListener('change', () => this.runHypothesis());

    const pairedCb = document.getElementById('hypoIsPaired');
    pairedCb?.addEventListener('change', (e) => {
      const isPaired = e.target.checked;
      const modeBadge = document.getElementById('hypoDesignModeBadge');
      if (modeBadge) {
        modeBadge.innerText = isPaired 
          ? 'Paired / Dependent Samples (Pre vs. Post / Matched)' 
          : 'Independent Cohorts (Two Separate Groups)';
        modeBadge.className = `badge ${isPaired ? 'badge-sig' : 'badge-neutral'}`;
      }
      const sampleBtn = document.getElementById('hypoSampleBtn');
      if (sampleBtn) {
        sampleBtn.innerText = isPaired ? 'Load Pre/Post ICP' : 'Load Independent Cohorts';
      }
      const nameAInput = document.getElementById('hypoNameA');
      const nameBInput = document.getElementById('hypoNameB');
      if (nameAInput && nameBInput) {
        if (isPaired && (nameAInput.value === 'Cohort 1 (Control)' || nameAInput.value === 'Standard Resection' || nameAInput.value === 'Group A')) {
          nameAInput.value = 'Pre-Infusion';
          nameBInput.value = 'Post-Infusion';
        } else if (!isPaired && (nameAInput.value === 'Pre-Infusion' || nameAInput.value === 'Pre-Intervention')) {
          nameAInput.value = 'Standard Resection';
          nameBInput.value = 'Supramarginal Resection';
        }
      }
      if (document.getElementById('hypoGroupA')?.value && document.getElementById('hypoGroupB')?.value) {
        this.runHypothesis();
      }
    });

    document.getElementById('hypoSampleBtn')?.addEventListener('click', () => {
      const isPaired = document.getElementById('hypoIsPaired')?.checked || false;
      if (isPaired) {
        document.getElementById('hypoNameA').value = 'Pre-Infusion';
        document.getElementById('hypoNameB').value = 'Post-Infusion';
        const s = DataParser.samples.icpDynamics;
        const ga = s.groupA.data || s.groupA;
        const gb = s.groupB.data || s.groupB;
        document.getElementById('hypoGroupA').value = ga.join(', ');
        document.getElementById('hypoGroupB').value = gb.join(', ');
      } else {
        document.getElementById('hypoNameA').value = 'Standard Resection';
        document.getElementById('hypoNameB').value = 'Supramarginal Resection';
        document.getElementById('hypoGroupA').value = [14.2, 15.1, 13.8, 16.5, 14.9, 15.8, 17.2, 13.5, 15.0, 14.6, 16.1, 14.8, 15.4, 16.0].join(', ');
        document.getElementById('hypoGroupB').value = [22.4, 28.1, 18.9, 31.5, 24.0, 19.8, 35.2, 26.7, 21.3, 29.4, 33.1, 20.5].join(', ');
      }
      this.runHypothesis();
    });

    const designModal = document.getElementById('hypoDesignModal');
    document.getElementById('hypoDesignInfoBtn')?.addEventListener('click', () => {
      if (designModal) designModal.style.display = 'flex';
    });
    document.getElementById('hypoDesignModalClose')?.addEventListener('click', () => {
      if (designModal) designModal.style.display = 'none';
    });
    document.getElementById('hypoDesignModalSetIndependent')?.addEventListener('click', () => {
      if (pairedCb) { pairedCb.checked = false; pairedCb.dispatchEvent(new Event('change')); }
      if (designModal) designModal.style.display = 'none';
    });
    document.getElementById('hypoDesignModalSetPaired')?.addEventListener('click', () => {
      if (pairedCb) { pairedCb.checked = true; pairedCb.dispatchEvent(new Event('change')); }
      if (designModal) designModal.style.display = 'none';
    });
    designModal?.addEventListener('click', (e) => {
      if (e.target === designModal) designModal.style.display = 'none';
    });

    // 3. ANOVA & Multi-Group Analysis
    document.getElementById('anovaComputeBtn')?.addEventListener('click', () => this.runAnova());
    document.getElementById('anovaAddGroupBtn')?.addEventListener('click', () => this.addAnovaGroup());
    
    const anovaPairedCb = document.getElementById('anovaIsPaired');
    anovaPairedCb?.addEventListener('change', () => {
      this.updateAnovaDesignUI();
      this.runAnova();
    });

    document.getElementById('anovaTestType')?.addEventListener('change', () => this.runAnova());
    document.getElementById('anovaErrorBarMode')?.addEventListener('change', () => this.runAnova());

    // ANOVA Presets
    document.getElementById('anovaSampleBtn')?.addEventListener('click', () => this.loadAnovaPreset('sample3'));
    document.getElementById('anovaSampleWelchBtn')?.addEventListener('click', () => this.loadAnovaPreset('sampleWelch'));
    document.getElementById('anovaSampleRMBtn')?.addEventListener('click', () => this.loadAnovaPreset('sampleRM'));
    document.getElementById('anovaSampleSkewBtn')?.addEventListener('click', () => this.loadAnovaPreset('sampleSkew'));

    // ANOVA Study Design Modal
    const anovaModal = document.getElementById('anovaDesignModal');
    document.getElementById('anovaDesignInfoBtn')?.addEventListener('click', () => {
      if (anovaModal) anovaModal.style.display = 'flex';
    });
    document.getElementById('anovaDesignModalClose')?.addEventListener('click', () => {
      if (anovaModal) anovaModal.style.display = 'none';
    });
    document.getElementById('anovaDesignModalSetIndependent')?.addEventListener('click', () => {
      if (anovaPairedCb) { anovaPairedCb.checked = false; anovaPairedCb.dispatchEvent(new Event('change')); }
      if (anovaModal) anovaModal.style.display = 'none';
    });
    document.getElementById('anovaDesignModalSetPaired')?.addEventListener('click', () => {
      if (anovaPairedCb) { anovaPairedCb.checked = true; anovaPairedCb.dispatchEvent(new Event('change')); }
      if (anovaModal) anovaModal.style.display = 'none';
    });
    anovaModal?.addEventListener('click', (e) => {
      if (e.target === anovaModal) anovaModal.style.display = 'none';
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

    // Student's t Convergence Controls
    document.getElementById('tConvNRange')?.addEventListener('input', (e) => {
      this.runTConvergence(parseInt(e.target.value));
    });

    document.querySelectorAll('.btn-t-preset').forEach(btn => {
      btn.addEventListener('click', (e) => {
        const n = parseInt(e.currentTarget.dataset.n);
        if (n) this.runTConvergence(n);
      });
    });

    document.getElementById('tConvAnimateBtn')?.addEventListener('click', () => {
      this.animateTConvergence();
    });

    document.getElementById('tConvResetBtn')?.addEventListener('click', () => {
      if (this.tConvAnimationTimer) {
        clearInterval(this.tConvAnimationTimer);
        this.tConvAnimationTimer = null;
        const btn = document.getElementById('tConvAnimateBtn');
        if (btn) btn.innerText = '▶ Animate Convergence';
      }
      this.runTConvergence(4);
    });

    document.getElementById('tConvShowTailArea')?.addEventListener('change', () => {
      this.runTConvergence();
    });

    // Section 4: Two-Sample Overlap, Group SD/SEM & Alpha Boundary Events
    const syncGroup1FromSD = () => {
      const sd1 = parseFloat(document.getElementById('overlapSD1Range')?.value) || 2.5;
      const n1 = parseInt(document.getElementById('overlapN1Range')?.value) || 16;
      const sem1 = sd1 / Math.sqrt(n1);
      const sem1El = document.getElementById('overlapSEM1Range');
      if (sem1El) sem1El.value = sem1.toFixed(3);
      if (document.getElementById('overlapLinkGroupsCheck')?.checked) {
        const sd2El = document.getElementById('overlapSD2Range');
        const sem2El = document.getElementById('overlapSEM2Range');
        if (sd2El) sd2El.value = sd1;
        if (sem2El) sem2El.value = sem1.toFixed(3);
      }
      this.runTwoSampleOverlap();
    };

    const syncGroup1FromSEM = () => {
      const sd1 = parseFloat(document.getElementById('overlapSD1Range')?.value) || 2.5;
      const sem1 = parseFloat(document.getElementById('overlapSEM1Range')?.value) || 0.625;
      const n1 = Math.max(2, Math.min(500, Math.round(Math.pow(sd1 / sem1, 2))));
      const n1El = document.getElementById('overlapN1Range');
      if (n1El) n1El.value = n1;
      if (document.getElementById('overlapLinkGroupsCheck')?.checked) {
        const sem2El = document.getElementById('overlapSEM2Range');
        const n2El = document.getElementById('overlapN2Range');
        if (sem2El) sem2El.value = sem1.toFixed(3);
        if (n2El) n2El.value = n1;
      }
      this.runTwoSampleOverlap();
    };

    const syncGroup1FromN = () => {
      const sd1 = parseFloat(document.getElementById('overlapSD1Range')?.value) || 2.5;
      const n1 = parseInt(document.getElementById('overlapN1Range')?.value) || 16;
      const sem1 = sd1 / Math.sqrt(n1);
      const sem1El = document.getElementById('overlapSEM1Range');
      if (sem1El) sem1El.value = sem1.toFixed(3);
      if (document.getElementById('overlapLinkGroupsCheck')?.checked) {
        const n2El = document.getElementById('overlapN2Range');
        const sem2El = document.getElementById('overlapSEM2Range');
        if (n2El) n2El.value = n1;
        if (sem2El) sem2El.value = sem1.toFixed(3);
      }
      this.runTwoSampleOverlap();
    };

    const syncGroup2FromSD = () => {
      const linkCheck = document.getElementById('overlapLinkGroupsCheck');
      if (linkCheck) linkCheck.checked = false;
      const sd2 = parseFloat(document.getElementById('overlapSD2Range')?.value) || 2.5;
      const n2 = parseInt(document.getElementById('overlapN2Range')?.value) || 16;
      const sem2 = sd2 / Math.sqrt(n2);
      const sem2El = document.getElementById('overlapSEM2Range');
      if (sem2El) sem2El.value = sem2.toFixed(3);
      this.runTwoSampleOverlap();
    };

    const syncGroup2FromSEM = () => {
      const linkCheck = document.getElementById('overlapLinkGroupsCheck');
      if (linkCheck) linkCheck.checked = false;
      const sd2 = parseFloat(document.getElementById('overlapSD2Range')?.value) || 2.5;
      const sem2 = parseFloat(document.getElementById('overlapSEM2Range')?.value) || 0.625;
      const n2 = Math.max(2, Math.min(500, Math.round(Math.pow(sd2 / sem2, 2))));
      const n2El = document.getElementById('overlapN2Range');
      if (n2El) n2El.value = n2;
      this.runTwoSampleOverlap();
    };

    const syncGroup2FromN = () => {
      const linkCheck = document.getElementById('overlapLinkGroupsCheck');
      if (linkCheck) linkCheck.checked = false;
      const sd2 = parseFloat(document.getElementById('overlapSD2Range')?.value) || 2.5;
      const n2 = parseInt(document.getElementById('overlapN2Range')?.value) || 16;
      const sem2 = sd2 / Math.sqrt(n2);
      const sem2El = document.getElementById('overlapSEM2Range');
      if (sem2El) sem2El.value = sem2.toFixed(3);
      this.runTwoSampleOverlap();
    };

    document.getElementById('overlapSD1Range')?.addEventListener('input', syncGroup1FromSD);
    document.getElementById('overlapSEM1Range')?.addEventListener('input', syncGroup1FromSEM);
    document.getElementById('overlapN1Range')?.addEventListener('input', syncGroup1FromN);

    document.getElementById('overlapSD2Range')?.addEventListener('input', syncGroup2FromSD);
    document.getElementById('overlapSEM2Range')?.addEventListener('input', syncGroup2FromSEM);
    document.getElementById('overlapN2Range')?.addEventListener('input', syncGroup2FromN);

    document.getElementById('overlapDeltaRange')?.addEventListener('input', () => this.runTwoSampleOverlap());
    document.getElementById('overlapAlphaRange')?.addEventListener('input', () => this.runTwoSampleOverlap());

    document.getElementById('overlapLinkGroupsCheck')?.addEventListener('change', (e) => {
      if (e.target.checked) {
        const sd1 = document.getElementById('overlapSD1Range')?.value;
        const sem1 = document.getElementById('overlapSEM1Range')?.value;
        const n1 = document.getElementById('overlapN1Range')?.value;
        const sd2El = document.getElementById('overlapSD2Range');
        const sem2El = document.getElementById('overlapSEM2Range');
        const n2El = document.getElementById('overlapN2Range');
        if (sd2El && sd1) sd2El.value = sd1;
        if (sem2El && sem1) sem2El.value = sem1;
        if (n2El && n1) n2El.value = n1;
        this.runTwoSampleOverlap();
      }
    });

    document.querySelectorAll('.btn-overlap-alpha').forEach(btn => {
      btn.addEventListener('click', (e) => {
        const alpha = parseFloat(e.currentTarget.dataset.alpha);
        const range = document.getElementById('overlapAlphaRange');
        if (range && !isNaN(alpha)) {
          range.value = alpha;
          document.querySelectorAll('.btn-overlap-alpha').forEach(b => {
            b.style.borderColor = '';
            b.style.color = '';
            b.style.fontWeight = '';
          });
          e.currentTarget.style.borderColor = '#f59e0b';
          e.currentTarget.style.color = '#f59e0b';
          e.currentTarget.style.fontWeight = '600';
          this.runTwoSampleOverlap();
        }
      });
    });

    document.querySelectorAll('.btn-overlap-mode').forEach(btn => {
      btn.addEventListener('click', (e) => {
        const mode = e.currentTarget.dataset.mode;
        this.currentOverlapMode = mode;
        document.querySelectorAll('.btn-overlap-mode').forEach(b => {
          b.classList.remove('btn-primary');
          b.classList.add('btn-secondary');
        });
        e.currentTarget.classList.remove('btn-secondary');
        e.currentTarget.classList.add('btn-primary');
        this.runTwoSampleOverlap();
      });
    });

    document.getElementById('overlapAnimateBtn')?.addEventListener('click', () => {
      this.animateOverlapSeparation();
    });

    document.getElementById('overlapResetBtn')?.addEventListener('click', () => {
      if (this.overlapAnimationTimer) {
        clearInterval(this.overlapAnimationTimer);
        this.overlapAnimationTimer = null;
        const btn = document.getElementById('overlapAnimateBtn');
        if (btn) btn.innerText = '▶ Animate';
      }
      const dRange = document.getElementById('overlapDeltaRange');
      const aRange = document.getElementById('overlapAlphaRange');
      const sd1Range = document.getElementById('overlapSD1Range');
      const sem1Range = document.getElementById('overlapSEM1Range');
      const n1Range = document.getElementById('overlapN1Range');
      const sd2Range = document.getElementById('overlapSD2Range');
      const sem2Range = document.getElementById('overlapSEM2Range');
      const n2Range = document.getElementById('overlapN2Range');
      const linkCheck = document.getElementById('overlapLinkGroupsCheck');
      if (dRange) dRange.value = 2.0;
      if (aRange) aRange.value = 0.050;
      if (sd1Range) sd1Range.value = 2.5;
      if (sem1Range) sem1Range.value = 0.625;
      if (n1Range) n1Range.value = 16;
      if (sd2Range) sd2Range.value = 2.5;
      if (sem2Range) sem2Range.value = 0.625;
      if (n2Range) n2Range.value = 16;
      if (linkCheck) linkCheck.checked = false;
      this.currentOverlapMode = 'means';
      document.querySelectorAll('.btn-overlap-mode').forEach(b => {
        b.classList.remove('btn-primary');
        b.classList.add('btn-secondary');
        if (b.dataset.mode === 'means') {
          b.classList.remove('btn-secondary');
          b.classList.add('btn-primary');
        }
      });
      document.querySelectorAll('.btn-overlap-alpha').forEach(b => {
        b.style.borderColor = '';
        b.style.color = '';
        b.style.fontWeight = '';
        if (b.dataset.alpha === '0.05') {
          b.style.borderColor = '#f59e0b';
          b.style.color = '#f59e0b';
          b.style.fontWeight = '600';
        }
      });
      this.runTwoSampleOverlap();
    });

    // Teaching Simulation 5: Power Simulation Listeners
    document.getElementById('powerSDRange')?.addEventListener('input', (e) => {
      this.runPowerSimulation({ sd: parseFloat(e.target.value) }, 'sd');
    });

    document.getElementById('powerSEMRange')?.addEventListener('input', (e) => {
      this.runPowerSimulation({ sem: parseFloat(e.target.value) }, 'sem');
    });

    document.getElementById('powerNRange')?.addEventListener('input', (e) => {
      this.runPowerSimulation({ n: parseInt(e.target.value) }, 'n');
    });

    document.getElementById('powerPowerRange')?.addEventListener('input', (e) => {
      this.runPowerSimulation({ power: parseFloat(e.target.value) }, 'power');
    });

    document.getElementById('powerBetaRange')?.addEventListener('input', (e) => {
      this.runPowerSimulation({ beta: parseFloat(e.target.value) }, 'beta');
    });

    document.getElementById('powerDeltaRange')?.addEventListener('input', (e) => {
      this.runPowerSimulation({ delta: parseFloat(e.target.value) }, 'delta');
    });

    document.getElementById('powerAlphaRange')?.addEventListener('input', (e) => {
      this.runPowerSimulation({ alpha: parseFloat(e.target.value) }, 'alpha');
    });

    document.querySelectorAll('.btn-power-preset').forEach(btn => {
      btn.addEventListener('click', (e) => {
        document.querySelectorAll('.btn-power-preset').forEach(b => b.classList.remove('active'));
        e.currentTarget.classList.add('active');
        const pwr = parseFloat(e.currentTarget.dataset.power);
        this.runPowerSimulation({ power: pwr }, 'power');
      });
    });

    document.querySelectorAll('.btn-power-effect-preset').forEach(btn => {
      btn.addEventListener('click', (e) => {
        document.querySelectorAll('.btn-power-effect-preset').forEach(b => b.classList.remove('active'));
        e.currentTarget.classList.add('active');
        const d = parseFloat(e.currentTarget.dataset.d);
        const sd = parseFloat(document.getElementById('powerSDRange')?.value) || 4.0;
        const delta = d * sd;
        this.runPowerSimulation({ delta }, 'delta');
      });
    });

    document.querySelectorAll('.btn-power-mode').forEach(btn => {
      btn.addEventListener('click', (e) => {
        document.querySelectorAll('.btn-power-mode').forEach(b => {
          b.classList.remove('btn-primary', 'active');
          b.classList.add('btn-secondary');
        });
        e.currentTarget.classList.remove('btn-secondary');
        e.currentTarget.classList.add('btn-primary', 'active');
        this.currentPowerMode = e.currentTarget.dataset.mode;
        this.runPowerSimulation();
      });
    });

    document.getElementById('powerAnimateBtn')?.addEventListener('click', () => {
      this.animatePowerGain();
    });

    document.getElementById('powerResetBtn')?.addEventListener('click', () => {
      if (this.powerAnimationTimer) {
        clearInterval(this.powerAnimationTimer);
        this.powerAnimationTimer = null;
        const btn = document.getElementById('powerAnimateBtn');
        if (btn) btn.innerText = '▶ Animate Power Gain';
      }
      const sdRange = document.getElementById('powerSDRange');
      const semRange = document.getElementById('powerSEMRange');
      const nRange = document.getElementById('powerNRange');
      const pRange = document.getElementById('powerPowerRange');
      const bRange = document.getElementById('powerBetaRange');
      const dRange = document.getElementById('powerDeltaRange');
      const aRange = document.getElementById('powerAlphaRange');
      if (sdRange) sdRange.value = 4.0;
      if (semRange) semRange.value = 0.500;
      if (nRange) nRange.value = 64;
      if (pRange) pRange.value = 0.80;
      if (bRange) bRange.value = 0.20;
      if (dRange) dRange.value = 2.0;
      if (aRange) aRange.value = 0.050;
      this.currentPowerMode = 'distributions';
      document.querySelectorAll('.btn-power-mode').forEach(b => {
        b.classList.remove('btn-primary', 'active');
        b.classList.add('btn-secondary');
        if (b.dataset.mode === 'distributions') {
          b.classList.remove('btn-secondary');
          b.classList.add('btn-primary', 'active');
        }
      });
      document.querySelectorAll('.btn-power-preset').forEach(b => {
        b.classList.toggle('active', b.dataset.power === '0.80');
      });
      document.querySelectorAll('.btn-power-effect-preset').forEach(b => {
        b.classList.toggle('active', b.dataset.d === '0.50');
      });
      this.runPowerSimulation({ power: 0.80, sd: 4.0, delta: 2.0, alpha: 0.05 }, 'power');
    });

    // Bayesian Simulation (Section 6) Events
    document.getElementById('bayesPriorRange')?.addEventListener('input', (e) => {
      this.runBayesianSimulation({ prior: parseFloat(e.target.value) / 100 });
    });

    document.getElementById('bayesSampleRange')?.addEventListener('input', (e) => {
      this.runBayesianSimulation({ sampleSize: parseInt(e.target.value) });
    });

    document.getElementById('bayesLikelihoodRange')?.addEventListener('input', (e) => {
      this.runBayesianSimulation({ likelihood: parseFloat(e.target.value) / 100 });
    });

    document.getElementById('bayesFalsePosRange')?.addEventListener('input', (e) => {
      this.runBayesianSimulation({ falsePositive: parseFloat(e.target.value) / 100 });
    });

    // Preset buttons
    document.querySelectorAll('.teaching-bayes-preset-btn').forEach(btn => {
      btn.addEventListener('click', (e) => {
        document.querySelectorAll('.teaching-bayes-preset-btn').forEach(b => b.classList.remove('active'));
        e.currentTarget.classList.add('active');
        const presetKey = e.currentTarget.dataset.preset;
        this.currentBayesPreset = presetKey;
        let p = 0.0476, lik = 0.40, fp = 0.10, n = 210;
        if (presetKey === 'steve') {
          p = 0.0476; lik = 0.40; fp = 0.10; n = 210;
        } else if (presetKey === 'disease') {
          p = 0.01; lik = 0.95; fp = 0.05; n = 1000;
        } else if (presetKey === 'science') {
          p = 0.10; lik = 0.80; fp = 0.05; n = 500;
        } else if (presetKey === 'equal') {
          p = 0.20; lik = 0.50; fp = 0.50; n = 200;
        } else if (presetKey === 'fair_coin') {
          p = 0.50; lik = 0.75; fp = 0.25; n = 100;
        }
        this.runBayesianSimulation({ prior: p, likelihood: lik, falsePositive: fp, sampleSize: n, preset: presetKey });
      });
    });

    // View mode pills
    document.querySelectorAll('.teaching-bayes-view-btn').forEach(btn => {
      btn.addEventListener('click', (e) => {
        document.querySelectorAll('.teaching-bayes-view-btn').forEach(b => {
          b.classList.remove('btn-primary', 'active');
          b.classList.add('btn-secondary');
        });
        e.currentTarget.classList.remove('btn-secondary');
        e.currentTarget.classList.add('btn-primary', 'active');
        this.currentBayesView = e.currentTarget.dataset.view;
        this.runBayesianSimulation();
      });
    });

    // Animate restriction button
    document.getElementById('bayesAnimateBtn')?.addEventListener('click', () => {
      this.animateBayesRestriction();
    });

    // Sequential evidence button (Today's posterior becomes tomorrow's prior)
    document.getElementById('bayesNextEvidenceBtn')?.addEventListener('click', () => {
      if (this.results?.teaching?.bayes) {
        const nextPrior = this.results.teaching.bayes.posterior;
        const priorSlider = document.getElementById('bayesPriorRange');
        if (priorSlider) priorSlider.value = (nextPrior * 100).toFixed(1);
        this.runBayesianSimulation({ prior: nextPrior });
      }
    });

    // Reset button
    document.getElementById('bayesResetBtn')?.addEventListener('click', () => {
      if (this.bayesAnimationTimer) {
        clearInterval(this.bayesAnimationTimer);
        this.bayesAnimationTimer = null;
        const btn = document.getElementById('bayesAnimateBtn');
        if (btn) btn.innerText = '▶ Animate Restriction';
      }
      this.currentBayesPreset = 'steve';
      this.currentBayesView = 'square';
      document.querySelectorAll('.teaching-bayes-preset-btn').forEach(b => {
        b.classList.toggle('active', b.dataset.preset === 'steve');
      });
      document.querySelectorAll('.teaching-bayes-view-btn').forEach(b => {
        b.classList.remove('btn-primary', 'active');
        b.classList.add('btn-secondary');
        if (b.dataset.view === 'square') {
          b.classList.remove('btn-secondary');
          b.classList.add('btn-primary', 'active');
        }
      });
      this.runBayesianSimulation({ prior: 0.0476, likelihood: 0.40, falsePositive: 0.10, sampleSize: 210, preset: 'steve' });
    });

    // Double-Click Inline Number Editing Configuration
    const editableConfigs = [
      // Section 6: Bayesian
      { span: 'bayesPriorVal', range: 'bayesPriorRange' },
      { span: 'bayesStatusPriorVal', range: 'bayesPriorRange' },
      { span: 'bayesSampleVal', range: 'bayesSampleRange' },
      { span: 'bayesLikelihoodVal', range: 'bayesLikelihoodRange' },
      { span: 'bayesStatusLikelihoodVal', range: 'bayesLikelihoodRange' },
      { span: 'bayesFalsePosVal', range: 'bayesFalsePosRange' },
      { span: 'bayesStatusFalsePosVal', range: 'bayesFalsePosRange' },

      // Teaching: Sample Size Distribution
      { span: 'teachingNVal', range: 'teachingNRange' },

      // Teaching: CLT
      { span: 'cltNVal', range: 'cltNRange' },

      // Teaching: Student's t
      { span: 'tConvNVal', range: 'tConvNRange' },

      // Teaching: Overlap
      { span: 'overlapDeltaVal', range: 'overlapDeltaRange' },
      { span: 'overlapAlphaVal', range: 'overlapAlphaRange' },
      { span: 'overlapSD1Val', range: 'overlapSD1Range' },
      { span: 'overlapSEM1Val', range: 'overlapSEM1Range' },
      { span: 'overlapN1Val', range: 'overlapN1Range' },
      { span: 'overlapSD2Val', range: 'overlapSD2Range' },
      { span: 'overlapSEM2Val', range: 'overlapSEM2Range' },
      { span: 'overlapN2Val', range: 'overlapN2Range' },

      // Teaching: Power Simulation
      { span: 'powerSDVal', range: 'powerSDRange' },
      { span: 'powerSEMVal', range: 'powerSEMRange' },
      { span: 'powerNVal', range: 'powerNRange' },
      { span: 'powerPowerVal', range: 'powerPowerRange' },
      { span: 'powerBetaVal', range: 'powerBetaRange' },
      { span: 'powerDeltaVal', range: 'powerDeltaRange' },
      { span: 'powerAlphaVal', range: 'powerAlphaRange' }
    ];

    editableConfigs.forEach(cfg => {
      this.setupEditableNumber(cfg.span, cfg.range, cfg.options);
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

  setupEditableNumber(spanTarget, rangeTarget, options = {}) {
    const spanEl = typeof spanTarget === 'string' ? document.getElementById(spanTarget) : spanTarget;
    const rangeEl = typeof rangeTarget === 'string' ? document.getElementById(rangeTarget) : rangeTarget;
    if (!spanEl || !rangeEl) return;

    spanEl.classList.add('editable-number-badge');
    if (!spanEl.title) {
      spanEl.title = 'Double-click to edit value directly';
    }

    spanEl.addEventListener('dblclick', (e) => {
      e.preventDefault();
      e.stopPropagation();

      if (spanEl.dataset.isEditing === 'true') return;
      spanEl.dataset.isEditing = 'true';

      // Read current numeric value from the linked range input
      const currentVal = rangeEl.value;
      const min = rangeEl.min !== '' ? rangeEl.min : (options.min !== undefined ? options.min : '');
      const max = rangeEl.max !== '' ? rangeEl.max : (options.max !== undefined ? options.max : '');
      const step = rangeEl.step !== '' ? rangeEl.step : (options.step !== undefined ? options.step : 'any');

      // Create inline input - use step 'any' to allow arbitrary decimals like 11.8 or 4.76
      const input = document.createElement('input');
      input.type = 'number';
      input.className = 'editable-number-input';
      if (min !== '') input.min = min;
      if (max !== '') input.max = max;
      input.step = 'any';
      input.value = currentVal;

      // Match styles with span
      const computed = window.getComputedStyle(spanEl);
      input.style.fontSize = computed.fontSize;
      input.style.fontWeight = computed.fontWeight || '700';
      input.style.color = computed.color || '#00d2ff';
      input.style.textAlign = computed.textAlign === 'center' ? 'center' : 'right';
      const initialWidth = Math.max(68, spanEl.offsetWidth + 18);
      input.style.width = `${initialWidth}px`;

      // Hide span and insert input right beside it
      const prevDisplay = spanEl.style.display;
      spanEl.style.display = 'none';
      spanEl.parentNode.insertBefore(input, spanEl.nextSibling);

      // Focus and select all
      setTimeout(() => {
        input.focus();
        input.select();
      }, 10);

      let isFinished = false;
      const finishEdit = (commit) => {
        if (isFinished) return;
        isFinished = true;

        input.remove();
        spanEl.style.display = prevDisplay;
        spanEl.dataset.isEditing = 'false';

        if (commit) {
          const cleanedStr = input.value.toString().replace(/[^0-9.-]/g, '');
          let val = parseFloat(cleanedStr);
          if (!isNaN(val)) {
            const minNum = parseFloat(min);
            const maxNum = parseFloat(max);
            if (!isNaN(minNum) && val < minNum) val = minNum;
            if (!isNaN(maxNum) && val > maxNum) val = maxNum;

            // Prevent native <input type="range"> from snapping custom decimals to coarse step
            const prevStep = rangeEl.step;
            rangeEl.step = 'any';
            rangeEl.value = val;
            rangeEl.dispatchEvent(new Event('input', { bubbles: true }));
            rangeEl.dispatchEvent(new Event('change', { bubbles: true }));

            // When user later drags slider by mouse/touch, revert to original slider step
            const restoreSliderStep = () => {
              if (prevStep) rangeEl.step = prevStep;
              rangeEl.removeEventListener('pointerdown', restoreSliderStep);
              rangeEl.removeEventListener('touchstart', restoreSliderStep);
              rangeEl.removeEventListener('mousedown', restoreSliderStep);
            };
            rangeEl.addEventListener('pointerdown', restoreSliderStep, { once: true });
            rangeEl.addEventListener('touchstart', restoreSliderStep, { once: true });
            rangeEl.addEventListener('mousedown', restoreSliderStep, { once: true });

            if (typeof options.onCommit === 'function') {
              options.onCommit(val);
            }
          }
        }
      };

      input.addEventListener('keydown', (ev) => {
        if (ev.key === 'Enter') {
          ev.preventDefault();
          finishEdit(true);
        } else if (ev.key === 'Escape') {
          ev.preventDefault();
          finishEdit(false);
        }
      });

      input.addEventListener('blur', () => {
        finishEdit(true);
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
    this.initTeachingModule();
    this.loadPsmSample();
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
    const selectedTest = document.getElementById('hypoTestType')?.value || 'auto';
    const isPaired = document.getElementById('hypoIsPaired')?.checked || false;

    const a = DataParser.parseSeries(rawA);
    const b = DataParser.parseSeries(rawB);

    // Evaluate statistical assumptions
    const assumptions = Hypothesis.evaluateAssumptions(a, b, isPaired);
    if (assumptions.error) {
      alert(assumptions.error);
      return;
    }

    let testToRun = selectedTest;
    if (selectedTest === 'auto') {
      testToRun = assumptions.recommendedTest;
    }

    let res;
    if (testToRun === 'student') res = Hypothesis.independentTTest(a, b);
    else if (testToRun === 'welch') res = Hypothesis.welchTTest(a, b);
    else if (testToRun === 'paired') res = Hypothesis.pairedTTest(a, b);
    else if (testToRun === 'wilcoxon') res = Hypothesis.wilcoxonSignedRank(a, b);
    else res = Hypothesis.mannWhitneyUTest(a, b);

    if (res.error) {
      alert(res.error);
      return;
    }

    res.groupA = Object.assign(res.groupA || Descriptive.calculate(a), { name: nameA });
    res.groupB = Object.assign(res.groupB || Descriptive.calculate(b), { name: nameB });
    res.assumptions = assumptions;
    res.isPaired = isPaired;
    res.selectedTest = selectedTest;
    res.actualTest = testToRun;

    document.getElementById('hypoStat').innerText = (res.statistic !== undefined ? res.statistic : (res.zScore || 0)).toFixed(2);
    document.getElementById('hypoPVal').innerText = Exporter.formatP(res.pValue);
    const badge = document.getElementById('hypoPValBadge');
    if (badge) {
      badge.className = `badge ${res.isSignificant ? 'badge-sig' : 'badge-ns'}`;
      badge.innerText = res.isSignificant ? 'Significant (p < .05)' : 'Not Significant (ns)';
    }
    document.getElementById('hypoEffect').innerText = (res.cohensD !== undefined ? res.cohensD.toFixed(2) : (res.rankBiserial !== undefined ? res.rankBiserial.toFixed(2) : '0.00'));
    document.getElementById('hypoDiff').innerText = res.meanDiff !== undefined ? ((res.meanDiff >= 0 ? '+' : '') + res.meanDiff.toFixed(2)) : (res.medianDiff !== undefined ? ((res.medianDiff >= 0 ? '+' : '') + res.medianDiff.toFixed(2)) : 'N/A');

    // Update Assumptions & Decision Engine Panel
    const recBadge = document.getElementById('hypoRecommendationBadge');
    if (recBadge) {
      recBadge.innerText = `Recommended: ${assumptions.recommendedTestName}`;
      recBadge.className = 'badge badge-sig';
    }

    const decisionTextEl = document.getElementById('hypoDecisionText');
    if (decisionTextEl) {
      let decisionHtml = `<strong>Recommended Test:</strong> ${assumptions.recommendedTestName}<br>`;
      decisionHtml += `<span>${assumptions.rationale}</span>`;
      if (selectedTest !== 'auto') {
        if (selectedTest === assumptions.recommendedTest) {
          decisionHtml += `<div style="margin-top: 0.35rem; color: var(--emerald-primary); font-weight: 600;">✓ Manual test selection (${res.testName}) perfectly matches the statistical recommendation.</div>`;
        } else {
          decisionHtml += `<div style="margin-top: 0.35rem; color: var(--amber-primary); font-weight: 600;">⚠️ Advisory Note: Executed test (${res.testName}) was manually selected, differing from the assumption-recommended test (${assumptions.recommendedTestName}).</div>`;
        }
      } else {
        decisionHtml += `<div style="margin-top: 0.35rem; color: var(--cyan-primary); font-weight: 600;">⚡ Automatically executed: ${res.testName}</div>`;
      }
      decisionTextEl.innerHTML = decisionHtml;
    }

    // Update Normality Diagnostics
    const normBadge = document.getElementById('hypoNormalityBadge');
    const normDetails = document.getElementById('hypoNormalityDetails');
    if (isPaired) {
      if (normBadge) {
        normBadge.className = `badge ${assumptions.normality.isNormal ? 'badge-sig' : 'badge-ns'}`;
        normBadge.innerText = assumptions.normality.isNormal ? 'Normal Differences (Parametric)' : 'Non-Normal Differences (Non-Parametric)';
      }
      if (normDetails) {
        normDetails.innerHTML = `
          <strong>Paired Differences (&Delta; = Post &minus; Pre, n = ${assumptions.n}):</strong><br>
          Mean &Delta; = ${assumptions.statsDiff.mean.toFixed(2)}, SD = ${assumptions.statsDiff.sd.toFixed(2)}<br>
          Skewness = ${assumptions.statsDiff.skewness.toFixed(2)} (${assumptions.statsDiff.skewnessInterpretation})<br>
          Excess Kurtosis = ${assumptions.statsDiff.kurtosis.toFixed(2)} (${assumptions.statsDiff.kurtosisInterpretation})<br>
          Jarque-Bera Test: JB = ${assumptions.normality.statistic.toFixed(2)}, p = ${Exporter.formatP(assumptions.normality.pValue)}<br>
          <span style="font-weight: 600; color: ${assumptions.normality.isNormal ? 'var(--emerald-primary)' : 'var(--amber-primary)'};">${assumptions.normality.interpretation}</span>
        `;
      }
    } else {
      if (normBadge) {
        normBadge.className = `badge ${assumptions.normality.isParametric ? 'badge-sig' : 'badge-ns'}`;
        normBadge.innerText = assumptions.normality.isParametric ? 'Both Cohorts Normal (Parametric)' : 'Normality Violated (Non-Parametric)';
      }
      if (normDetails) {
        normDetails.innerHTML = `
          <strong>${nameA} (n = ${assumptions.statsA.n}):</strong> Skew = ${assumptions.statsA.skewness.toFixed(2)}, Kurt = ${assumptions.statsA.kurtosis.toFixed(2)}, JB = ${assumptions.normality.jbA.statistic.toFixed(2)} (${Exporter.formatP(assumptions.normality.jbA.pValue)}) [${assumptions.normality.normA ? 'Normal' : 'Skewed'}]<br>
          <strong>${nameB} (n = ${assumptions.statsB.n}):</strong> Skew = ${assumptions.statsB.skewness.toFixed(2)}, Kurt = ${assumptions.statsB.kurtosis.toFixed(2)}, JB = ${assumptions.normality.jbB.statistic.toFixed(2)} (${Exporter.formatP(assumptions.normality.jbB.pValue)}) [${assumptions.normality.normB ? 'Normal' : 'Skewed'}]<br>
          <span style="font-weight: 600; color: ${assumptions.normality.isParametric ? 'var(--emerald-primary)' : 'var(--amber-primary)'};">${assumptions.normality.interpretation}</span>
        `;
      }
    }

    // Update Variance Equality Diagnostics
    const varBadge = document.getElementById('hypoVarianceBadge');
    const varDetails = document.getElementById('hypoVarianceDetails');
    if (isPaired) {
      if (varBadge) {
        varBadge.className = 'badge badge-neutral';
        varBadge.innerText = 'N/A (Paired Design)';
      }
      if (varDetails) {
        varDetails.innerHTML = `
          <strong>Paired Repeated Measures:</strong><br>
          Between-cohort homoscedasticity is not required for paired analysis because the evaluation is performed on within-subject difference scores (&Delta;<sub>i</sub> = Post<sub>i</sub> &minus; Pre<sub>i</sub>), removing inter-subject variance. Sphericity is naturally satisfied with 2 repeated measures.
        `;
      }
    } else {
      const eq = assumptions.varianceEquality.equalVariance;
      if (varBadge) {
        varBadge.className = `badge ${eq ? 'badge-sig' : 'badge-ns'}`;
        varBadge.innerText = eq ? 'Equal Variances (Homoscedastic)' : 'Unequal Variances (Heteroscedastic)';
      }
      if (varDetails) {
        varDetails.innerHTML = `
          <strong>F-Test of Equal Variances:</strong><br>
          ${nameA} Variance s₁² = ${assumptions.statsA.variance.toFixed(2)} | ${nameB} Variance s₂² = ${assumptions.statsB.variance.toFixed(2)}<br>
          Variance Ratio F(${assumptions.varianceEquality.df1}, ${assumptions.varianceEquality.df2}) = ${assumptions.varianceEquality.fStat.toFixed(2)}, p = ${Exporter.formatP(assumptions.varianceEquality.pValue)}<br>
          <span style="font-weight: 600; color: ${eq ? 'var(--emerald-primary)' : 'var(--amber-primary)'};">${assumptions.varianceEquality.interpretation}</span>
        `;
      }
    }

    // Update Pipeline Flow Badges
    const flowDesign = document.getElementById('hypoFlowDesign');
    const flowNorm = document.getElementById('hypoFlowNorm');
    const flowVar = document.getElementById('hypoFlowVar');
    const flowTest = document.getElementById('hypoFlowTest');
    if (flowDesign) flowDesign.innerText = isPaired ? 'Paired Samples' : 'Independent Samples';
    if (flowNorm) flowNorm.innerText = (isPaired ? assumptions.normality.isNormal : assumptions.normality.isParametric) ? 'Parametric (Normal)' : 'Non-Parametric';
    if (flowVar) flowVar.innerText = isPaired ? 'Within-Subject' : (assumptions.varianceEquality.equalVariance ? 'Equal Variance' : 'Unequal Variance');
    if (flowTest) {
      flowTest.innerText = res.testName;
      flowTest.className = `badge ${res.isSignificant ? 'badge-sig' : 'badge-neutral'}`;
    }

    // Comprehensive Clinical / Publication APA Narrative
    let report = `A ${res.testName} was conducted to compare ${nameA} and ${nameB}.\n\n`;
    report += `Diagnostic Assumption Testing:\n`;
    if (isPaired) {
      report += `• Sample Design: Paired / repeated measures (n = ${assumptions.n} paired pairs).\n`;
      report += `• Normality of Within-Subject Differences: Jarque-Bera JB = ${assumptions.normality.statistic.toFixed(2)}, ${Exporter.formatP(assumptions.normality.pValue)} (Skewness = ${assumptions.normality.skewness.toFixed(2)}, Kurtosis = ${assumptions.normality.kurtosis.toFixed(2)}). The difference distribution was determined to be ${assumptions.normality.isNormal ? 'normally distributed' : 'non-normally distributed'}.\n`;
      report += `• Homoscedasticity: Not applicable for paired design (within-subject differencing removes inter-subject variance).\n`;
    } else {
      report += `• Sample Design: Independent two-cohort comparison (${nameA}: n = ${assumptions.statsA.n}; ${nameB}: n = ${assumptions.statsB.n}).\n`;
      report += `• Normality Assessment: ${nameA} (JB = ${assumptions.normality.jbA.statistic.toFixed(2)}, ${Exporter.formatP(assumptions.normality.jbA.pValue)}, Skew = ${assumptions.statsA.skewness.toFixed(2)}); ${nameB} (JB = ${assumptions.normality.jbB.statistic.toFixed(2)}, ${Exporter.formatP(assumptions.normality.jbB.pValue)}, Skew = ${assumptions.statsB.skewness.toFixed(2)}). Distribution: ${assumptions.normality.isParametric ? 'Parametric (Normal)' : 'Non-Parametric (Skewed/Deviated)'}.\n`;
      report += `• Equality of Variances: F-test F(${assumptions.varianceEquality.df1}, ${assumptions.varianceEquality.df2}) = ${assumptions.varianceEquality.fStat.toFixed(2)}, ${Exporter.formatP(assumptions.varianceEquality.pValue)}, confirming ${assumptions.varianceEquality.equalVariance ? 'equal variances (homoscedasticity)' : 'unequal variances (heteroscedasticity)'}.\n`;
    }
    report += `• Decision Rationale: ${assumptions.rationale}\n\n`;
    report += `Inferential Test Results:\n`;
    const statLabel = res.testName.includes('Mann-Whitney') ? 'U' : (res.testName.includes('Wilcoxon') ? 'W' : 't');
    const dfLabel = res.df !== undefined ? `(${res.df.toFixed(1)})` : '';
    const effectLabel = res.cohensD !== undefined ? `Cohen's d = ${res.cohensD.toFixed(2)}` : `Rank-Biserial r = ${(res.rankBiserial || 0).toFixed(2)}`;
    report += `There was a ${res.isSignificant ? 'statistically significant' : 'non-significant'} difference between ${nameA} (M = ${res.groupA.mean.toFixed(2)}, SD = ${res.groupA.sd.toFixed(2)}) and ${nameB} (M = ${res.groupB.mean.toFixed(2)}, SD = ${res.groupB.sd.toFixed(2)}): ${statLabel}${dfLabel} = ${(res.statistic !== undefined ? res.statistic : res.zScore || 0).toFixed(2)}, ${Exporter.formatP(res.pValue)}, ${effectLabel}.`;
    if (res.ci95) {
      report += ` 95% Confidence Interval: [${res.ci95[0].toFixed(2)}, ${res.ci95[1].toFixed(2)}].`;
    }

    document.getElementById('hypoReportText').innerText = report;

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
    this.results = this.results || {}; this.results.hypothesis = Object.assign(res, { assumptions, nameA, nameB, isPaired });
  }

  initAnovaGroups() {
    const container = document.getElementById('anovaGroupsContainer');
    if (!container) return;

    container.addEventListener('input', (e) => {
      if (e.target.classList.contains('anova-group-data')) {
        const card = e.target.closest('.anova-group-card');
        const countBadge = card?.querySelector('.anova-group-count');
        if (countBadge) {
          const count = DataParser.parseSeries(e.target.value).length;
          countBadge.innerText = `n = ${count}`;
        }
      }
    });

    container.addEventListener('click', (e) => {
      const btn = e.target.closest('.btn-remove-anova-group');
      if (btn) {
        const card = btn.closest('.anova-group-card');
        if (card) this.removeAnovaGroup(card);
      }
    });

    this.updateAnovaRemoveButtons();
    this.updateAnovaDesignUI();
  }

  addAnovaGroup(name = '', dataStr = '') {
    const container = document.getElementById('anovaGroupsContainer');
    if (!container) return;

    const currentCards = container.querySelectorAll('.anova-group-card');
    const k = currentCards.length + 1;
    const cohortName = name.trim() || `Cohort ${k}`;
    const count = dataStr ? DataParser.parseSeries(dataStr).length : 0;

    const card = document.createElement('div');
    card.className = 'anova-group-card';
    card.style.cssText = 'background: var(--bg-surface-elevated); border: 1px solid var(--border-subtle); border-radius: 8px; padding: 0.75rem 0.85rem;';
    card.innerHTML = `
      <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 0.4rem;">
        <div style="display: flex; align-items: center; gap: 0.5rem; flex: 1;">
          <span class="anova-group-index" style="font-weight: 700; color: var(--cyan-primary); font-size: 0.82rem; min-width: 60px;">Cohort ${k}:</span>
          <input type="text" class="form-control anova-group-name" value="${cohortName}" style="font-size: 0.82rem; padding: 0.25rem 0.5rem; height: 28px; font-weight: 600;" placeholder="Cohort Name">
        </div>
        <div style="display: flex; align-items: center; gap: 0.5rem;">
          <span class="badge badge-neutral anova-group-count" style="font-size: 0.72rem;">n = ${count}</span>
          <button type="button" class="btn btn-secondary btn-sm btn-remove-anova-group" style="padding: 0.15rem 0.45rem; font-size: 0.72rem; height: 26px; color: var(--rose-primary); border-color: rgba(244, 63, 94, 0.3);" title="Remove this cohort">✕</button>
        </div>
      </div>
      <textarea class="form-control anova-group-data" rows="2" placeholder="Comma or newline-separated values...">${dataStr}</textarea>
    `;

    container.appendChild(card);
    this.updateAnovaRemoveButtons();
    this.runAnova();
  }

  removeAnovaGroup(cardEl) {
    if (!cardEl) return;
    const container = document.getElementById('anovaGroupsContainer');
    const cards = container?.querySelectorAll('.anova-group-card');
    if (!cards || cards.length <= 2) {
      alert('Multi-group comparison requires at least 2 cohorts.');
      return;
    }

    cardEl.remove();

    // Re-index cohort headers
    const remaining = container.querySelectorAll('.anova-group-card');
    remaining.forEach((c, idx) => {
      const idxSpan = c.querySelector('.anova-group-index');
      if (idxSpan) idxSpan.innerText = `Cohort ${idx + 1}:`;
    });

    this.updateAnovaRemoveButtons();
    this.runAnova();
  }

  updateAnovaRemoveButtons() {
    const cards = document.querySelectorAll('#anovaGroupsContainer .anova-group-card');
    const canRemove = cards.length > 2;
    cards.forEach(card => {
      const btn = card.querySelector('.btn-remove-anova-group');
      if (btn) btn.disabled = !canRemove;
    });
  }

  updateAnovaDesignUI() {
    const isPaired = document.getElementById('anovaIsPaired')?.checked || false;
    const modeBadge = document.getElementById('anovaDesignModeBadge');
    const step1 = document.getElementById('anovaPipelineStep1');

    if (modeBadge) {
      if (isPaired) {
        modeBadge.className = 'badge badge-sig';
        modeBadge.innerText = 'Repeated Measures (Matched Within-Subjects)';
      } else {
        modeBadge.className = 'badge badge-neutral';
        modeBadge.innerText = 'Independent Cohorts (Between-Subjects)';
      }
    }

    if (step1) {
      step1.className = isPaired ? 'badge badge-sig' : 'badge badge-neutral';
      step1.innerText = `1. Design: ${isPaired ? 'Paired (RM)' : 'Independent'}`;
    }
  }

  loadAnovaPreset(presetKey) {
    const presets = {
      sample3: {
        isPaired: false,
        testType: 'auto',
        groups: [
          { name: 'Conservative', data: '7.2, 6.8, 7.5, 6.9, 8.1, 7.0, 7.4, 6.5, 7.9, 7.1' },
          { name: 'Orthotic Helmet', data: '4.1, 3.8, 4.5, 3.9, 4.8, 3.6, 4.2, 3.5, 4.0, 3.7' },
          { name: 'Endoscopic Strip', data: '2.5, 2.8, 2.2, 2.6, 3.1, 2.4, 2.9, 2.1, 2.7, 2.3' }
        ]
      },
      sampleWelch: {
        isPaired: false,
        testType: 'auto',
        groups: [
          { name: 'Cohort A (Small Var)', data: '10.1, 10.3, 10.0, 10.2, 9.9, 10.4, 10.1, 9.8, 10.2, 10.0' },
          { name: 'Cohort B (Mod Var)', data: '12.4, 11.1, 13.5, 10.8, 14.2, 11.9, 13.0, 12.1, 11.5, 13.8' },
          { name: 'Cohort C (Large Var)', data: '15.2, 8.5, 19.4, 11.1, 22.0, 14.3, 7.8, 18.6, 12.0, 20.5' },
          { name: 'Cohort D (High Var)', data: '18.0, 29.5, 9.2, 35.1, 14.8, 27.2, 8.1, 31.4, 12.5, 25.8' }
        ]
      },
      sampleRM: {
        isPaired: true,
        testType: 'auto',
        groups: [
          { name: 'Baseline (T0)', data: '22.4, 25.1, 19.8, 27.3, 23.5, 26.2, 21.9, 24.8, 20.5, 23.9' },
          { name: 'Week 2 (T1)', data: '19.1, 22.0, 17.5, 24.2, 20.8, 23.1, 18.9, 21.5, 18.0, 20.7' },
          { name: 'Week 6 (T2)', data: '15.3, 18.4, 14.1, 20.5, 17.0, 19.2, 15.6, 17.9, 14.8, 17.2' },
          { name: 'Month 3 (T3)', data: '12.1, 14.8, 11.2, 16.9, 13.5, 15.4, 12.4, 14.1, 11.5, 13.8' }
        ]
      },
      sampleSkew: {
        isPaired: false,
        testType: 'auto',
        groups: [
          { name: 'Standard Care', data: '1.2, 1.4, 1.1, 1.3, 1.5, 1.2, 1.4, 8.5, 12.3, 19.8' },
          { name: 'Modified Protocol', data: '2.1, 2.3, 2.0, 2.4, 2.2, 2.5, 15.1, 22.4, 28.0, 35.2' },
          { name: 'Novel Intervention', data: '5.5, 5.8, 5.2, 5.9, 6.1, 5.4, 32.0, 45.6, 58.2, 72.1' }
        ]
      }
    };

    const preset = presets[presetKey] || presets.sample3;
    const container = document.getElementById('anovaGroupsContainer');
    if (!container) return;

    container.innerHTML = '';
    preset.groups.forEach((g, idx) => {
      const k = idx + 1;
      const count = DataParser.parseSeries(g.data).length;
      const card = document.createElement('div');
      card.className = 'anova-group-card';
      card.style.cssText = 'background: var(--bg-surface-elevated); border: 1px solid var(--border-subtle); border-radius: 8px; padding: 0.75rem 0.85rem;';
      card.innerHTML = `
        <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 0.4rem;">
          <div style="display: flex; align-items: center; gap: 0.5rem; flex: 1;">
            <span class="anova-group-index" style="font-weight: 700; color: var(--cyan-primary); font-size: 0.82rem; min-width: 60px;">Cohort ${k}:</span>
            <input type="text" class="form-control anova-group-name" value="${g.name}" style="font-size: 0.82rem; padding: 0.25rem 0.5rem; height: 28px; font-weight: 600;" placeholder="Cohort Name">
          </div>
          <div style="display: flex; align-items: center; gap: 0.5rem;">
            <span class="badge badge-neutral anova-group-count" style="font-size: 0.72rem;">n = ${count}</span>
            <button type="button" class="btn btn-secondary btn-sm btn-remove-anova-group" style="padding: 0.15rem 0.45rem; font-size: 0.72rem; height: 26px; color: var(--rose-primary); border-color: rgba(244, 63, 94, 0.3);" title="Remove this cohort">✕</button>
          </div>
        </div>
        <textarea class="form-control anova-group-data" rows="2" placeholder="Comma or newline-separated values...">${g.data}</textarea>
      `;
      container.appendChild(card);
    });

    const pairedCb = document.getElementById('anovaIsPaired');
    if (pairedCb) pairedCb.checked = !!preset.isPaired;

    const testSelect = document.getElementById('anovaTestType');
    if (testSelect) testSelect.value = preset.testType || 'auto';

    this.updateAnovaRemoveButtons();
    this.updateAnovaDesignUI();
    this.runAnova();
  }

  getAnovaGroups() {
    const container = document.getElementById('anovaGroupsContainer');
    if (!container) return [];

    const cards = container.querySelectorAll('.anova-group-card');
    const groups = [];

    cards.forEach((card, idx) => {
      const nameInput = card.querySelector('.anova-group-name');
      const dataArea = card.querySelector('.anova-group-data');
      const countBadge = card.querySelector('.anova-group-count');

      const name = nameInput?.value.trim() || `Cohort ${idx + 1}`;
      const rawData = dataArea?.value || '';
      const parsed = DataParser.parseSeries(rawData);

      if (countBadge) {
        countBadge.innerText = `n = ${parsed.length}`;
      }
      groups.push({ name, data: parsed });
    });

    return groups;
  }

  runAnova() {
    const groups = this.getAnovaGroups();
    const isPaired = document.getElementById('anovaIsPaired')?.checked || false;
    const testType = document.getElementById('anovaTestType')?.value || 'auto';

    this.updateAnovaDesignUI();

    if (!groups || groups.length < 2) {
      alert('Please provide at least 2 cohorts for analysis.');
      return;
    }

    // Validate cohort data
    for (let i = 0; i < groups.length; i++) {
      if (!groups[i].data || groups[i].data.length === 0) {
        const msg = `Cohort "${groups[i].name}" has no valid numerical data.`;
        document.getElementById('anovaReportText').innerText = msg;
        return;
      }
    }

    const res = Anova.test(groups, testType, isPaired);
    if (res.error) {
      document.getElementById('anovaReportText').innerText = `Analysis Error: ${res.error}`;
      return;
    }

    // Update Metric Cards
    const fLabel = document.getElementById('anovaFLabel');
    const fVal = document.getElementById('anovaF');
    const pVal = document.getElementById('anovaP');
    const pBadge = document.getElementById('anovaPValBadge');
    const etaLabel = document.getElementById('anovaEtaLabel');
    const etaVal = document.getElementById('anovaEta');
    const omegaLabel = document.getElementById('anovaOmegaLabel');
    const omegaVal = document.getElementById('anovaOmega');

    // 1. Test Statistic
    if (fLabel) {
      if (res.testKey === 'kruskal') fLabel.innerText = 'Kruskal-Wallis (H)';
      else if (res.testKey === 'friedman') fLabel.innerText = 'Friedman (Q / χ²ᵣ)';
      else if (res.testKey === 'welch') fLabel.innerText = "Welch's F-Test";
      else if (res.testKey === 'rm_anova') fLabel.innerText = 'RM-ANOVA (F)';
      else fLabel.innerText = 'Fisher ANOVA (F)';
    }
    if (fVal) {
      const stat = res.statistic !== undefined ? res.statistic : (res.fStatistic || res.hStatistic || res.qStatistic || 0);
      fVal.innerText = isFinite(stat) ? stat.toFixed(2) : '--';
    }

    // 2. p-value & Significance badge
    if (pVal) pVal.innerText = Exporter.formatP(res.pValue);
    if (pBadge) {
      pBadge.className = res.isSignificant ? 'badge badge-sig' : 'badge badge-ns';
      pBadge.innerText = res.isSignificant ? 'Significant (p < .05)' : 'Not Significant (ns)';
    }

    // 3. Effect Size
    if (etaLabel) {
      if (res.testKey === 'kruskal') etaLabel.innerText = 'Epsilon-Squared (ε²)';
      else if (res.testKey === 'friedman') etaLabel.innerText = "Kendall's W";
      else if (res.testKey === 'welch') etaLabel.innerText = 'Estimated ω²';
      else etaLabel.innerText = 'Partial η² / ω²';
    }
    if (etaVal) {
      if (res.testKey === 'kruskal' && res.effectSize?.epsilonSquared !== undefined) {
        etaVal.innerText = res.effectSize.epsilonSquared.toFixed(3);
      } else if (res.testKey === 'friedman' && res.effectSize?.kendallsW !== undefined) {
        etaVal.innerText = res.effectSize.kendallsW.toFixed(3);
      } else if (res.etaSquared !== undefined) {
        etaVal.innerText = res.etaSquared.toFixed(3);
      } else if (res.omegaSquared !== undefined) {
        etaVal.innerText = res.omegaSquared.toFixed(3);
      } else {
        etaVal.innerText = '--';
      }
    }

    // 4. Design & Degrees of Freedom
    if (omegaLabel) {
      omegaLabel.innerText = res.isPaired ? 'Paired Design & df' : 'Independent Design & df';
    }
    if (omegaVal) {
      if (res.dfBetween !== undefined && res.dfWithin !== undefined) {
        omegaVal.innerText = `df: (${res.dfBetween}, ${typeof res.dfWithin === 'number' ? res.dfWithin.toFixed(1) : res.dfWithin})`;
      } else if (res.df !== undefined) {
        omegaVal.innerText = `df = ${res.df}`;
      } else {
        omegaVal.innerText = `k = ${res.k}`;
      }
    }

    // 5. Assumptions Diagnostic Card
    const asm = res.assumptions;
    if (asm) {
      const recBadge = document.getElementById('anovaRecommendationBadge');
      if (recBadge) {
        recBadge.innerText = `Recommended: ${asm.recommendedTestName}`;
        recBadge.className = 'badge badge-sig';
      }

      const decisionBanner = document.getElementById('anovaDecisionBanner');
      const decisionText = document.getElementById('anovaDecisionText');
      if (decisionText) {
        let bannerHtml = '';
        if (res.userOverride) {
          bannerHtml = `<span class="badge badge-warn" style="font-size: 0.72rem; margin-bottom: 0.35rem; display: inline-block;">Manual Selection</span><br>` +
            `User opted to execute <strong>${res.testName}</strong>. Based on data diagnostics, the statistically optimal test is <strong>${asm.recommendedTestName}</strong>.<br>` +
            `<span style="color: var(--text-dim); font-size: 0.82rem;">${asm.rationale}</span>`;
          if (decisionBanner) decisionBanner.style.borderLeftColor = 'var(--gold-primary)';
        } else {
          bannerHtml = `<span class="badge badge-sig" style="font-size: 0.72rem; margin-bottom: 0.35rem; display: inline-block;">Automated Recommendation Executed</span><br>` +
            `Executed <strong>${res.testName}</strong>.<br>` +
            `<span style="color: var(--text-main);">${asm.rationale}</span>`;
          if (decisionBanner) decisionBanner.style.borderLeftColor = 'var(--cyan-primary)';
        }
        decisionText.innerHTML = bannerHtml;
      }

      // Normality Breakdown
      const normBadge = document.getElementById('anovaNormalityBadge');
      const normDetails = document.getElementById('anovaNormalityDetails');
      if (normBadge) {
        normBadge.className = asm.isNormal ? 'badge badge-sig' : 'badge badge-warn';
        normBadge.innerText = asm.isNormal ? 'Parametric (All Normal)' : 'Non-Parametric (Skewed)';
      }
      if (normDetails && asm.normalityTests) {
        normDetails.innerHTML = asm.normalityTests.map(n => 
          `• <strong>${n.group}</strong>: n=${n.n}, JB=${n.jbStat.toFixed(2)}, p=${Exporter.formatP(n.pValue)} (${n.isNormal ? '<span style="color: var(--emerald-primary);">Normal</span>' : '<span style="color: var(--rose-primary);">Skewed, p < .05</span>'})`
        ).join('<br>');
      }

      // Variance Homogeneity / Sphericity Breakdown
      const varBadge = document.getElementById('anovaVarianceBadge');
      const varDetails = document.getElementById('anovaVarianceDetails');
      if (varBadge) {
        if (isPaired) {
          varBadge.className = 'badge badge-sig';
          varBadge.innerText = asm.sphericity ? `Sphericity ε̂ = ${asm.sphericity.epsilon.toFixed(3)}` : 'Repeated Measures';
        } else {
          varBadge.className = asm.isHomoscedastic ? 'badge badge-sig' : 'badge badge-warn';
          varBadge.innerText = asm.isHomoscedastic ? 'Equal Variances' : 'Unequal Variances';
        }
      }
      if (varDetails) {
        if (isPaired) {
          if (asm.sphericity) {
            varDetails.innerHTML = `Greenhouse-Geisser correction factor: <strong>ε̂ = ${asm.sphericity.epsilon.toFixed(3)}</strong>.<br>` +
              (asm.sphericity.isSpherical 
                ? `<span style="color: var(--emerald-primary);">Sphericity assumption reasonably met (ε̂ ≈ 1.0).</span>` 
                : `<span style="color: var(--gold-primary);">Sphericity violated (ε̂ < 0.75). Degrees of freedom adjusted via Greenhouse-Geisser.</span>`);
          } else {
            varDetails.innerText = 'Matched repeated measures design: within-subject correlation structure preserved.';
          }
        } else if (asm.leveneTest) {
          const lev = asm.leveneTest;
          varDetails.innerHTML = `Brown-Forsythe Levene's Test: <strong>F(${lev.df1}, ${lev.df2}) = ${lev.fStat.toFixed(2)}, p = ${Exporter.formatP(lev.pValue)}</strong>.<br>` +
            (asm.isHomoscedastic 
              ? `<span style="color: var(--emerald-primary);">Homoscedasticity confirmed (p ≥ .05). Residual variances across cohorts are equal.</span>`
              : `<span style="color: var(--rose-primary);">Heteroscedasticity detected (p < .05). Residual variances differ significantly; Welch's robust F recommended.</span>`);
        }
      }

      // Pipeline step badges
      const pStep1 = document.getElementById('anovaPipelineStep1');
      const pStep2 = document.getElementById('anovaPipelineStep2');
      const pStep3 = document.getElementById('anovaPipelineStep3');
      const pStepFinal = document.getElementById('anovaPipelineStepFinal');
      if (pStep1) {
        pStep1.className = isPaired ? 'badge badge-sig' : 'badge badge-neutral';
        pStep1.innerText = `1. Design: ${isPaired ? 'Paired (RM)' : 'Independent'}`;
      }
      if (pStep2) {
        pStep2.className = asm.isNormal ? 'badge badge-sig' : 'badge badge-warn';
        pStep2.innerText = `2. Normality: ${asm.isNormal ? 'Parametric' : 'Non-Parametric'}`;
      }
      if (pStep3) {
        if (isPaired) {
          pStep3.className = 'badge badge-sig';
          pStep3.innerText = `3. Sphericity: ε̂ = ${asm.sphericity ? asm.sphericity.epsilon.toFixed(2) : 'N/A'}`;
        } else {
          pStep3.className = asm.isHomoscedastic ? 'badge badge-sig' : 'badge badge-warn';
          pStep3.innerText = `3. Variances: ${asm.isHomoscedastic ? 'Equal' : 'Unequal'}`;
        }
      }
      if (pStepFinal) {
        pStepFinal.className = 'badge badge-sig';
        pStepFinal.innerText = `Selected: ${res.testName}`;
      }
    }

    // 6. Post-Hoc Pairwise Table
    const postHocTitle = document.getElementById('anovaPostHocTitle');
    const postHocSub = document.getElementById('anovaPostHocSubtitle');
    const statCol = document.getElementById('anovaPostHocStatCol');
    if (postHocTitle) postHocTitle.innerText = `🔬 Post-Hoc Pairwise Contrasts (${res.postHocMethod || 'Pairwise'})`;
    if (postHocSub) postHocSub.innerText = `Contrasts between individual cohorts with family-wise error rate control (${res.postHocMethod || 'Contrasts'})`;
    if (statCol) {
      if (res.testKey === 'kruskal') statCol.innerText = "Dunn's z-Stat";
      else if (res.testKey === 'friedman') statCol.innerText = 'Wilcoxon W';
      else if (res.testKey === 'welch') statCol.innerText = 'Games-Howell t';
      else if (res.testKey === 'rm_anova') statCol.innerText = 'Paired t-Stat';
      else statCol.innerText = 'Tukey q (t)';
    }

    const tbody = document.getElementById('anovaPostHocBody');
    if (tbody) {
      if (res.pairwise && res.pairwise.length > 0) {
        tbody.innerHTML = res.pairwise.map(p => {
          const diffVal = p.meanDiff !== undefined ? p.meanDiff : (p.diff !== undefined ? p.diff : 0);
          const seVal = p.seDiff !== undefined ? p.seDiff.toFixed(2) : (p.se !== undefined ? p.se.toFixed(2) : '--');
          const statText = p.qStatistic !== undefined ? `q = ${p.qStatistic.toFixed(2)} (t = ${p.tStatistic.toFixed(2)})` :
                           p.tStatistic !== undefined ? `t = ${p.tStatistic.toFixed(2)}` :
                           p.zStatistic !== undefined ? `z = ${p.zStatistic.toFixed(2)}` :
                           p.wStatistic !== undefined ? `W = ${p.wStatistic.toFixed(2)}` : '--';
          const ciText = p.ci95 ? `[${p.ci95[0].toFixed(2)}, ${p.ci95[1].toFixed(2)}]` : '--';
          const esText = p.cohensD !== undefined ? `d = ${p.cohensD.toFixed(2)}` :
                         p.r !== undefined ? `r = ${p.r.toFixed(2)}` : '--';
          return `
            <tr>
              <td style="font-weight: 600; color: var(--text-main);">${p.comparison}</td>
              <td>${diffVal >= 0 ? '+' : ''}${diffVal.toFixed(2)}</td>
              <td>${seVal}</td>
              <td>${statText}</td>
              <td style="font-weight: 600; color: ${p.isSignificant ? 'var(--cyan-primary)' : 'var(--text-muted)'};">${Exporter.formatP(p.pValue)}</td>
              <td>${ciText}</td>
              <td>${esText}</td>
              <td>
                <span class="badge ${p.isSignificant ? 'badge-sig' : 'badge-ns'}">
                  ${p.isSignificant ? 'Significant (p < .05)' : 'Not Significant (ns)'}
                </span>
              </td>
            </tr>
          `;
        }).join('');
      } else {
        tbody.innerHTML = `<tr><td colspan="8" style="text-align: center; color: var(--text-dim);">No pairwise contrasts calculated (omnibus effect not significant or single cohort).</td></tr>`;
      }
    }

    // 7. Clinical / Academic APA Summary Narrative
    let report = '';
    if (res.testKey === 'one_way') {
      report = `A one-way between-subjects ANOVA was conducted across ${res.k} cohorts (N = ${res.totalN}). `;
      report += `There was a ${res.isSignificant ? 'statistically significant' : 'non-significant'} omnibus effect: F(${res.dfBetween}, ${res.dfWithin}) = ${res.fStatistic.toFixed(2)}, ${Exporter.formatP(res.pValue)}, η² = ${res.etaSquared.toFixed(3)}, ω² = ${res.omegaSquared.toFixed(3)}.\n\n`;
    } else if (res.testKey === 'welch') {
      report = `A Welch's robust one-way ANOVA (adjusting for heteroscedasticity) was conducted across ${res.k} cohorts (N = ${res.totalN}). `;
      report += `There was a ${res.isSignificant ? 'statistically significant' : 'non-significant'} omnibus effect: Welch's F(${res.dfBetween}, ${res.dfWithin.toFixed(2)}) = ${res.fStatistic.toFixed(2)}, ${Exporter.formatP(res.pValue)}, estimated ω² = ${res.omegaSquared.toFixed(3)}.\n\n`;
    } else if (res.testKey === 'kruskal') {
      report = `A non-parametric Kruskal-Wallis H test was conducted across ${res.k} cohorts (N = ${res.totalN}). `;
      report += `There was a ${res.isSignificant ? 'statistically significant' : 'non-significant'} omnibus rank difference: H(${res.df}) = ${res.hStatistic.toFixed(2)}, ${Exporter.formatP(res.pValue)}, ε² = ${res.effectSize.epsilonSquared.toFixed(3)}.\n\n`;
    } else if (res.testKey === 'rm_anova') {
      report = `A one-way repeated measures ANOVA was conducted across ${res.k} conditions (N = ${res.nSubjects} subjects). `;
      report += `Greenhouse-Geisser sphericity correction: ε̂ = ${res.epsilon.toFixed(3)}. Omnibus effect: F(${res.dfTreatment.toFixed(2)}, ${res.dfError.toFixed(2)}) = ${res.fStatistic.toFixed(2)}, ${Exporter.formatP(res.pValue)}, partial η² = ${res.partialEtaSquared.toFixed(3)}.\n\n`;
    } else if (res.testKey === 'friedman') {
      report = `A non-parametric Friedman rank sum test was conducted across ${res.k} repeated conditions (N = ${res.n} subjects). `;
      report += `Omnibus rank difference: Q(${res.df}) = ${res.qStatistic.toFixed(2)}, ${Exporter.formatP(res.pValue)}, Kendall's W = ${res.effectSize.kendallsW.toFixed(3)}.\n\n`;
    }

    // Append post-hoc summary
    if (res.pairwise && res.pairwise.length > 0) {
      report += `Post-hoc contrasts (${res.postHocMethod}) revealed:\n`;
      res.pairwise.forEach(p => {
        const diffVal = p.meanDiff !== undefined ? p.meanDiff : (p.diff !== undefined ? p.diff : 0);
        const statVal = p.qStatistic !== undefined ? `q = ${p.qStatistic.toFixed(2)}` :
                        p.tStatistic !== undefined ? `t = ${p.tStatistic.toFixed(2)}` :
                        p.zStatistic !== undefined ? `z = ${p.zStatistic.toFixed(2)}` :
                        p.wStatistic !== undefined ? `W = ${p.wStatistic.toFixed(2)}` : '';
        if (p.isSignificant) {
          report += `• ${p.comparison}: Statistically significant difference (Δ = ${diffVal.toFixed(2)}, ${statVal ? statVal + ', ' : ''}${Exporter.formatP(p.pValue)}).\n`;
        } else {
          report += `• ${p.comparison}: No statistically significant difference (Δ = ${diffVal.toFixed(2)}, ${statVal ? statVal + ', ' : ''}${Exporter.formatP(p.pValue)}, ns).\n`;
        }
      });
    }

    // Append assumption decision rationale
    if (asm) {
      report += `\nMethodological Rationale: ${asm.rationale}`;
    }

    document.getElementById('anovaReportText').innerText = report;

    // 8. Dispersion Plot
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
    this.results = this.results || {};
    this.results.anova = res;
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
      exportData = Object.assign({}, res, {
        reportText: document.getElementById('anovaReportText')?.innerText
      });
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
        tConv: res.tConv || Teaching.tConvergence.getMetrics(4),
        overlap: res.overlap || Teaching.significanceOverlap.getMetrics(),
        power: res.power || Teaching.powerSimulation.getMetrics(),
        reportText: document.getElementById('teachingReportText')?.innerText
      };
    } else if (tabId === 'teaching-bayesian') {
      exportData = {
        bayes: res.bayes || Teaching.bayesianSimulation.getMetrics(),
        reportText: document.getElementById('teachingBayesReportText')?.innerText || document.getElementById('bayesPedagogyText')?.innerText
      };
    } else if (tabId === 'propensity') {
      exportData = {
        name: 'Clinical Observational Cohort',
        treatmentCol: res.treatmentCol,
        outcomeCol: res.outcomeCol,
        covariateCols: res.covariateCols,
        totalN: res.totalN,
        completeN: res.completeN,
        nMatchedPairs: res.nMatchedPairs,
        matchedTreatedN: res.matchedTreatedN,
        matchedControlN: res.matchedControlN,
        logisticRegression: res.logisticRegression,
        caliper: res.caliper,
        balance: res.balance,
        outcome: res.outcome,
        reportText: document.getElementById('psmReportText')?.innerText
      };
    }

    Exporter.exportTabToDocx(tabId, exportData);
  }

  // --- Teaching & Simulation Module ---

  initTeachingModule() {
    this.renderTeachingParams();
    this.runTeachingDistribution();
    this.updateCltUI(Teaching.clt.getSummary());
    this.runTConvergence(4);
    this.runTwoSampleOverlap();
    this.runPowerSimulation();
    this.runBayesianSimulation();
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
    this.runTConvergence();
  }

  runTConvergence(sampleSize = null) {
    const rangeEl = document.getElementById('tConvNRange');
    if (sampleSize !== null && rangeEl) {
      rangeEl.value = sampleSize;
    }
    const n = parseInt(rangeEl?.value) || 4;
    const valEl = document.getElementById('tConvNVal');
    if (valEl) {
      valEl.innerText = `n = ${n} (ν = ${n - 1})`;
    }

    const metrics = Teaching.tConvergence.getMetrics(n);

    // Update UI Cards
    const dfEl = document.getElementById('tConvDf');
    const dfSubEl = document.getElementById('tConvDfSub');
    if (dfEl) dfEl.innerText = metrics.df;
    if (dfSubEl) dfSubEl.innerText = `Sample size n = ${metrics.sampleSize}`;

    const critEl = document.getElementById('tConvCrit');
    const critSubEl = document.getElementById('tConvCritSub');
    if (critEl) critEl.innerText = metrics.tCrit.toFixed(3);
    if (critSubEl) critSubEl.innerText = `vs z = 1.960 (${metrics.critDiffPct >= 0 ? '+' : ''}${metrics.critDiffPct.toFixed(1)}%)`;

    const peakEl = document.getElementById('tConvPeak');
    const peakSubEl = document.getElementById('tConvPeakSub');
    if (peakEl) peakEl.innerText = metrics.tPeak.toFixed(4);
    if (peakSubEl) peakSubEl.innerText = `Normal: 0.3989 (${metrics.peakDiffPct.toFixed(1)}%)`;

    const tailEl = document.getElementById('tConvTailProb');
    const tailSubEl = document.getElementById('tConvTailProbSub');
    if (tailEl) {
      tailEl.innerText = `${(metrics.tailProb * 100).toFixed(1)}%`;
      tailEl.style.color = metrics.df < 30 ? '#ef4444' : '#22c55e';
    }
    if (tailSubEl) tailSubEl.innerText = metrics.df < 30 ? 'Normal: 5.0% (Type I Risk)' : 'Normal: 5.0% (Matched)';

    const maxDiffEl = document.getElementById('tConvMaxDiff');
    if (maxDiffEl) maxDiffEl.innerText = metrics.maxDiscrepancy.toFixed(4);

    const kurtEl = document.getElementById('tConvKurt');
    const kurtSubEl = document.getElementById('tConvKurtSub');
    if (kurtEl) {
      kurtEl.innerText = metrics.excessKurtosis === Infinity ? '∞' : metrics.excessKurtosis.toFixed(2);
    }
    if (kurtSubEl) {
      kurtSubEl.innerText = metrics.df <= 4 ? 'Fat-tailed (ν ≤ 4)' : `Excess Kurtosis: 6/(ν-4)`;
    }

    const pedagogyEl = document.getElementById('tConvPedagogyText');
    if (pedagogyEl) pedagogyEl.innerText = metrics.clinicalNote;

    const chartTitleEl = document.getElementById('tConvChartTitle');
    if (chartTitleEl) {
      chartTitleEl.innerText = `Student's t(ν = ${metrics.df}) Density Curve vs Standard Normal N(0, 1)`;
    }

    const showTailArea = document.getElementById('tConvShowTailArea')?.checked !== false;

    // Render Canvas
    if (this.engines['teachingTCanvas']) {
      Plots.renderTConvergence(this.engines['teachingTCanvas'], metrics, {
        showTailArea,
        title: `Student's t(ν = ${metrics.df}) vs Standard Normal N(0, 1)`
      });
    }

    // Cache
    if (!this.results.teaching) this.results.teaching = {};
    this.results.teaching.tConv = metrics;
  }

  animateTConvergence() {
    const btn = document.getElementById('tConvAnimateBtn');
    if (this.tConvAnimationTimer) {
      clearInterval(this.tConvAnimationTimer);
      this.tConvAnimationTimer = null;
      if (btn) btn.innerText = '▶ Animate Convergence';
      return;
    }

    const frames = [2, 3, 4, 5, 7, 10, 15, 20, 25, 31, 45, 61, 80, 100, 121, 150];
    let currentIndex = 0;
    const currentN = parseInt(document.getElementById('tConvNRange')?.value) || 4;
    const startIdx = frames.findIndex(f => f >= currentN);
    if (startIdx >= 0 && startIdx < frames.length - 1) currentIndex = startIdx;

    if (btn) btn.innerText = '⏸ Pause Animation';

    this.tConvAnimationTimer = setInterval(() => {
      currentIndex++;
      if (currentIndex >= frames.length) {
        clearInterval(this.tConvAnimationTimer);
        this.tConvAnimationTimer = null;
        if (btn) btn.innerText = '▶ Replay Animation';
        return;
      }
      this.runTConvergence(frames[currentIndex]);
    }, 450);
  }

  runTwoSampleOverlap(overrideParams = {}) {
    const delta = overrideParams.delta !== undefined
      ? overrideParams.delta
      : (parseFloat(document.getElementById('overlapDeltaRange')?.value) || 2.0);
    const alpha = overrideParams.alpha !== undefined
      ? overrideParams.alpha
      : (parseFloat(document.getElementById('overlapAlphaRange')?.value) || 0.05);

    const sd1 = overrideParams.sd1 !== undefined
      ? overrideParams.sd1
      : (parseFloat(document.getElementById('overlapSD1Range')?.value) || 2.5);
    const sem1 = overrideParams.sem1 !== undefined
      ? overrideParams.sem1
      : (parseFloat(document.getElementById('overlapSEM1Range')?.value) || 0.625);
    const n1 = overrideParams.n1 !== undefined
      ? overrideParams.n1
      : (parseInt(document.getElementById('overlapN1Range')?.value) || 16);

    const sd2 = overrideParams.sd2 !== undefined
      ? overrideParams.sd2
      : (parseFloat(document.getElementById('overlapSD2Range')?.value) || 2.5);
    const sem2 = overrideParams.sem2 !== undefined
      ? overrideParams.sem2
      : (parseFloat(document.getElementById('overlapSEM2Range')?.value) || 0.625);
    const n2 = overrideParams.n2 !== undefined
      ? overrideParams.n2
      : (parseInt(document.getElementById('overlapN2Range')?.value) || 16);

    const viewMode = this.currentOverlapMode || 'means';

    const metrics = Teaching.significanceOverlap.getMetrics({
      mean1: 10.0,
      delta,
      sd1,
      sd2,
      sem1,
      sem2,
      n1,
      n2,
      alpha,
      viewMode
    });

    // Update Slider Displays
    const deltaValEl = document.getElementById('overlapDeltaVal');
    if (deltaValEl) deltaValEl.innerText = metrics.delta.toFixed(2);

    const alphaValEl = document.getElementById('overlapAlphaVal');
    if (alphaValEl) alphaValEl.innerText = `α = ${metrics.alpha.toFixed(3)} (z = ${metrics.zCrit.toFixed(3)})`;

    // Group 1 Displays
    const g1SummaryEl = document.getElementById('overlapG1Summary');
    if (g1SummaryEl) g1SummaryEl.innerText = `SD₁ = ${metrics.sd1.toFixed(2)} | SEM₁ = ${metrics.sem1.toFixed(3)}`;

    const sd1ValEl = document.getElementById('overlapSD1Val');
    if (sd1ValEl) sd1ValEl.innerText = metrics.sd1.toFixed(2);

    const sem1ValEl = document.getElementById('overlapSEM1Val');
    if (sem1ValEl) sem1ValEl.innerText = metrics.sem1.toFixed(3);

    const n1ValEl = document.getElementById('overlapN1Val');
    if (n1ValEl) n1ValEl.innerText = `n₁ = ${metrics.n1}`;

    // Group 2 Displays
    const g2SummaryEl = document.getElementById('overlapG2Summary');
    if (g2SummaryEl) g2SummaryEl.innerText = `SD₂ = ${metrics.sd2.toFixed(2)} | SEM₂ = ${metrics.sem2.toFixed(3)}`;

    const sd2ValEl = document.getElementById('overlapSD2Val');
    if (sd2ValEl) sd2ValEl.innerText = metrics.sd2.toFixed(2);

    const sem2ValEl = document.getElementById('overlapSEM2Val');
    if (sem2ValEl) sem2ValEl.innerText = metrics.sem2.toFixed(3);

    const n2ValEl = document.getElementById('overlapN2Val');
    if (n2ValEl) n2ValEl.innerText = `n₂ = ${metrics.n2}`;

    // Update Metric Cards
    const statusValEl = document.getElementById('overlapStatusValue');
    const pValSubEl = document.getElementById('overlapPValueSub');
    if (statusValEl) {
      statusValEl.innerText = metrics.isSignificant ? 'SIGNIFICANT' : 'NOT SIGNIFICANT';
      statusValEl.style.color = metrics.isSignificant ? '#10b981' : '#ef4444';
    }
    if (pValSubEl) {
      pValSubEl.innerText = `p = ${metrics.pValue < 0.0001 ? '< 0.0001' : metrics.pValue.toFixed(4)} ${metrics.isSignificant ? '<' : '≥'} α = ${metrics.alpha.toFixed(3)}`;
    }

    const deltaDispEl = document.getElementById('overlapDeltaDisplay');
    const deltaCritSubEl = document.getElementById('overlapDeltaCritSub');
    if (deltaDispEl) deltaDispEl.innerText = `Δ = ${metrics.delta.toFixed(2)}`;
    if (deltaCritSubEl) deltaCritSubEl.innerText = `Δcrit = ${metrics.deltaCrit.toFixed(2)} (Boundary)`;

    const g1MetricVal = document.getElementById('overlapG1MetricValue');
    const g1MetricSub = document.getElementById('overlapG1MetricSub');
    if (g1MetricVal) g1MetricVal.innerText = `${metrics.sd1.toFixed(2)} | ${metrics.sem1.toFixed(3)}`;
    if (g1MetricSub) g1MetricSub.innerText = `n₁ = ${metrics.n1} subjects`;

    const g2MetricVal = document.getElementById('overlapG2MetricValue');
    const g2MetricSub = document.getElementById('overlapG2MetricSub');
    if (g2MetricVal) g2MetricVal.innerText = `${metrics.sd2.toFixed(2)} | ${metrics.sem2.toFixed(3)}`;
    if (g2MetricSub) g2MetricSub.innerText = `n₂ = ${metrics.n2} subjects`;

    const meansOVLEl = document.getElementById('overlapMeansOVLValue');
    const patientOVLSubEl = document.getElementById('overlapPatientOVLSub');
    if (meansOVLEl) {
      meansOVLEl.innerText = `${(metrics.meansOVL * 100).toFixed(1)}%`;
      meansOVLEl.style.color = metrics.isSignificant ? '#10b981' : '#ef4444';
    }
    if (patientOVLSubEl) patientOVLSubEl.innerText = `Patient Overlap: ${(metrics.patientOVL * 100).toFixed(1)}%`;

    const critValEl = document.getElementById('overlapCritValue');
    const alphaSubEl = document.getElementById('overlapAlphaSub');
    if (critValEl) critValEl.innerText = `t = ${metrics.tCrit.toFixed(3)}`;
    if (alphaSubEl) alphaSubEl.innerText = `SE_diff: ${metrics.seDiff.toFixed(3)} | df: ${metrics.df.toFixed(1)}`;

    // Update Pedagogical Text
    const pedaEl = document.getElementById('overlapPedagogyText');
    if (pedaEl) pedaEl.innerText = metrics.explanation;

    // Update Chart Title
    const titleEl = document.getElementById('overlapChartTitle');
    if (titleEl) {
      if (viewMode === 'patients') {
        titleEl.innerText = `Individual Patient Populations (SD): Biological Overlap vs Mean Significance`;
      } else if (viewMode === 'dual') {
        titleEl.innerText = `Dual Overlay: Patient Biological Spread (SD) vs Inferential Mean Precision (SEM)`;
      } else if (viewMode === 'null') {
        titleEl.innerText = `Null Hypothesis Difference Test: H₀ (Δ=0) vs Observed Separation`;
      } else {
        titleEl.innerText = `Sampling Distributions of Means (SEM): Overlap & Significance Boundary`;
      }
    }

    // Render Canvas
    if (this.engines['teachingOverlapCanvas']) {
      Plots.renderTwoSampleOverlap(this.engines['teachingOverlapCanvas'], metrics);
    }

    // Cache
    if (!this.results.teaching) this.results.teaching = {};
    this.results.teaching.overlap = metrics;
  }

  animateOverlapSeparation() {
    const btn = document.getElementById('overlapAnimateBtn');
    if (this.overlapAnimationTimer) {
      clearInterval(this.overlapAnimationTimer);
      this.overlapAnimationTimer = null;
      if (btn) btn.innerText = '▶ Animate Separation';
      return;
    }

    const deltas = [0.0, 0.4, 0.8, 1.2, 1.5, 1.8, 2.0, 2.3, 2.7, 3.2, 3.8, 4.5];
    let currentIndex = 0;
    const currentDelta = parseFloat(document.getElementById('overlapDeltaRange')?.value) || 0;
    const startIdx = deltas.findIndex(d => d >= currentDelta);
    if (startIdx >= 0 && startIdx < deltas.length - 1) currentIndex = startIdx;

    if (btn) btn.innerText = '⏸ Pause Animation';

    this.overlapAnimationTimer = setInterval(() => {
      currentIndex++;
      if (currentIndex >= deltas.length) {
        clearInterval(this.overlapAnimationTimer);
        this.overlapAnimationTimer = null;
        if (btn) btn.innerText = '▶ Replay Animation';
        return;
      }
      const nextDelta = deltas[currentIndex];
      const dRange = document.getElementById('overlapDeltaRange');
      if (dRange) dRange.value = nextDelta;
      this.runTwoSampleOverlap({ delta: nextDelta });
    }, 500);
  }

  // --- Teaching Simulation 5: Power Simulation Methods ---

  runPowerSimulation(overrideParams = {}, lastChanged = 'power') {
    const sd = overrideParams.sd !== undefined
      ? overrideParams.sd
      : (parseFloat(document.getElementById('powerSDRange')?.value) || 4.0);
    const sem = overrideParams.sem !== undefined
      ? overrideParams.sem
      : (parseFloat(document.getElementById('powerSEMRange')?.value) || 0.50);
    const n = overrideParams.n !== undefined
      ? overrideParams.n
      : (parseInt(document.getElementById('powerNRange')?.value) || 64);
    const power = overrideParams.power !== undefined
      ? overrideParams.power
      : (parseFloat(document.getElementById('powerPowerRange')?.value) || 0.80);
    const beta = overrideParams.beta !== undefined
      ? overrideParams.beta
      : (parseFloat(document.getElementById('powerBetaRange')?.value) || 0.20);
    const delta = overrideParams.delta !== undefined
      ? overrideParams.delta
      : (parseFloat(document.getElementById('powerDeltaRange')?.value) || 2.0);
    const alpha = overrideParams.alpha !== undefined
      ? overrideParams.alpha
      : (parseFloat(document.getElementById('powerAlphaRange')?.value) || 0.05);

    const viewMode = this.currentPowerMode || 'distributions';

    const metrics = Teaching.powerSimulation.getMetrics({
      sd,
      sem,
      n,
      power,
      beta,
      delta,
      alpha,
      viewMode,
      lastChanged
    });

    // Synchronize UI Slider Values without trigger loop
    const sdInput = document.getElementById('powerSDRange');
    const semInput = document.getElementById('powerSEMRange');
    const nInput = document.getElementById('powerNRange');
    const powerInput = document.getElementById('powerPowerRange');
    const betaInput = document.getElementById('powerBetaRange');
    const deltaInput = document.getElementById('powerDeltaRange');
    const alphaInput = document.getElementById('powerAlphaRange');

    if (sdInput && lastChanged !== 'sd') sdInput.value = metrics.sd.toFixed(1);
    if (semInput && lastChanged !== 'sem') semInput.value = metrics.sem.toFixed(3);
    if (nInput && lastChanged !== 'n') nInput.value = metrics.n;
    if (powerInput && lastChanged !== 'power') powerInput.value = metrics.power.toFixed(2);
    if (betaInput && lastChanged !== 'beta') betaInput.value = metrics.beta.toFixed(2);
    if (deltaInput && lastChanged !== 'delta') deltaInput.value = metrics.delta.toFixed(1);
    if (alphaInput && lastChanged !== 'alpha') alphaInput.value = metrics.alpha.toFixed(3);

    // Update Slider Display Labels
    const sdValEl = document.getElementById('powerSDVal');
    const semValEl = document.getElementById('powerSEMVal');
    const nValEl = document.getElementById('powerNVal');
    const powerValEl = document.getElementById('powerPowerVal');
    const betaValEl = document.getElementById('powerBetaVal');
    const deltaValEl = document.getElementById('powerDeltaVal');
    const alphaValEl = document.getElementById('powerAlphaVal');

    if (sdValEl) sdValEl.innerText = metrics.sd.toFixed(2);
    if (semValEl) semValEl.innerText = metrics.sem.toFixed(3);
    if (nValEl) nValEl.innerText = `n = ${metrics.n} (N = ${metrics.totalN})`;
    if (powerValEl) {
      powerValEl.innerText = `${(metrics.power * 100).toFixed(1)}%`;
      powerValEl.style.color = metrics.power >= 0.80 ? '#10b981' : (metrics.power >= 0.60 ? '#f59e0b' : '#ef4444');
    }
    if (betaValEl) {
      betaValEl.innerText = `${(metrics.beta * 100).toFixed(1)}%`;
      betaValEl.style.color = metrics.beta <= 0.20 ? '#10b981' : (metrics.beta <= 0.40 ? '#f59e0b' : '#ef4444');
    }
    if (deltaValEl) deltaValEl.innerText = metrics.delta.toFixed(2);
    if (alphaValEl) alphaValEl.innerText = metrics.alpha.toFixed(3);

    const dSub = document.getElementById('powerCohensDSub');
    if (dSub) dSub.innerText = `Cohen's d = ${metrics.cohensD.toFixed(2)}`;
    const zSub = document.getElementById('powerZCritSub');
    if (zSub) zSub.innerText = `z_crit = ${metrics.zCrit.toFixed(3)}`;

    const g1Summary = document.getElementById('powerG1Summary');
    if (g1Summary) g1Summary.innerText = `SD = ${metrics.sd.toFixed(2)} | SEM = ${metrics.sem.toFixed(3)}`;
    const g2Summary = document.getElementById('powerG2Summary');
    if (g2Summary) g2Summary.innerText = `Power = ${(metrics.power * 100).toFixed(1)}% | β = ${(metrics.beta * 100).toFixed(1)}%`;

    // Update 6 Status Metric Cards
    const statusPwrVal = document.getElementById('powerStatusPowerVal');
    const statusPwrSub = document.getElementById('powerStatusPowerSub');
    if (statusPwrVal) {
      statusPwrVal.innerText = `${(metrics.power * 100).toFixed(1)}%`;
      statusPwrVal.style.color = metrics.power >= 0.80 ? '#10b981' : (metrics.power >= 0.60 ? '#f59e0b' : '#ef4444');
    }
    if (statusPwrSub) {
      statusPwrSub.innerText = `Target: ≥ 80% (${metrics.powerRating})`;
    }

    const statusBetaVal = document.getElementById('powerStatusBetaVal');
    if (statusBetaVal) {
      statusBetaVal.innerText = `${(metrics.beta * 100).toFixed(1)}%`;
      statusBetaVal.style.color = metrics.beta <= 0.20 ? '#10b981' : (metrics.beta <= 0.40 ? '#f59e0b' : '#ef4444');
    }

    const statusNVal = document.getElementById('powerStatusNVal');
    const statusNSub = document.getElementById('powerStatusNSub');
    if (statusNVal) statusNVal.innerText = `n = ${metrics.n}`;
    if (statusNSub) statusNSub.innerText = `Total N = ${metrics.totalN} subjects`;

    const statusDeltaVal = document.getElementById('powerStatusDeltaVal');
    const statusDeltaSub = document.getElementById('powerStatusDeltaSub');
    if (statusDeltaVal) statusDeltaVal.innerText = `Δ = ${metrics.delta.toFixed(2)}`;
    if (statusDeltaSub) {
      const dLabel = metrics.cohensD >= 0.8 ? 'Large effect' : (metrics.cohensD >= 0.5 ? 'Medium effect' : 'Small effect');
      statusDeltaSub.innerText = `Cohen's d = ${metrics.cohensD.toFixed(2)} (${dLabel})`;
    }

    const statusSEMVal = document.getElementById('powerStatusSEMVal');
    const statusSEMSub = document.getElementById('powerStatusSEMSub');
    if (statusSEMVal) statusSEMVal.innerText = `SEM = ${metrics.sem.toFixed(3)}`;
    if (statusSEMSub) statusSEMSub.innerText = `SE_diff: ${metrics.seDiff.toFixed(3)} (σ·√(2/n))`;

    const statusCutoffVal = document.getElementById('powerStatusCutoffVal');
    const statusCutoffSub = document.getElementById('powerStatusCutoffSub');
    if (statusCutoffVal) statusCutoffVal.innerText = `xcrit = ${metrics.xCrit.toFixed(3)}`;
    if (statusCutoffSub) statusCutoffSub.innerText = `α = ${metrics.alpha.toFixed(3)} (Type I: ${(metrics.alpha * 100).toFixed(1)}%)`;

    // Update Pedagogical Text
    const pedaEl = document.getElementById('powerPedagogyText');
    if (pedaEl) pedaEl.innerText = metrics.explanation;

    // Update Chart Title
    const titleEl = document.getElementById('powerChartTitle');
    if (titleEl) {
      if (viewMode === 'curve') {
        titleEl.innerText = `Statistical Power Curve: Power (1 − β) vs Sample Size n (Δ = ${metrics.delta.toFixed(2)}, SD = ${metrics.sd.toFixed(2)})`;
      } else if (viewMode === 'matrix') {
        titleEl.innerText = `2×2 Decision Error Matrix: True State vs Clinical Statistical Decision`;
      } else {
        titleEl.innerText = `Dual Distribution: Null Hypothesis H₀ vs True Effect H₁ (Power & Beta Shading)`;
      }
    }

    // Render Canvas
    if (this.engines['teachingPowerCanvas']) {
      Plots.renderPowerSimulation(this.engines['teachingPowerCanvas'], metrics);
    }

    // Cache
    if (!this.results.teaching) this.results.teaching = {};
    this.results.teaching.power = metrics;
  }

  animatePowerGain() {
    const btn = document.getElementById('powerAnimateBtn');
    if (this.powerAnimationTimer) {
      clearInterval(this.powerAnimationTimer);
      this.powerAnimationTimer = null;
      if (btn) btn.innerText = '▶ Animate Power Gain';
      return;
    }

    const sampleSizes = [10, 16, 24, 34, 46, 64, 86, 112, 144];
    let currentIndex = 0;
    const currentN = parseInt(document.getElementById('powerNRange')?.value) || 64;
    const startIdx = sampleSizes.findIndex(n => n >= currentN);
    if (startIdx >= 0 && startIdx < sampleSizes.length - 1) currentIndex = startIdx;

    if (btn) btn.innerText = '⏸ Pause Animation';

    this.powerAnimationTimer = setInterval(() => {
      currentIndex++;
      if (currentIndex >= sampleSizes.length) {
        clearInterval(this.powerAnimationTimer);
        this.powerAnimationTimer = null;
        if (btn) btn.innerText = '▶ Replay Animation';
        return;
      }
      const nextN = sampleSizes[currentIndex];
      const nRange = document.getElementById('powerNRange');
      if (nRange) nRange.value = nextN;
      this.runPowerSimulation({ n: nextN }, 'n');
    }, 550);
  }

  runBayesianSimulation(overrideParams = {}) {
    // Read input values with fallbacks
    const priorInput = document.getElementById('bayesPriorRange');
    const sampleInput = document.getElementById('bayesSampleRange');
    const likInput = document.getElementById('bayesLikelihoodRange');
    const fpInput = document.getElementById('bayesFalsePosRange');

    let prior = overrideParams.prior !== undefined 
      ? overrideParams.prior 
      : (priorInput ? parseFloat(priorInput.value) / 100 : 0.0476);
    let sampleSize = overrideParams.sampleSize !== undefined 
      ? overrideParams.sampleSize 
      : (sampleInput ? parseInt(sampleInput.value) : 210);
    let likelihood = overrideParams.likelihood !== undefined 
      ? overrideParams.likelihood 
      : (likInput ? parseFloat(likInput.value) / 100 : 0.40);
    let falsePositive = overrideParams.falsePositive !== undefined 
      ? overrideParams.falsePositive 
      : (fpInput ? parseFloat(fpInput.value) / 100 : 0.10);

    const viewMode = this.currentBayesView || 'square';
    const preset = overrideParams.preset || this.currentBayesPreset || 'steve';

    // Smart percentage formatter preserving exact user-entered precision (e.g. 11.8% instead of rounding to 12%)
    const formatPercentSmart = (rate) => {
      if (rate === undefined || rate === null || isNaN(rate)) return '0.0%';
      const pct = rate * 100;
      const rounded1 = Math.round(pct * 10) / 10;
      if (Math.abs(pct - rounded1) < 1e-4) {
        return `${rounded1.toFixed(1)}%`;
      }
      const rounded2 = Math.round(pct * 100) / 100;
      if (Math.abs(pct - rounded2) < 1e-4) {
        return `${rounded2.toFixed(2)}%`;
      }
      return `${parseFloat(pct.toFixed(3))}%`;
    };

    const getSliderDisplayVal = (rate) => {
      const pct = rate * 100;
      const rounded1 = Math.round(pct * 10) / 10;
      if (Math.abs(pct - rounded1) < 1e-4) {
        return rounded1.toFixed(1);
      }
      const rounded2 = Math.round(pct * 100) / 100;
      if (Math.abs(pct - rounded2) < 1e-4) {
        return rounded2.toFixed(2);
      }
      return parseFloat(pct.toFixed(3)).toString();
    };

    // Synchronize UI slider values and numeric display labels
    if (priorInput && overrideParams.prior !== undefined) {
      priorInput.value = getSliderDisplayVal(prior);
    }
    if (sampleInput && overrideParams.sampleSize !== undefined) {
      sampleInput.value = sampleSize;
    }
    if (likInput && overrideParams.likelihood !== undefined) {
      likInput.value = getSliderDisplayVal(likelihood);
    }
    if (fpInput && overrideParams.falsePositive !== undefined) {
      fpInput.value = getSliderDisplayVal(falsePositive);
    }

    const priorValEl = document.getElementById('bayesPriorVal');
    if (priorValEl) priorValEl.innerText = formatPercentSmart(prior);

    const sampleValEl = document.getElementById('bayesSampleVal');
    if (sampleValEl) sampleValEl.innerText = `N = ${sampleSize} people`;

    const likValEl = document.getElementById('bayesLikelihoodVal');
    if (likValEl) likValEl.innerText = formatPercentSmart(likelihood);

    const fpValEl = document.getElementById('bayesFalsePosVal');
    if (fpValEl) fpValEl.innerText = formatPercentSmart(falsePositive);

    const priorSumEl = document.getElementById('bayesPriorSummary');
    if (priorSumEl) {
      const priorOddsRatio = (1 - prior) / prior;
      priorSumEl.innerText = `Prior: ${formatPercentSmart(prior)} (1:${priorOddsRatio.toFixed(1)})`;
    }

    const likSumEl = document.getElementById('bayesLikelihoodSummary');
    if (likSumEl) {
      likSumEl.innerText = `P(E|H): ${formatPercentSmart(likelihood)} | P(E|¬H): ${formatPercentSmart(falsePositive)}`;
    }

    // Calculate metrics
    const metrics = Teaching.bayesianSimulation.getMetrics({
      prior,
      likelihood,
      falsePositive,
      sampleSize,
      viewMode,
      preset
    });

    // Update real-time Metric Cards
    const statusPriorVal = document.getElementById('bayesStatusPriorVal');
    const statusPriorSub = document.getElementById('bayesStatusPriorSub');
    if (statusPriorVal) statusPriorVal.innerText = formatPercentSmart(metrics.prior);
    if (statusPriorSub) statusPriorSub.innerText = `Prior Odds: 1 : ${(1 / metrics.priorOdds).toFixed(1)}`;

    const statusLikVal = document.getElementById('bayesStatusLikelihoodVal');
    const statusLikSub = document.getElementById('bayesStatusLikelihoodSub');
    if (statusLikVal) statusLikVal.innerText = formatPercentSmart(metrics.likelihood);
    if (statusLikSub) statusLikSub.innerText = `True Positive Probability`;

    const statusFPVal = document.getElementById('bayesStatusFalsePosVal');
    const statusFPSub = document.getElementById('bayesStatusFalsePosSub');
    if (statusFPVal) statusFPVal.innerText = formatPercentSmart(metrics.falsePositive);
    if (statusFPSub) statusFPSub.innerText = `False Positive Probability`;

    const statusPEVal = document.getElementById('bayesStatusPEvidenceVal');
    const statusPESub = document.getElementById('bayesStatusPEvidenceSub');
    if (statusPEVal) statusPEVal.innerText = `${(metrics.pEvidence * 100).toFixed(2)}%`;
    if (statusPESub) statusPESub.innerText = `Marginal / Area (${metrics.countTotalE} of ${metrics.sampleSize})`;

    const statusBFVal = document.getElementById('bayesStatusBFVal');
    const statusBFSub = document.getElementById('bayesStatusBFSub');
    if (statusBFVal) statusBFVal.innerText = `${metrics.bayesFactor.toFixed(2)}×`;
    if (statusBFSub) statusBFSub.innerText = metrics.evidenceRating;

    const statusPostVal = document.getElementById('bayesStatusPosteriorVal');
    const statusPostSub = document.getElementById('bayesStatusPosteriorSub');
    if (statusPostVal) statusPostVal.innerText = `${(metrics.posterior * 100).toFixed(1)}%`;
    if (statusPostSub) {
      const sign = metrics.beliefShift >= 0 ? '+' : '';
      statusPostSub.innerText = `Belief Shift: ${sign}${(metrics.beliefShift * 100).toFixed(1)}%`;
    }

    // Update Pedagogical Text
    const pedaEl = document.getElementById('bayesPedagogyText');
    if (pedaEl) pedaEl.innerText = metrics.explanation;

    const bayesReportEl = document.getElementById('teachingBayesReportText');
    if (bayesReportEl) bayesReportEl.innerText = metrics.explanation;

    // Update Chart Title
    const titleEl = document.getElementById('teachingBayesChartTitle');
    if (titleEl) {
      if (viewMode === 'sample') {
        titleEl.innerText = `Representative Sample (N = ${metrics.sampleSize}): Natural Counts Frequency Tree (${metrics.countHAndE} vs ${metrics.countNotHAndE})`;
      } else if (viewMode === 'odds') {
        titleEl.innerText = `Odds Form & Bayes Factor: Prior Odds (1 : ${(1 / metrics.priorOdds).toFixed(1)}) × LR (${metrics.bayesFactor.toFixed(2)}×) = Posterior Odds (${metrics.posteriorOdds.toFixed(2)})`;
      } else if (viewMode === 'sequential') {
        titleEl.innerText = `Sequential Bayesian Updating: Compounding Evidence Trajectory (P₀ = ${(metrics.prior * 100).toFixed(1)}% → P₄ = ${(metrics.trajectory[4].p * 100).toFixed(1)}%)`;
      } else {
        titleEl.innerText = `3Blue1Brown 1×1 Unit Square: Visualizing Bayes' Theorem as Proportions of Area (P(H|E) = ${(metrics.posterior * 100).toFixed(1)}%)`;
      }
    }

    // Render Canvas
    if (this.engines['teachingBayesCanvas']) {
      Plots.renderBayesianSimulation(this.engines['teachingBayesCanvas'], metrics);
    }

    // Cache
    if (!this.results.teaching) this.results.teaching = {};
    this.results.teaching.bayes = metrics;
  }

  animateBayesRestriction() {
    const btn = document.getElementById('bayesAnimateBtn');
    if (this.bayesAnimationTimer) {
      clearInterval(this.bayesAnimationTimer);
      this.bayesAnimationTimer = null;
      if (btn) btn.innerText = '▶ Animate Restriction';
      return;
    }

    // Cycle through view modes
    const views = ['square', 'sample', 'odds', 'sequential'];
    let idx = views.indexOf(this.currentBayesView || 'square');
    if (btn) btn.innerText = '⏸ Pause Animation';

    this.bayesAnimationTimer = setInterval(() => {
      idx = (idx + 1) % views.length;
      const nextView = views[idx];
      this.currentBayesView = nextView;
      document.querySelectorAll('.teaching-bayes-view-btn').forEach(b => {
        const isActive = b.dataset.view === nextView;
        b.classList.toggle('btn-primary', isActive);
        b.classList.toggle('active', isActive);
        b.classList.toggle('btn-secondary', !isActive);
      });
      this.runBayesianSimulation();
    }, 1200);
  }

  // ==========================================
  // RANDOMISER WORKFLOW & EVENT HANDLERS
  // ==========================================
  initRandomiser() {
    this.randomiserState = {
      mode: localStorage.getItem('statis_gravity_randomiser_mode') || 'simple',
      simple: {
        history: [],
        nextId: 1,
        labelA: 'Group A (Treatment)',
        labelB: 'Group B (Control)'
      },
      block: {
        history: [],
        nextId: 1,
        blockSize: 4,
        targetN: 40,
        labelA: 'Group A (Intervention)',
        labelB: 'Group B (Control)',
        activeBlock: null
      }
    };

    try {
      const savedSimple = localStorage.getItem('statis_gravity_randomiser_simple');
      if (savedSimple) {
        const parsed = JSON.parse(savedSimple);
        if (parsed && Array.isArray(parsed.history)) {
          this.randomiserState.simple.history = parsed.history;
          this.randomiserState.simple.nextId = parsed.nextId || (parsed.history.length + 1);
          if (parsed.labelA) this.randomiserState.simple.labelA = parsed.labelA;
          if (parsed.labelB) this.randomiserState.simple.labelB = parsed.labelB;
        }
      }
      const savedBlock = localStorage.getItem('statis_gravity_randomiser_block');
      if (savedBlock) {
        const parsed = JSON.parse(savedBlock);
        if (parsed && Array.isArray(parsed.history)) {
          this.randomiserState.block.history = parsed.history;
          this.randomiserState.block.nextId = parsed.nextId || (parsed.history.length + 1);
          if (parsed.blockSize) this.randomiserState.block.blockSize = parsed.blockSize;
          if (parsed.targetN) this.randomiserState.block.targetN = parsed.targetN;
          if (parsed.labelA) this.randomiserState.block.labelA = parsed.labelA;
          if (parsed.labelB) this.randomiserState.block.labelB = parsed.labelB;
          if (parsed.activeBlock) this.randomiserState.block.activeBlock = parsed.activeBlock;
        }
      }
    } catch (err) {
      console.warn('Could not load randomiser state from localStorage:', err);
    }

    // Sync form input fields
    const sLabelA = document.getElementById('randomiserSimpleLabelA');
    const sLabelB = document.getElementById('randomiserSimpleLabelB');
    const sNextId = document.getElementById('randomiserSimpleNextId');
    if (sLabelA) sLabelA.value = this.randomiserState.simple.labelA;
    if (sLabelB) sLabelB.value = this.randomiserState.simple.labelB;
    if (sNextId) sNextId.value = this.randomiserState.simple.nextId;

    const bSize = document.getElementById('randomiserBlockSizeSelect');
    const bTargetN = document.getElementById('randomiserBlockTargetN');
    const bLabelA = document.getElementById('randomiserBlockLabelA');
    const bLabelB = document.getElementById('randomiserBlockLabelB');
    if (bSize) bSize.value = String(this.randomiserState.block.blockSize);
    if (bTargetN) bTargetN.value = this.randomiserState.block.targetN;
    if (bLabelA) bLabelA.value = this.randomiserState.block.labelA;
    if (bLabelB) bLabelB.value = this.randomiserState.block.labelB;

    // Event listeners for Simple inputs
    sLabelA?.addEventListener('input', (e) => {
      this.randomiserState.simple.labelA = e.target.value.trim() || 'Group A';
      this.saveRandomiserState('simple');
    });
    sLabelB?.addEventListener('input', (e) => {
      this.randomiserState.simple.labelB = e.target.value.trim() || 'Group B';
      this.saveRandomiserState('simple');
    });
    sNextId?.addEventListener('change', (e) => {
      const val = parseInt(e.target.value, 10);
      if (!isNaN(val) && val >= 1) {
        this.randomiserState.simple.nextId = val;
        const btnId = document.getElementById('randomiserSimpleBtnIdText');
        if (btnId) btnId.innerText = String(val);
        this.saveRandomiserState('simple');
      }
    });

    // Event listeners for Block inputs
    bSize?.addEventListener('change', (e) => {
      const newSize = parseInt(e.target.value, 10);
      this.randomiserState.block.blockSize = newSize;
      if (!this.randomiserState.block.activeBlock || this.randomiserState.block.activeBlock.currentIndex === 0) {
        this.randomiserState.block.activeBlock = null;
      }
      this.updateRandomiserBlockUI();
      this.saveRandomiserState('block');
    });
    bTargetN?.addEventListener('input', (e) => {
      const val = parseInt(e.target.value, 10);
      if (!isNaN(val) && val >= 2) {
        this.randomiserState.block.targetN = val;
        this.updateRandomiserBlockUI();
        this.saveRandomiserState('block');
      }
    });
    bLabelA?.addEventListener('input', (e) => {
      this.randomiserState.block.labelA = e.target.value.trim() || 'Group A';
      this.saveRandomiserState('block');
    });
    bLabelB?.addEventListener('input', (e) => {
      this.randomiserState.block.labelB = e.target.value.trim() || 'Group B';
      this.saveRandomiserState('block');
    });

    // Sub-Navigation mode switching
    const modeSimpleBtn = document.getElementById('randomiserModeSimpleBtn');
    const modeBlockBtn = document.getElementById('randomiserModeBlockBtn');
    const simplePanel = document.getElementById('randomiserSimplePanel');
    const blockPanel = document.getElementById('randomiserBlockPanel');

    const setMode = (mode) => {
      this.randomiserState.mode = mode;
      localStorage.setItem('statis_gravity_randomiser_mode', mode);

      if (mode === 'simple') {
        modeSimpleBtn?.classList.add('active');
        modeSimpleBtn?.setAttribute('aria-selected', 'true');
        modeBlockBtn?.classList.remove('active');
        modeBlockBtn?.setAttribute('aria-selected', 'false');
        if (simplePanel) simplePanel.style.display = 'block';
        if (blockPanel) blockPanel.style.display = 'none';
      } else {
        modeBlockBtn?.classList.add('active');
        modeBlockBtn?.setAttribute('aria-selected', 'true');
        modeSimpleBtn?.classList.remove('active');
        modeSimpleBtn?.setAttribute('aria-selected', 'false');
        if (simplePanel) simplePanel.style.display = 'none';
        if (blockPanel) blockPanel.style.display = 'block';
      }
      this.renderRandomiserAuditTable();
    };

    modeSimpleBtn?.addEventListener('click', () => setMode('simple'));
    modeBlockBtn?.addEventListener('click', () => setMode('block'));

    // Action Buttons
    document.getElementById('randomiserSimpleGenerateBtn')?.addEventListener('click', () => {
      this.generateSimpleAllocationAction();
    });

    document.getElementById('randomiserBlockAssignBtn')?.addEventListener('click', () => {
      this.assignBlockAllocationAction();
    });

    document.getElementById('randomiserExportCsvBtn')?.addEventListener('click', () => {
      this.exportRandomiserCSV();
    });

    document.getElementById('randomiserCopyTableBtn')?.addEventListener('click', () => {
      this.copyRandomiserAuditTable();
    });

    // Reset Modal Handlers
    const resetModal = document.getElementById('randomiserResetModal');
    document.getElementById('randomiserResetBtn')?.addEventListener('click', () => {
      if (resetModal) resetModal.classList.remove('hidden');
    });
    document.getElementById('randomiserResetCancelBtn')?.addEventListener('click', () => {
      if (resetModal) resetModal.classList.add('hidden');
    });
    document.getElementById('randomiserResetConfirmBtn')?.addEventListener('click', () => {
      this.resetRandomiserSession();
      if (resetModal) resetModal.classList.add('hidden');
    });

    // Initial View Setup
    setMode(this.randomiserState.mode);
    this.updateRandomiserSimpleUI();
    this.updateRandomiserBlockUI();
    this.renderRandomiserAuditTable();
  }

  generateSimpleAllocationAction() {
    const pid = this.randomiserState.simple.nextId;
    const record = Randomiser.generateSimpleAllocation({
      participantId: pid,
      labelA: this.randomiserState.simple.labelA,
      labelB: this.randomiserState.simple.labelB
    });

    this.randomiserState.simple.history.push(record);
    this.randomiserState.simple.nextId += 1;

    const sNextIdEl = document.getElementById('randomiserSimpleNextId');
    if (sNextIdEl) sNextIdEl.value = this.randomiserState.simple.nextId;

    this.saveRandomiserState('simple');
    this.updateRandomiserSimpleUI(record);
    this.renderRandomiserAuditTable();
  }

  updateRandomiserSimpleUI(lastRecord = null) {
    const sNextId = this.randomiserState.simple.nextId;
    const btnIdText = document.getElementById('randomiserSimpleBtnIdText');
    if (btnIdText) btnIdText.innerText = String(sNextId);

    const history = this.randomiserState.simple.history;
    const rec = lastRecord || (history.length > 0 ? history[history.length - 1] : null);

    const idleView = document.getElementById('randomiserSimpleIdleView');
    const activeView = document.getElementById('randomiserSimpleActiveView');
    const revealCard = document.getElementById('randomiserSimpleRevealCard');
    const badge = document.getElementById('randomiserSimpleBadge');
    const resId = document.getElementById('randomiserSimpleResultId');
    const resNum = document.getElementById('randomiserSimpleResultNum');
    const resParity = document.getElementById('randomiserSimpleResultParity');
    const resTime = document.getElementById('randomiserSimpleResultTime');
    const tsText = document.getElementById('randomiserSimpleTimestampText');

    if (rec) {
      if (idleView) idleView.style.display = 'none';
      if (activeView) activeView.style.display = 'block';
      if (resId) resId.innerText = `#${rec.participantId}`;
      if (badge) {
        badge.innerText = rec.groupLabel;
        badge.className = `allocation-badge-large ${rec.groupKey === 'A' ? 'badge-group-a' : 'badge-group-b'}`;
      }
      if (revealCard) {
        revealCard.className = `allocation-reveal-card ${rec.groupKey === 'A' ? 'revealed-a' : 'revealed-b'}`;
      }
      if (resNum) resNum.innerText = String(rec.randomNumber);
      if (resParity) {
        resParity.innerText = rec.parity;
        resParity.style.color = rec.groupKey === 'A' ? 'var(--cyan-primary)' : 'var(--emerald-primary)';
      }
      if (resTime) resTime.innerText = rec.displayTime || '';
      if (tsText) tsText.innerText = `Last assigned at ${rec.displayTime || ''}`;
    } else {
      if (idleView) idleView.style.display = 'block';
      if (activeView) activeView.style.display = 'none';
      if (revealCard) revealCard.className = 'allocation-reveal-card';
      if (tsText) tsText.innerText = 'No allocation yet';
    }

    // Summary metrics
    const summary = Randomiser.computeSummary(history);
    const totalEl = document.getElementById('randomiserSimpleTotalN');
    const countAEl = document.getElementById('randomiserSimpleCountA');
    const countBEl = document.getElementById('randomiserSimpleCountB');
    const ratioEl = document.getElementById('randomiserSimpleRatio');
    const pctAEl = document.getElementById('randomiserSimplePctA');
    const pctBEl = document.getElementById('randomiserSimplePctB');
    const barA = document.getElementById('randomiserSimpleBarA');
    const barB = document.getElementById('randomiserSimpleBarB');

    if (totalEl) totalEl.innerText = String(summary.total);
    if (countAEl) countAEl.innerText = `${summary.countA} (${summary.pctA.toFixed(0)}%)`;
    if (countBEl) countBEl.innerText = `${summary.countB} (${summary.pctB.toFixed(0)}%)`;
    if (ratioEl) ratioEl.innerText = summary.ratioStr;
    if (pctAEl) pctAEl.innerText = `${summary.pctA.toFixed(0)}%`;
    if (pctBEl) pctBEl.innerText = `${summary.pctB.toFixed(0)}%`;

    const barAWidth = summary.total === 0 ? 50 : Math.max(5, Math.min(95, summary.pctA));
    const barBWidth = summary.total === 0 ? 50 : (100 - barAWidth);
    if (barA) barA.style.width = `${barAWidth}%`;
    if (barB) barB.style.width = `${barBWidth}%`;
  }

  assignBlockAllocationAction() {
    const pid = this.randomiserState.block.nextId;
    const res = Randomiser.assignNextInBlock(
      this.randomiserState.block.activeBlock,
      pid,
      {
        blockSize: this.randomiserState.block.blockSize,
        labelA: this.randomiserState.block.labelA,
        labelB: this.randomiserState.block.labelB
      }
    );

    this.randomiserState.block.activeBlock = res.updatedBlock;
    this.randomiserState.block.history.push(res.allocationRecord);
    this.randomiserState.block.nextId += 1;

    this.saveRandomiserState('block');
    this.updateRandomiserBlockUI(res.allocationRecord);
    this.renderRandomiserAuditTable();
  }

  updateRandomiserBlockUI(lastRecord = null) {
    const bNextId = this.randomiserState.block.nextId;
    const btnIdText = document.getElementById('randomiserBlockBtnIdText');
    if (btnIdText) btnIdText.innerText = `Participant #${bNextId}`;

    const history = this.randomiserState.block.history;
    const rec = lastRecord || (history.length > 0 ? history[history.length - 1] : null);

    const idleView = document.getElementById('randomiserBlockIdleView');
    const activeView = document.getElementById('randomiserBlockActiveView');
    const revealCard = document.getElementById('randomiserBlockRevealCard');
    const badge = document.getElementById('randomiserBlockBadge');
    const resId = document.getElementById('randomiserBlockResultId');
    const resBlockNum = document.getElementById('randomiserBlockResultBlockNum');
    const resSlot = document.getElementById('randomiserBlockResultSlot');
    const resTime = document.getElementById('randomiserBlockResultTime');

    if (rec) {
      if (idleView) idleView.style.display = 'none';
      if (activeView) activeView.style.display = 'block';
      if (resId) resId.innerText = `#${rec.participantId}`;
      if (badge) {
        badge.innerText = rec.groupLabel;
        badge.className = `allocation-badge-large ${rec.groupKey === 'A' ? 'badge-group-a' : 'badge-group-b'}`;
      }
      if (revealCard) {
        revealCard.className = `allocation-reveal-card ${rec.groupKey === 'A' ? 'revealed-a' : 'revealed-b'}`;
      }
      if (resBlockNum) resBlockNum.innerText = String(rec.blockNumber);
      if (resSlot) resSlot.innerText = `${rec.slotInBlock} of ${rec.blockSize}`;
      if (resTime) resTime.innerText = rec.displayTime || '';
    } else {
      if (idleView) idleView.style.display = 'block';
      if (activeView) activeView.style.display = 'none';
      if (revealCard) revealCard.className = 'allocation-reveal-card';
    }

    // Block slots visual tracker
    const activeBlock = this.randomiserState.block.activeBlock;
    const bSize = this.randomiserState.block.blockSize;
    const slotsContainer = document.getElementById('randomiserBlockSlotsContainer');
    const progressText = document.getElementById('randomiserBlockProgressText');
    const statusBadge = document.getElementById('randomiserBlockStatusBadge');

    if (slotsContainer) {
      slotsContainer.innerHTML = '';
      const currentBlockNum = activeBlock ? activeBlock.blockNumber : (Math.floor(history.length / bSize) + 1);
      if (statusBadge) statusBadge.innerText = `Block #${currentBlockNum} Active`;

      if (activeBlock && activeBlock.slots) {
        const assignedCount = activeBlock.currentIndex;
        if (progressText) progressText.innerText = `Slot ${assignedCount} of ${activeBlock.blockSize} assigned`;

        activeBlock.slots.forEach(slot => {
          const pill = document.createElement('div');
          if (slot.assigned) {
            pill.className = `block-slot-pill ${slot.groupKey === 'A' ? 'assigned-a' : 'assigned-b'}`;
            pill.innerHTML = `✓ Slot ${slot.slotIndex}: <strong>${slot.groupKey}</strong> (P#${slot.participantId})`;
          } else if (slot.slotIndex === activeBlock.currentIndex + 1) {
            pill.className = 'block-slot-pill concealed-slot active-slot';
            pill.innerHTML = `★ Next Slot ${slot.slotIndex}: 🔒 Concealed`;
          } else {
            pill.className = 'block-slot-pill concealed-slot';
            pill.innerHTML = `Slot ${slot.slotIndex}: 🔒 Pending`;
          }
          slotsContainer.appendChild(pill);
        });
      } else {
        if (progressText) progressText.innerText = `Slot 0 of ${bSize} assigned`;
        for (let i = 1; i <= bSize; i++) {
          const pill = document.createElement('div');
          pill.className = 'block-slot-pill concealed-slot';
          pill.innerHTML = `Slot ${i}: 🔒 Concealed`;
          slotsContainer.appendChild(pill);
        }
      }
    }

    // Summary metrics
    const summary = Randomiser.computeSummary(history);
    const totalBlockEl = document.getElementById('randomiserBlockTotalN');
    const countABlockEl = document.getElementById('randomiserBlockCountA');
    const countBBlockEl = document.getElementById('randomiserBlockCountB');
    const completedBlocksEl = document.getElementById('randomiserBlockCompletedCount');
    const pctABlockEl = document.getElementById('randomiserBlockPctA');
    const pctBBlockEl = document.getElementById('randomiserBlockPctB');
    const barABlock = document.getElementById('randomiserBlockBarA');
    const barBBlock = document.getElementById('randomiserBlockBarB');

    if (totalBlockEl) totalBlockEl.innerText = `${summary.total} / ${this.randomiserState.block.targetN}`;
    if (countABlockEl) countABlockEl.innerText = `${summary.countA} (${summary.pctA.toFixed(0)}%)`;
    if (countBBlockEl) countBBlockEl.innerText = `${summary.countB} (${summary.pctB.toFixed(0)}%)`;
    const completedBlocks = Math.floor(summary.total / bSize);
    if (completedBlocksEl) completedBlocksEl.innerText = String(completedBlocks);
    if (pctABlockEl) pctABlockEl.innerText = `${summary.pctA.toFixed(0)}%`;
    if (pctBBlockEl) pctBBlockEl.innerText = `${summary.pctB.toFixed(0)}%`;

    const barAWidth = summary.total === 0 ? 50 : Math.max(5, Math.min(95, summary.pctA));
    const barBWidth = summary.total === 0 ? 50 : (100 - barAWidth);
    if (barABlock) barABlock.style.width = `${barAWidth}%`;
    if (barBBlock) barBBlock.style.width = `${barBWidth}%`;
  }

  renderRandomiserAuditTable() {
    const mode = this.randomiserState.mode;
    const history = this.randomiserState[mode].history;
    const tbody = document.getElementById('randomiserAuditTableBody');
    const countBadge = document.getElementById('randomiserAuditCountBadge');

    if (!tbody) return;

    if (countBadge) {
      countBadge.innerText = `${history.length} Record${history.length === 1 ? '' : 's'}`;
    }

    if (history.length === 0) {
      tbody.innerHTML = `
        <tr id="randomiserAuditEmptyRow">
          <td colspan="6" style="text-align: center; padding: 2.5rem 1rem; color: var(--text-dim);">
            <div style="font-size: 1.8rem; margin-bottom: 0.4rem;">🎲</div>
            <div style="font-weight: 600; color: var(--text-muted);">No participants randomized yet in ${mode === 'simple' ? 'Simple' : 'Block'} mode</div>
            <div style="font-size: 0.78rem; margin-top: 0.2rem;">Click "Generate / Assign Next Participant" above to initiate sequence allocation.</div>
          </td>
        </tr>
      `;
      return;
    }

    tbody.innerHTML = '';
    // Render in reverse chronological order (newest first)
    for (let i = history.length - 1; i >= 0; i--) {
      const item = history[i];
      const tr = document.createElement('tr');
      tr.style.borderBottom = '1px solid var(--border-subtle)';

      const drawCol = mode === 'simple'
        ? `<span class="badge ${item.parity === 'Odd' ? 'badge-sig' : 'badge-neutral'}">${item.randomNumber} (${item.parity})</span>`
        : `<span class="badge badge-sig">Block ${item.blockNumber} (Slot ${item.slotInBlock}/${item.blockSize})</span>`;

      tr.innerHTML = `
        <td style="padding: 0.6rem 0.85rem; font-weight: 800; color: var(--text-main);">#${item.participantId}</td>
        <td style="padding: 0.6rem 0.85rem; font-family: var(--font-mono); font-size: 0.78rem; color: var(--text-muted);">${item.displayTime || (item.timestamp ? item.timestamp.slice(11, 19) : '--')}</td>
        <td style="padding: 0.6rem 0.85rem; font-size: 0.8rem; color: var(--text-muted);">${item.method || (mode === 'simple' ? 'Simple (1-100)' : 'Block Permuted')}</td>
        <td style="padding: 0.6rem 0.85rem;">${drawCol}</td>
        <td style="padding: 0.6rem 0.85rem; font-weight: 800; color: ${item.groupKey === 'A' ? 'var(--cyan-primary)' : 'var(--emerald-primary)'}; font-size: 1.05rem;">${item.groupKey}</td>
        <td style="padding: 0.6rem 0.85rem; font-weight: 600; color: var(--text-main);">${item.groupLabel}</td>
      `;
      tbody.appendChild(tr);
    }
  }

  saveRandomiserState(type) {
    try {
      if (type === 'simple') {
        localStorage.setItem('statis_gravity_randomiser_simple', JSON.stringify({
          history: this.randomiserState.simple.history,
          nextId: this.randomiserState.simple.nextId,
          labelA: this.randomiserState.simple.labelA,
          labelB: this.randomiserState.simple.labelB
        }));
      } else if (type === 'block') {
        localStorage.setItem('statis_gravity_randomiser_block', JSON.stringify({
          history: this.randomiserState.block.history,
          nextId: this.randomiserState.block.nextId,
          blockSize: this.randomiserState.block.blockSize,
          targetN: this.randomiserState.block.targetN,
          labelA: this.randomiserState.block.labelA,
          labelB: this.randomiserState.block.labelB,
          activeBlock: this.randomiserState.block.activeBlock
        }));
      }
    } catch (err) {
      console.warn('Error saving randomiser state to localStorage:', err);
    }
  }

  exportRandomiserCSV() {
    const mode = this.randomiserState.mode;
    const history = this.randomiserState[mode].history;

    if (!history || history.length === 0) {
      alert(`No randomized records to export in ${mode === 'simple' ? 'Simple' : 'Block'} Randomization mode.`);
      return;
    }

    const csvContent = Randomiser.exportToCSV(history, mode);
    const blob = new Blob([csvContent], { type: 'text/csv;charset=utf-8;' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    const now = new Date();
    const dateStr = now.toISOString().slice(0, 10);
    a.download = `statis_gravity_${mode}_randomisation_audit_${dateStr}.csv`;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(url);
  }

  copyRandomiserAuditTable() {
    const mode = this.randomiserState.mode;
    const history = this.randomiserState[mode].history;
    const copyBtn = document.getElementById('randomiserCopyTableBtn');

    if (!history || history.length === 0) {
      alert('Audit table is empty. Randomize at least one participant first.');
      return;
    }

    const headers = mode === 'simple'
      ? ['Participant ID', 'Timestamp', 'Method', 'Random Draw (1-100)', 'Parity', 'Group Code', 'Group Label']
      : ['Participant ID', 'Timestamp', 'Method', 'Block Number', 'Block Size', 'Slot in Block', 'Group Code', 'Group Label'];

    const rows = history.map(item => mode === 'simple'
      ? [item.participantId, item.timestamp, item.method || 'Simple', item.randomNumber, item.parity, item.groupKey, item.groupLabel]
      : [item.participantId, item.timestamp, item.method || 'Block', item.blockNumber, item.blockSize, item.slotInBlock, item.groupKey, item.groupLabel]
    );

    const tsv = [headers.join('\t'), ...rows.map(r => r.join('\t'))].join('\n');

    if (navigator.clipboard && navigator.clipboard.writeText) {
      navigator.clipboard.writeText(tsv).then(() => {
        if (copyBtn) {
          const orig = copyBtn.innerText;
          copyBtn.innerText = '✓ Copied!';
          setTimeout(() => { copyBtn.innerText = orig; }, 1800);
        }
      });
    } else {
      const ta = document.createElement('textarea');
      ta.value = tsv;
      ta.style.position = 'fixed';
      ta.style.opacity = '0';
      document.body.appendChild(ta);
      ta.select();
      document.execCommand('copy');
      document.body.removeChild(ta);
      if (copyBtn) {
        const orig = copyBtn.innerText;
        copyBtn.innerText = '✓ Copied!';
        setTimeout(() => { copyBtn.innerText = orig; }, 1800);
      }
    }
  }

  resetRandomiserSession() {
    const mode = this.randomiserState.mode;
    if (mode === 'simple') {
      this.randomiserState.simple.history = [];
      this.randomiserState.simple.nextId = 1;
      const sNextIdEl = document.getElementById('randomiserSimpleNextId');
      if (sNextIdEl) sNextIdEl.value = '1';
      localStorage.removeItem('statis_gravity_randomiser_simple');
      this.updateRandomiserSimpleUI();
    } else {
      this.randomiserState.block.history = [];
      this.randomiserState.block.nextId = 1;
      this.randomiserState.block.activeBlock = null;
      localStorage.removeItem('statis_gravity_randomiser_block');
      this.updateRandomiserBlockUI();
    }
    this.renderRandomiserAuditTable();
  }

    // ==========================================
    // 14. PROPENSITY SCORE MATCHING (PSM) WORKFLOW
    // ==========================================
    initPsm() {
      this.psmCurrentScriptLang = 'python';
      this.psmActiveView = 'love';

      document.getElementById('psmSampleBtn')?.addEventListener('click', () => {
        this.loadPsmSample();
      });

      document.getElementById('psmClearBtn')?.addEventListener('click', () => {
        const input = document.getElementById('psmCsvInput');
        if (input) input.value = '';
        this.populatePsmVariables();
      });

      document.getElementById('psmComputeBtn')?.addEventListener('click', () => {
        this.runPsm();
      });

      const loveBtn = document.getElementById('psmViewLoveBtn');
      const overlapBtn = document.getElementById('psmViewOverlapBtn');
      const saveBtn = document.getElementById('psmSavePlotBtn');
      const loveContainer = document.getElementById('psmLovePlotContainer');
      const overlapContainer = document.getElementById('psmOverlapContainer');

      loveBtn?.addEventListener('click', () => {
        this.psmActiveView = 'love';
        loveBtn.classList.replace('btn-secondary', 'btn-primary');
        overlapBtn?.classList.replace('btn-primary', 'btn-secondary');
        if (loveContainer) loveContainer.style.display = 'block';
        if (overlapContainer) overlapContainer.style.display = 'none';
        if (saveBtn) saveBtn.dataset.canvasId = 'psmLovePlotCanvas';
        this.renderPsmPlots();
      });

      overlapBtn?.addEventListener('click', () => {
        this.psmActiveView = 'overlap';
        overlapBtn.classList.replace('btn-secondary', 'btn-primary');
        loveBtn?.classList.replace('btn-primary', 'btn-secondary');
        if (loveContainer) loveContainer.style.display = 'none';
        if (overlapContainer) overlapContainer.style.display = 'block';
        if (saveBtn) saveBtn.dataset.canvasId = 'psmOverlapCanvas';
        this.renderPsmPlots();
      });

      const pyBtn = document.getElementById('psmScriptLangPy');
      const rBtn = document.getElementById('psmScriptLangR');
      const stataBtn = document.getElementById('psmScriptLangStata');

      pyBtn?.addEventListener('click', () => this.switchPsmScript('python'));
      rBtn?.addEventListener('click', () => this.switchPsmScript('r'));
      stataBtn?.addEventListener('click', () => this.switchPsmScript('stata'));

      document.getElementById('psmCopyScriptBtn')?.addEventListener('click', () => {
        this.copyPsmScript();
      });

      document.getElementById('psmDownloadScriptBtn')?.addEventListener('click', () => {
        this.downloadPsmScript();
      });

      const csvInput = document.getElementById('psmCsvInput');
      csvInput?.addEventListener('input', () => {
        this.populatePsmVariables();
      });
      csvInput?.addEventListener('change', () => {
        this.populatePsmVariables();
      });

      document.getElementById('psmTreatmentSelect')?.addEventListener('change', () => {
        this.updatePsmCovariatesList();
      });
      document.getElementById('psmOutcomeSelect')?.addEventListener('change', () => {
        this.updatePsmCovariatesList();
      });
    }

    loadPsmSample() {
      const records = Psm.getSampleClinicalDataset();
      if (!records || records.length === 0) return;

      const headers = Object.keys(records[0]);
      const rows = records.map(r => headers.map(h => r[h]).join(','));
      const csv = [headers.join(','), ...rows].join('\\n');

      const input = document.getElementById('psmCsvInput');
      if (input) {
        input.value = csv;
      }
      this.populatePsmVariables();
      this.runPsm();
    }

    populatePsmVariables() {
      const csvText = document.getElementById('psmCsvInput')?.value || '';
      const lines = csvText.trim().split(/\\r?\\n/).filter(l => l.trim().length > 0);
      const badge = document.getElementById('psmDatasetBadge');
      if (lines.length < 2) {
        if (badge) badge.innerText = '0 rows';
        return;
      }

      if (badge) {
        badge.innerText = `N = ${lines.length - 1} patients`;
      }

      const delimiter = lines[0].includes('\\t') ? '\\t' : (lines[0].includes(';') ? ';' : ',');
      const headers = lines[0].split(delimiter).map(h => h.trim().replace(/^["']|["']$/g, ''));

      const tSelect = document.getElementById('psmTreatmentSelect');
      const oSelect = document.getElementById('psmOutcomeSelect');

      const currentT = tSelect?.value;
      const currentO = oSelect?.value;

      if (tSelect) {
        tSelect.innerHTML = '';
        headers.forEach(h => {
          const opt = document.createElement('option');
          opt.value = h;
          opt.textContent = h;
          if (h.toLowerCase().includes('treat') || h === 'treatment_col' || h === 'group') {
            opt.selected = true;
          }
          tSelect.appendChild(opt);
        });
        if (currentT && headers.includes(currentT)) tSelect.value = currentT;
      }

      if (oSelect) {
        oSelect.innerHTML = '';
        headers.forEach(h => {
          const opt = document.createElement('option');
          opt.value = h;
          opt.textContent = h;
          if (h.toLowerCase().includes('outcom') || h === 'outcome_col' || h.toLowerCase().includes('los')) {
            opt.selected = true;
          }
          oSelect.appendChild(opt);
        });
        if (currentO && headers.includes(currentO)) oSelect.value = currentO;
      }

      this.updatePsmCovariatesList(headers);
    }

    updatePsmCovariatesList(headers) {
      if (!headers) {
        const csvText = document.getElementById('psmCsvInput')?.value || '';
        const lines = csvText.trim().split(/\\r?\\n/).filter(l => l.trim().length > 0);
        if (lines.length < 2) return;
        const delimiter = lines[0].includes('\\t') ? '\\t' : (lines[0].includes(';') ? ';' : ',');
        headers = lines[0].split(delimiter).map(h => h.trim().replace(/^["']|["']$/g, ''));
      }

      const tVal = document.getElementById('psmTreatmentSelect')?.value;
      const oVal = document.getElementById('psmOutcomeSelect')?.value;
      const container = document.getElementById('psmCovariatesContainer');
      if (!container) return;

      const nonCovariates = new Set([tVal, oVal, 'patient_id', 'id', 'subject_id', 'ID', '_rowId']);
      const availableCovariates = headers.filter(h => !nonCovariates.has(h));

      const checkedSet = new Set();
      container.querySelectorAll('input[type="checkbox"]:checked').forEach(cb => checkedSet.add(cb.value));

      container.innerHTML = '';
      availableCovariates.forEach(cov => {
        const label = document.createElement('label');
        label.style.display = 'inline-flex';
        label.style.alignItems = 'center';
        label.style.gap = '0.3rem';
        label.style.padding = '0.2rem 0.5rem';
        label.style.background = 'var(--bg-surface)';
        label.style.border = '1px solid var(--border-subtle)';
        label.style.borderRadius = '4px';
        label.style.fontSize = '0.74rem';
        label.style.cursor = 'pointer';

        const cb = document.createElement('input');
        cb.type = 'checkbox';
        cb.value = cov;
        cb.classList.add('psm-cov-cb');
        cb.style.accentColor = 'var(--cyan-primary)';
        if (checkedSet.has(cov) || checkedSet.size === 0) {
          cb.checked = true;
        }

        cb.addEventListener('change', () => this.updatePsmCovariatesCount());

        label.appendChild(cb);
        label.appendChild(document.createTextNode(cov));
        container.appendChild(label);
      });

      this.updatePsmCovariatesCount();
    }

    updatePsmCovariatesCount() {
      const container = document.getElementById('psmCovariatesContainer');
      const countEl = document.getElementById('psmCovariatesCount');
      if (container && countEl) {
        const count = container.querySelectorAll('input[type="checkbox"]:checked').length;
        countEl.innerText = count.toString();
      }
    }

    runPsm() {
      const csvText = document.getElementById('psmCsvInput')?.value || '';
      const records = Psm.parseClinicalCsv(csvText);
      if (records.length === 0) {
        alert('Please provide valid clinical dataset rows in CSV or TSV format.');
        return;
      }

      const treatmentCol = document.getElementById('psmTreatmentSelect')?.value;
      const outcomeCol = document.getElementById('psmOutcomeSelect')?.value;

      const container = document.getElementById('psmCovariatesContainer');
      const selectedCovariates = [];
      container?.querySelectorAll('input[type="checkbox"]:checked').forEach(cb => {
        selectedCovariates.push(cb.value);
      });

      if (!treatmentCol || !outcomeCol) {
        alert('Please select both a treatment variable and an outcome variable.');
        return;
      }

      if (selectedCovariates.length === 0) {
        alert('Please select at least one confounding baseline covariate for matching.');
        return;
      }

      const caliperMultiplier = parseFloat(document.getElementById('psmCaliperMultiplier')?.value) || 0.20;
      const enforceCommonSupport = document.getElementById('psmCommonSupportCb')?.checked ?? true;

      const analysis = Psm.executeAnalysis(records, {
        treatmentCol,
        outcomeCol,
        covariateCols: selectedCovariates,
        caliperMultiplier,
        enforceCommonSupport
      });

      if (analysis.error) {
        alert(`Propensity Score Analysis Error: ${analysis.error}`);
        return;
      }

      this.psmLastAnalysis = analysis;
      if (!this.results) this.results = {};
      this.results.propensity = analysis;

      // Update Top Metrics
      const pairsEl = document.getElementById('psmMatchedPairs');
      const pairsBadge = document.getElementById('psmMatchedPairsBadge');
      if (pairsEl) pairsEl.innerText = `${analysis.nMatchedPairs} pairs`;
      if (pairsBadge) {
        const retainedPct = ((analysis.matchedTreatedN / analysis.unmatchedTreatedN) * 100).toFixed(0);
        pairsBadge.innerText = `${analysis.nMatchedPairs * 2} patients (${retainedPct}% treated matched)`;
      }

      const attEl = document.getElementById('psmAttEstimate');
      const attBadge = document.getElementById('psmAttBadge');
      const attCiEl = document.getElementById('psmAttCi');
      const attSeEl = document.getElementById('psmAttSe');
      const maxSmdEl = document.getElementById('psmMaxSmd');
      const balanceBadge = document.getElementById('psmBalanceBadge');

      const outcome = analysis.outcome;
      if (attEl) attEl.innerText = `${outcome.att >= 0 ? '+' : ''}${outcome.att.toFixed(3)}`;
      if (attBadge) {
        attBadge.innerText = outcome.isSignificant ? `p < .001 (Significant)` : `p = ${outcome.pValue.toFixed(3)}`;
        attBadge.className = outcome.isSignificant ? 'badge badge-sig' : 'badge badge-neutral';
      }

      if (attCiEl) attCiEl.innerText = `[${outcome.ci95[0].toFixed(3)}, ${outcome.ci95[1].toFixed(3)}]`;
      if (attSeEl) attSeEl.innerText = `SE = ${outcome.se.toFixed(3)}`;

      const maxSmd = analysis.balance.maxAbsSmdPost;
      if (maxSmdEl) maxSmdEl.innerText = `${maxSmd.toFixed(3)}`;
      if (balanceBadge) {
        if (maxSmd < 0.10) {
          balanceBadge.innerText = 'Excellent Balance (|SMD| < 0.10)';
          balanceBadge.className = 'badge badge-sig';
        } else if (maxSmd < 0.20) {
          balanceBadge.innerText = 'Acceptable Balance (|SMD| < 0.20)';
          balanceBadge.className = 'badge badge-warn';
        } else {
          balanceBadge.innerText = 'Residual Imbalance (|SMD| ≥ 0.20)';
          balanceBadge.className = 'badge badge-danger';
        }
      }

      // Populate Covariate Balance Table
      const bTableBody = document.getElementById('psmBalanceTableBody');
      if (bTableBody) {
        bTableBody.innerHTML = '';
        analysis.balance.balanceTable.forEach(row => {
          const tr = document.createElement('tr');
          const isBalanced = row.absSmdPost <= 0.10;
          const statusBadge = isBalanced
            ? `<span class="badge badge-sig" style="font-size: 0.70rem;">✓ Balanced (&lt; 0.10)</span>`
            : `<span class="badge badge-danger" style="font-size: 0.70rem;">⚠ Imbalance (&gt; 0.10)</span>`;

          tr.innerHTML = `
            <td style="font-weight: 600;">${row.covariate}</td>
            <td>${row.meanTreatedPre.toFixed(2)} vs ${row.meanControlPre.toFixed(2)}</td>
            <td style="color: #f43f5e; font-weight: 600;">${row.smdPre.toFixed(3)}</td>
            <td>${row.meanTreatedPost.toFixed(2)} vs ${row.meanControlPost.toFixed(2)}</td>
            <td style="color: #10b981; font-weight: 700;">${row.smdPost.toFixed(3)}</td>
            <td>${row.varRatioPost.toFixed(2)}</td>
            <td style="color: ${row.percentReduction >= 0 ? '#10b981' : '#f43f5e'}; font-weight: 600;">${row.percentReduction.toFixed(1)}%</td>
            <td>${statusBadge}</td>
          `;
          bTableBody.appendChild(tr);
        });
      }

      // Populate Outcome Comparison Table
      const oTableBody = document.getElementById('psmOutcomeTableBody');
      if (oTableBody) {
        oTableBody.innerHTML = '';
        const unadj = outcome.unadjustedDiff;
        const pUnadjStr = unadj.pValue < 0.001 ? '< .001' : unadj.pValue.toFixed(3);
        const pAttStr = outcome.pValue < 0.001 ? '< .001' : outcome.pValue.toFixed(3);

        const trUnadj = document.createElement('tr');
        trUnadj.innerHTML = `
          <td style="font-weight: 600; color: #f43f5e;">Unadjusted (Raw Observational)</td>
          <td>${unadj.diff >= 0 ? '+' : ''}${unadj.diff.toFixed(3)}</td>
          <td>${unadj.se.toFixed(3)}</td>
          <td>t = ${unadj.statistic.toFixed(2)}</td>
          <td>${pUnadjStr}</td>
          <td>[${unadj.ci95[0].toFixed(3)}, ${unadj.ci95[1].toFixed(3)}]</td>
          <td>Confounded by baseline clinical risk</td>
        `;

        const trAtt = document.createElement('tr');
        trAtt.style.background = 'rgba(16, 185, 129, 0.08)';
        trAtt.innerHTML = `
          <td style="font-weight: 700; color: #10b981;">Matched Causal ATT (1:1 NN Caliper)</td>
          <td style="font-weight: 700; color: #10b981;">${outcome.att >= 0 ? '+' : ''}${outcome.att.toFixed(3)}</td>
          <td>${outcome.se.toFixed(3)}</td>
          <td>${outcome.type === 'continuous' ? `t = ${outcome.statistic.toFixed(2)}` : `z = ${outcome.statistic.toFixed(2)}`}</td>
          <td style="font-weight: 700;">${pAttStr}</td>
          <td style="font-weight: 600;">[${outcome.ci95[0].toFixed(3)}, ${outcome.ci95[1].toFixed(3)}]</td>
          <td>Unconfounded treatment effect on treated</td>
        `;

        oTableBody.appendChild(trUnadj);
        oTableBody.appendChild(trAtt);
      }

      this.renderPsmPlots();
      this.updatePsmScript();

      const reportEl = document.getElementById('psmReportText');
      if (reportEl) {
        reportEl.innerText = analysis.reportText;
      }
    }

    renderPsmPlots() {
      if (!this.psmLastAnalysis) return;
      const analysis = this.psmLastAnalysis;

      if (this.psmActiveView === 'love') {
        const engine = this.engines['psmLovePlotCanvas'];
        if (engine && Plots.renderLovePlot) {
          Plots.renderLovePlot(engine, analysis.balance, {
            title: `Love Plot: Baseline Covariate Balance (N = ${analysis.nMatchedPairs} Pairs)`
          });
        }
      } else {
        const engine = this.engines['psmOverlapCanvas'];
        if (engine && Plots.renderPsmOverlapPlot) {
          Plots.renderPsmOverlapPlot(engine, analysis.overlap, {
            title: `Propensity Score Distribution & Common Support Overlap (Caliper = ${analysis.caliper.width.toFixed(4)})`
          });
        }
      }
    }

    switchPsmScript(lang) {
      this.psmCurrentScriptLang = lang;
      const pyBtn = document.getElementById('psmScriptLangPy');
      const rBtn = document.getElementById('psmScriptLangR');
      const stataBtn = document.getElementById('psmScriptLangStata');

      if (pyBtn) {
        pyBtn.className = lang === 'python' ? 'btn btn-primary btn-sm' : 'btn btn-secondary btn-sm';
      }
      if (rBtn) {
        rBtn.className = lang === 'r' ? 'btn btn-primary btn-sm' : 'btn btn-secondary btn-sm';
      }
      if (stataBtn) {
        stataBtn.className = lang === 'stata' ? 'btn btn-primary btn-sm' : 'btn btn-secondary btn-sm';
      }

      this.updatePsmScript();
    }

    updatePsmScript() {
      if (!this.psmLastAnalysis) return;
      const analysis = this.psmLastAnalysis;
      const codeBlock = document.getElementById('psmCodeBlock');
      if (!codeBlock) return;

      const opts = {
        treatmentCol: analysis.treatmentCol,
        outcomeCol: analysis.outcomeCol,
        covariateCols: analysis.covariateCols,
        caliperMultiplier: analysis.caliper.multiplier,
        enforceCommonSupport: analysis.caliper.enforceCommonSupport,
        isBinaryOutcome: analysis.outcome.type === 'binary'
      };

      let code = '';
      if (this.psmCurrentScriptLang === 'python') {
        code = Psm.generatePythonScript(opts);
      } else if (this.psmCurrentScriptLang === 'r') {
        code = Psm.generateRScript(opts);
      } else if (this.psmCurrentScriptLang === 'stata') {
        code = Psm.generateStataScript(opts);
      }

      codeBlock.innerText = code;
    }

    copyPsmScript() {
      const codeBlock = document.getElementById('psmCodeBlock');
      const copyBtn = document.getElementById('psmCopyScriptBtn');
      if (!codeBlock) return;

      const text = codeBlock.innerText;
      if (navigator.clipboard && navigator.clipboard.writeText) {
        navigator.clipboard.writeText(text).then(() => {
          if (copyBtn) {
            const orig = copyBtn.innerText;
            copyBtn.innerText = '✓ Copied!';
            setTimeout(() => { copyBtn.innerText = orig; }, 1800);
          }
        });
      } else {
        const ta = document.createElement('textarea');
        ta.value = text;
        ta.style.position = 'fixed';
        ta.style.opacity = '0';
        document.body.appendChild(ta);
        ta.select();
        document.execCommand('copy');
        document.body.removeChild(ta);
        if (copyBtn) {
          const orig = copyBtn.innerText;
          copyBtn.innerText = '✓ Copied!';
          setTimeout(() => { copyBtn.innerText = orig; }, 1800);
        }
      }
    }

    downloadPsmScript() {
      const codeBlock = document.getElementById('psmCodeBlock');
      if (!codeBlock) return;
      const text = codeBlock.innerText;

      const extMap = { python: 'py', r: 'R', stata: 'do' };
      const ext = extMap[this.psmCurrentScriptLang] || 'txt';
      const filename = `psm_clinical_analysis.${ext}`;

      const blob = new Blob([text], { type: 'text/plain;charset=utf-8' });
      const url = URL.createObjectURL(blob);
      const link = document.createElement('a');
      link.href = url;
      link.download = filename;
      document.body.appendChild(link);
      link.click();
      document.body.removeChild(link);
      setTimeout(() => URL.revokeObjectURL(url), 1000);
    }

}


