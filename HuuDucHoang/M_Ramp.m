clear; clc;
%% Self coupling values to compare
Wvals = [0, 1.0];
colors = [0 0 0; 0 0.447 0.741];

%% ----- Model + params (single neuron with autapse) ----------------------
N = 1;
params_base = setup_params( ...
    N, ...
    -65, 0, ...             % E_leak mean, sigma
    2.8, 0, ...             % g_NaP mean, sigma
    0.1, 0, ...             % g_syn mean, sigma (unused for N=1 if autapse via W)
    0.3, 0, ...             % phi mean, sigma
    85,  0, ...             % theta_O2 mean, sigma
    30,  0);                % sigma_O2 mean, sigma
params_base.W = 0;          % will overwrite per run

%% ----- Warm-up and ramp design -----------------------------------------
M0      = 8e-6;             % baseline metabolic demand
M1      = 2.3e-5;           % final demand (past the cliff)
T_warm  = 180e3;            % ms, warm-up at M0
T_ramp  = 240e3;            % ms, linear ramp duration
T_tail  =  60e3;            % ms, post-ramp observation
t_end   = T_warm + T_ramp + T_tail;

% M(t): piecewise linear ramp
Mfun = @(t) ...
    (t < T_warm).*M0 + ...
    (t >= T_warm & t < T_warm+T_ramp).*(M0 + (M1-M0).*(t-T_warm)/T_ramp) + ...
    (t >= T_warm+T_ramp).*M1;

%% ----- Initial conditions ----------------------------------------------
initsA   = [-58.5754 0.0006 0.7252 0.0010 2.2665 103.3461 102.2224];
v0       = repmat(initsA(1), N, 1);
n0       = repmat(initsA(2), N, 1);
h0       = repmat(initsA(3), N, 1);
s0       = zeros(N,1);
alpha0   = initsA(4);
voll0    = initsA(5);
PO2lung0 = initsA(6);
PO2blood0= initsA(7);
u0       = [v0; n0; h0; s0; alpha0; voll0; PO2lung0; PO2blood0];

%% ----- Integrate for each W and store results --------------------------
opts = odeset('RelTol',1e-8,'AbsTol',1e-8);

T_all   = cell(numel(Wvals),1);
U_all   = cell(numel(Wvals),1);
Mt_all  = cell(numel(Wvals),1);
dPO2_all= cell(numel(Wvals),1);

for k = 1:numel(Wvals)
    p        = params_base;
    p.W      = Wvals(k);
    p.Mfun   = Mfun;

    [t,U]    = ode15s(@(tt,u) closedloop_populationM_autapse(tt,u,p), [0 t_end], u0, opts);
    ts       = t/1000;                       % seconds
    M_t      = Mfun(t);

    % Extract
    V_all    = U(:,1:N);
    alpha    = U(:,4*N+1);
    vol_lung = U(:,4*N+2);
    PO2_lung = U(:,4*N+3);
    PO2_b    = U(:,4*N+4);

    % Derivative of PaO2
    dpo2     = gradient(PO2_b, mean(diff(ts)));

    % Stash
    T_all{k}    = ts;
    U_all{k}    = struct('V',V_all,'alpha',alpha,'vol',vol_lung,'PO2b',PO2_b);
    Mt_all{k}   = M_t;
    dPO2_all{k} = dpo2;
end

%% ----- Plots: superimposed traces for both W ---------------------------
figure('Color','w','Position',[100 100 1050 750]);

% Panel 1 — mean V (for N=1, it's just V)
ax1 = subplot(3,1,1); hold on; box on; grid on;
hV = gobjects(1,numel(Wvals));
for k = 1:numel(Wvals)
    ts    = T_all{k};
    Vmean = mean(U_all{k}.V,2);
    hV(k) = plot(ts, Vmean, '-', 'LineWidth', 1.4, 'Color', colors(k,:));
end
ylabel('\langle V\rangle (mV)');
title('Burst structure during a slow ramp of M');
xlim([0, T_all{end}(end)]);  % use last run’s span
legend(hV, arrayfun(@(w) sprintf('W=%.2f',w), Wvals, 'uni', 0), ...
       'Location','best'); legend boxoff;

% Panel 2 — motor pool (left) and lung volume (right), both runs
ax2 = subplot(3,1,2); hold on; box on; grid on;
yyaxis left
hA = gobjects(1,numel(Wvals));
for k = 1:numel(Wvals)
    ts = T_all{k};
    hA(k) = plot(ts, U_all{k}.alpha, '-', 'LineWidth', 1.2, 'Color', colors(k,:));
end
ylabel('\alpha')

yyaxis right
hVol = gobjects(1,numel(Wvals));
for k = 1:numel(Wvals)
    ts = T_all{k};
    hVol(k) = plot(ts, U_all{k}.vol, '--', 'LineWidth', 1.2, 'Color', colors(k,:));
end
ylabel('Vol_{lung}')
xlim([0, T_all{end}(end)]);
% Legends on each axis
yyaxis left
legend(hA, arrayfun(@(w) sprintf('\\alpha, W=%.2f',w), Wvals, 'uni', 0), ...
       'Location','northwest'); legend boxoff;
yyaxis right
legend(hVol, arrayfun(@(w) sprintf('vol_{lung}, W=%.2f',w), Wvals, 'uni', 0), ...
       'Location','northeast'); legend boxoff;

% Panel 3 — PaO2 (left) with M(t) (right axis, single ramp shown once)
ax3 = subplot(3,1,3); hold on; box on; grid on;
hP = gobjects(1,numel(Wvals));
for k = 1:numel(Wvals)
    ts = T_all{k};
    hP(k) = plot(ts, U_all{k}.PO2b, '-', 'LineWidth', 1.4, 'Color', colors(k,:));
end
ylabel('P_{aO_2} (mmHg)'); xlabel('Time (s)');
xlim([0, T_all{end}(end)]);
yyaxis right
plot(T_all{1}, Mt_all{1}*1e5, 'k:', 'LineWidth', 1.5);  % one ramp for reference
ylabel('M \times 10^{-5}');
yyaxis left
legend(hP, arrayfun(@(w) sprintf('P_{aO_2}, W=%.2f',w), Wvals, 'uni', 0), ...
       'Location','best'); legend boxoff;

linkaxes([ax1 ax2 ax3],'x');

%% ----- Optional: dPaO2/dt overlay in a separate figure -----------------
figure('Color','w','Position',[100 100 950 400]);
hold on; box on; grid on;
for k = 1:numel(Wvals)
    plot(T_all{k}, dPO2_all{k}, 'LineWidth', 1.2, 'Color', colors(k,:));
end
yline(0,'k:');
xlabel('Time (s)'); ylabel('dP_{a}O_2/dt (mmHg/s)');
title('Instantaneous PaO_2 slope'); 
legend(arrayfun(@(w) sprintf('W=%.2f',w), Wvals, 'uni', 0), 'Location','best'); legend boxoff;

%% ===== Cycle-aware alignment using vol_lung onsets (instead of alpha) ===
% Align W = Wvals(2) (target) onto W = Wvals(1) (reference)

refIdx    = 1;              % reference run: W = Wvals(1)
tgtIdx    = 2;              % target run to be shifted: W = Wvals(2)
mCand     = -6:6;           % integer cycle shifts to test
minPairs  = 5;              % need at least this many matched breaths
minGapSec = 1.0;            % debounce for onset detection (seconds)
fracThr   = 0.40;           % threshold = min + frac*(range) for upcrossings

% --- pull the signals ---
tsA = T_all{refIdx};  vL_A = U_all{refIdx}.vol;   % reference vol_lung
tsB = T_all{tgtIdx};  vL_B = U_all{tgtIdx}.vol;   % target   vol_lung

% --- detect breath onsets from vol_lung upcrossings with debounce ---
thrA = min(vL_A) + fracThr*(max(vL_A) - min(vL_A));
thrB = min(vL_B) + fracThr*(max(vL_B) - min(vL_B));
tOnA = detect_onsets_threshold(tsA, vL_A, thrA, minGapSec);  % seconds
tOnB = detect_onsets_threshold(tsB, vL_B, thrB, minGapSec);  % seconds

% --- choose best integer-cycle shift m and time shift tau* --------------
[m_best, tau_best, pairs_used] = align_onsets(tOnA, tOnB, mCand, minPairs);
fprintf('vol_lung alignment: m_best = %+d, tau* = %.3f s, matched pairs = %d\n', ...
        m_best, tau_best, size(pairs_used,1));

% --- shift target timeline for re-plotting against reference ------------
tsB_shift = tsB - tau_best;

%% Quick sanity plot: vol_lung onsets before/after alignment
figure('Color','w','Position',[60 60 980 360]);
subplot(1,2,1); hold on; box on; grid on;
plot(tsA, vL_A, 'k-', 'LineWidth', 1.0);
plot(tsB, vL_B, 'Color',[0 0.447 0.741], 'LineWidth', 1.0);
plot(tOnA, interp1(tsA, vL_A, tOnA), 'ko',  'MarkerFaceColor','k', 'MarkerSize',4);
plot(tOnB, interp1(tsB, vL_B, tOnB), 'o',   'Color',[0 0.447 0.741], 'MarkerSize',4);
title('vol_{lung} onsets (before alignment)'); xlabel('Time (s)'); ylabel('vol_{lung}');

subplot(1,2,2); hold on; box on; grid on;
plot(tsA, vL_A, 'k-', 'LineWidth', 1.0);
plot(tsB_shift, vL_B, 'Color',[0 0.447 0.741], 'LineWidth', 1.0);
plot(tOnA,            interp1(tsA, vL_A, tOnA), 'ko',  'MarkerFaceColor','k', 'MarkerSize',4);
plot(tOnB - tau_best, interp1(tsB, vL_B, tOnB), 'o',   'Color',[0 0.447 0.741], 'MarkerSize',4);
title(sprintf('vol_{lung} onsets (aligned: m=%+d, \\tau^*=%.3fs)', m_best, tau_best));
xlabel('Time (s)'); ylabel('vol_{lung}');

%% Re-plot aligned traces (V, vol_lung, PaO2) with M(t) on right axis
colRef = [0 0 0];
colTgt = [0 0.447 0.741];

figure('Color','w','Position',[100 100 1050 750]);

% Panel 1 — V (for N=1, mean(V)=V)
ax1 = subplot(3,1,1); hold on; box on; grid on;
plot(tsA, mean(U_all{refIdx}.V,2), '-', 'LineWidth', 1.4, 'Color', colRef);
plot(tsB_shift, mean(U_all{tgtIdx}.V,2), '-', 'LineWidth', 1.4, 'Color', colTgt);
ylabel('\langle V\rangle (mV)');
title(sprintf('Aligned by vol_{lung} (m=%+d, \\tau^*=%.3fs)', m_best, tau_best));
xlim([max([tsA(1), tsB_shift(1)]), min([tsA(end), tsB_shift(end)])]);
legend({sprintf('W=%.2f', Wvals(1)), sprintf('W=%.2f (shifted)', Wvals(2))}, ...
       'Location','best'); legend boxoff;

% Panel 2 — vol_lung (left) and P_{O2}^{lung} (right), both runs
ax2 = subplot(3,1,2); hold on; box on; grid on;
yyaxis left
plot(tsA, U_all{refIdx}.vol, '-', 'LineWidth', 1.2, 'Color', colRef);
plot(tsB_shift, U_all{tgtIdx}.vol, '-', 'LineWidth', 1.2, 'Color', colTgt);
ylabel('vol_{lung}');
yyaxis right
xlim([max([tsA(1), tsB_shift(1)]), min([tsA(end), tsB_shift(end)])]);

% Panel 3 — PaO2 with M(t) (right axis)
ax3 = subplot(3,1,3); hold on; box on; grid on;
plot(tsA, U_all{refIdx}.PO2b, '-', 'LineWidth', 1.4, 'Color', colRef);
plot(tsB_shift, U_all{tgtIdx}.PO2b, '-', 'LineWidth', 1.4, 'Color', colTgt);
ylabel('P_{aO_2} (mmHg)'); xlabel('Time (s)');
xlim([max([tsA(1), tsB_shift(1)]), min([tsA(end), tsB_shift(end)])]);
yyaxis right
plot(T_all{refIdx}, Mt_all{refIdx}*1e5, 'k:', 'LineWidth', 1.5);
ylabel('M \times 10^{-5}');
yyaxis left
legend({sprintf('W=%.2f', Wvals(1)), sprintf('W=%.2f (shifted)', Wvals(2))}, ...
       'Location','best'); legend boxoff;

linkaxes([ax1 ax2 ax3],'x');

%% ===== Helper functions ===
function tOn = detect_onsets_threshold(ts, sig, thr, minGapSec)
% Upward threshold crossings in a smooth signal (vol_lung).
% Debounce to keep only the first onset per breath.
above = sig > thr;
iUp   = find(diff(above) == 1) + 1;     % index immediately after up-crossing
tCand = ts(iUp);
tOn   = [];
for k = 1:numel(tCand)
    if isempty(tOn) || (tCand(k) - tOn(end) >= minGapSec)
        tOn(end+1) = tCand(k); %#ok<AGROW>
    end
end
end

function [m_best, tau_best, pairs_used] = align_onsets(tA, tB, mCand, minPairs)
% Try integer cycle shifts m, pair tB_k with tA_{k+m}, pick m that makes
% the timing differences most consistent (min IQR). tau* = median of diffs.
bestObj    = Inf;
m_best     = 0;
tau_best   = 0;
pairs_used = [];

for m = mCand
    kA = (1+m):(numel(tA)+m);   % target indices mapped onto reference
    kB = 1:numel(tB);
    mask = (kA>=1) & (kA<=numel(tA));
    kA = kA(mask);  kB = kB(mask);
    if numel(kA) < minPairs, continue; end

    d   = tB(kB) - tA(kA);      % time differences for paired breaths
    obj = iqr(d);               % robust spread (could also use var)
    if obj < bestObj
        bestObj    = obj;
        m_best     = m;
        tau_best   = median(d);
        pairs_used = [kA(:) kB(:)];
    end
end
end

%% Figure showing one way of shifting the time and PaO2 levels 
% so that the PaO2 curves match as closely as possible near the collapse
% point.

figure
subplot(3,1,1)
plot(tsA,U_all{1}.V)
hold on
plot(tsB-63,U_all{2}.V)
grid on
ylabel('Voltage (mV)')
set(gca,'FontSize',20)

subplot(3,1,2)
plot(tsA,vL_A)
hold on
plot(tsB-63,vL_B)
grid on
ylabel('Lung Volume (L)')
set(gca,'FontSize',20)

subplot(3,1,3)
plot(tsA,U_all{1}.PO2b)
hold on
plot(tsB-63,U_all{2}.PO2b+2)
grid on
xlabel('Time (sec)')
ylabel('PaO2 (mm Hg)')
set(gca,'FontSize',20)
