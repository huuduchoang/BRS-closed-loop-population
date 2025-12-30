% compare_max_vs_softplus_overlay.m
% Runs closedloop_populationM_autapse_ReLU twice then overlays MAX vs SOFTPLUS
% on the same axes for each variable.

clear; clc; close all;

%% Number of neurons and heterogeneity settings
N = 1;

% Heterogeneity (set sigmas to 0 for no variability)
Eleak_mu    = -65;   Eleak_sigma    = 0;
gNaP_mu     =   2.8; gNaP_sigma     = 0;
gsyn_mu     =   0;   gsyn_sigma     = 0;
phi_mu      =   0.3; phi_sigma      = 0;
thetaO2_mu  =  85;   thetaO2_sigma  = 0;
sigmaO2_mu  =  30;   sigmaO2_sigma  = 0;

% Build the params struct (all draws happen inside)
params = setup_params( ...
    N, ...
    Eleak_mu,    Eleak_sigma, ...
    gNaP_mu,     gNaP_sigma, ...
    gsyn_mu,     gsyn_sigma, ...
    phi_mu,      phi_sigma, ...
    thetaO2_mu,  thetaO2_sigma, ...
    sigmaO2_mu,  sigmaO2_sigma);

%% Fix M to be 8e-6
params.Mfun = [];
params.M    = 1.2e-5;
params.W    = 0;

%% Softplus steepness (choose this)
r_soft = 80000;

%% Initial conditions (same inits for each cell)
% [V n h alpha vollung PO2lung PO2blood]
initsA = [-58.5754 0.0006 0.7252 0.0010 2.2665 103.3461 102.2224];

v0       = repmat(initsA(1), N, 1);
n0       = repmat(initsA(2), N, 1);
h0       = repmat(initsA(3), N, 1);
s0       = zeros(N,1);            % synaptic gates
alpha0   = initsA(4);
voll0    = initsA(5);
PO2lung0 = initsA(6);
PO2blood0= initsA(7);

u0 = [v0; n0; h0; s0; alpha0; voll0; PO2lung0; PO2blood0];

%% Integrate with ode15s
tf   = 120000;   % ms
opts = odeset('RelTol',1e-13,'AbsTol',1e-10);

%% -------------------- RUN 1: MAX (ReLU) --------------------
params_max = params;   % IMPORTANT: do NOT set params_max.softplus_r
[t_max, U_max] = ode15s(@(t,u) closedloop_populationM_autapse_ReLU(t,u,params_max), ...
                        [0 tf], u0, opts);
time_s_max = t_max/1000;

%% -------------------- RUN 2: SOFTPLUS --------------------
params_sp = params;
params_sp.softplus_r = r_soft;
[t_sp, U_sp] = ode15s(@(t,u) closedloop_populationM_autapse_ReLU(t,u,params_sp), ...
                      [0 tf], u0, opts);
time_s_sp = t_sp/1000;

%% Extract per-neuron and loop signals (MAX)
V_max    = U_max(:,    1:N);
n_max    = U_max(:,N+1:2*N);
h_max    = U_max(:,2*N+1:3*N);
s_max    = U_max(:,3*N+1:4*N);
alpha_max    = U_max(:,4*N+1);
vol_lung_max = U_max(:,4*N+2);
PO2_lung_max = U_max(:,4*N+3);
PO2_b_max    = U_max(:,4*N+4);

%% Extract per-neuron and loop signals (SOFTPLUS)
V_sp    = U_sp(:,    1:N);
n_sp    = U_sp(:,N+1:2*N);
h_sp    = U_sp(:,2*N+1:3*N);
s_sp    = U_sp(:,3*N+1:4*N);
alpha_sp    = U_sp(:,4*N+1);
vol_lung_sp = U_sp(:,4*N+2);
PO2_lung_sp = U_sp(:,4*N+3);
PO2_b_sp    = U_sp(:,4*N+4);

%% Compute g_tonic for all neurons (MAX + SOFTPLUS)
phi_vec = params.phi;
th_vec  = params.thetaO2;
sg_vec  = params.sigmaO2;

gtonic_max = zeros(length(t_max),N);
gtonic_sp  = zeros(length(t_sp),N);
for i = 1:N
  gtonic_max(:,i) = phi_vec(i).*(1 - tanh((PO2_b_max - th_vec(i))./sg_vec(i)));
  gtonic_sp(:,i)  = phi_vec(i).*(1 - tanh((PO2_b_sp  - th_vec(i))./sg_vec(i)));
end

%% Colors for neuron overlays
if exist('distinguishable_colors','file')
  colors = distinguishable_colors(N);
else
  colors = lines(N);
end

%% -------------------- Figure 1: OVERLAY 2×3 closed-loop grid --------------------
figure('Position',[100 100 900 600],'Color','w');
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

% Panel 1: CPG (V)
nexttile(1), hold on
for i=1:N
  plot(time_s_max, V_max(:,i), 'Color', colors(i,:), 'LineWidth', 1.3, 'LineStyle','-');
  plot(time_s_sp,  V_sp(:,i),  'Color', colors(i,:), 'LineWidth', 1.3, 'LineStyle','--');
end
hold off
ylabel('V_i (mV)');
title(sprintf('CPG (V) | max (solid) vs softplus (dashed), r=%g', r_soft));
xlim([0 tf/1000]); box on; grid on

% Panel 2: Motor pool (α)
nexttile(2), hold on
plot(time_s_max, alpha_max, 'k', 'LineWidth', 1.5, 'LineStyle','-');
plot(time_s_sp,  alpha_sp,  'k', 'LineWidth', 1.5, 'LineStyle','--');
hold off
ylabel('\alpha'); title('Motor pool');
xlim([0 tf/1000]); box on; grid on

% Panel 3: Lung volume
nexttile(3), hold on
plot(time_s_max, vol_lung_max, 'k', 'LineWidth', 1.5, 'LineStyle','-');
plot(time_s_sp,  vol_lung_sp,  'k', 'LineWidth', 1.5, 'LineStyle','--');
hold off
ylabel('Vol_{lung}'); title('Lung Volume');
xlim([0 tf/1000]); box on; grid on

% Panel 4: Chemosensation (g_tonic)
nexttile(4), hold on
for i=1:N
  plot(time_s_max, gtonic_max(:,i), 'Color', colors(i,:), 'LineWidth', 1.3, 'LineStyle','-');
  plot(time_s_sp,  gtonic_sp(:,i),  'Color', colors(i,:), 'LineWidth', 1.3, 'LineStyle','--');
end
hold off
ylabel('g_{tonic,i} (nS)'); title('Chemosensation');
xlim([0 tf/1000]); box on; grid on

% Panel 5: Blood O2
nexttile(5), hold on
plot(time_s_max, PO2_b_max, 'k', 'LineWidth', 1.5, 'LineStyle','-');
plot(time_s_sp,  PO2_b_sp,  'k', 'LineWidth', 1.5, 'LineStyle','--');
hold off
ylabel('P_{aO_2} (mmHg)'); title('Blood Oxygen');
xlim([0 tf/1000]); box on; grid on

% Panel 6: Lung O2
nexttile(6), hold on
plot(time_s_max, PO2_lung_max, 'k', 'LineWidth', 1.5, 'LineStyle','-');
plot(time_s_sp,  PO2_lung_sp,  'k', 'LineWidth', 1.5, 'LineStyle','--');
hold off
ylabel('P_{O_2}^{lung} (mmHg)'); title('Lung Oxygen');
xlabel('Time (s)'); xlim([0 tf/1000]); box on; grid on

% Legend (put it once for the whole figure)
legend({'max (solid)','softplus (dashed)'}, 'Location','best');

%% -------------------- Figure 2: OVERLAY gating h and n --------------------
figure('Position',[120 120 600 400],'Color','w');
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

% Row 1: h_i
nexttile, hold on
for i=1:N
  plot(time_s_max, h_max(:,i), 'Color', colors(i,:), 'LineWidth', 1.3, 'LineStyle','-');
  plot(time_s_sp,  h_sp(:,i),  'Color', colors(i,:), 'LineWidth', 1.3, 'LineStyle','--');
end
hold off
ylabel('h_i'); title(sprintf('CPG inactivation (h) | max vs softplus r=%g', r_soft));
xlim([0 tf/1000]); box on; grid on

% Row 2: n_i
nexttile, hold on
for i=1:N
  plot(time_s_max, n_max(:,i), 'Color', colors(i,:), 'LineWidth', 1.3, 'LineStyle','-');
  plot(time_s_sp,  n_sp(:,i),  'Color', colors(i,:), 'LineWidth', 1.3, 'LineStyle','--');
end
hold off
ylabel('n_i'); title('CPG activation (n)');
xlabel('Time (s)'); xlim([0 tf/1000]); box on; grid on
