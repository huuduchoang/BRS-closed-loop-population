% test_N_neurons.m
%% Number of neurons and heterogeneity settings
N = 20; 
% seed = 1; % Stable Alternans (19 vs 9) (First Mode)
% seed = 2; % Normal Eupneic (constant 18) (Second Mode)
% seed = 3; % Normal Eupneic (constant 19, run for 200ms) (Second Mode)
% seed = 4; % Stable Alternans (17 vs 16) (First Mode)
% seed = 5; % Stable Alternans (20 vs 8) (First Mode)
% seed = 6; % Unstable/Noisy Alternans (19 vs 11/12/13)
% seed = 7; % Unstable/Noisy Alternans (Alternate 7~9 to 19)
% seed = 8; % Normal Eupneic (constant 19)
% seed = 9; % Stable Alternans (19 vs 4)
% seed = 10; % Nonlinear Response, Heart Failure? (Mode 3)
% seed = 11; % Unstable Alternans (17 vs 9)
% seed = 12; % Unstable Alternans (20 vs 11)
% seed = 13; % Noisy Alternans (Alternate 9 to 14~15) --> Recheck at longer
% runtime
% seed = 14; % Stable Alternans (19 vs 18)
% seed = 15; % Normal Eupneic (constant 19, occasional 20)
% seed = 16; % Noisy Alternans  (20 vs 9, occasional 8)
% seed = 17; % Normal Eupneic (constant 17)
% seed = 18; % Noisy Alternans (Alternate 13~14 to 15~16)
% seed = 19; % Stable Alternans (19 vs 3) 
seed = 20; % Nonlinear/Irregular (most of peaks and troughs of 15 vs 12)

rng(seed, 'twister');

% Heterogeneity (you can set any sigma to zero if you want no variability)

% Leads to period-halving pattern
% Reduce to 5% max --> Standard Period
% Eleak_mu    = -65;   Eleak_sigma    = 2;
% gNaP_mu     =   2.8; gNaP_sigma     = 0.1;
% gsyn_mu     =   0.1;   gsyn_sigma     = 0.005;
% phi_mu      =   0.3; phi_sigma      = 0.003;
% thetaO2_mu  =  85;   thetaO2_sigma  = 8.5;
% sigmaO2_mu  =  30;   sigmaO2_sigma  = 3;

% Eleak_mu    = -65;   Eleak_sigma    = 0;
% gNaP_mu     =   2.8; gNaP_sigma     = 0;
% gsyn_mu     =   0.1;   gsyn_sigma     = 0;
% phi_mu      =   0.3; phi_sigma      = 0;
% thetaO2_mu  =  85;   thetaO2_sigma  = 0;
% sigmaO2_mu  =  30;   sigmaO2_sigma  = 0;

Eleak_mu    = -65;   Eleak_sigma    = 2;
gNaP_mu     =   2.8; gNaP_sigma     = 0.1;
gsyn_mu     =   0.1;   gsyn_sigma     = 0.003;
phi_mu      =   0.3; phi_sigma      = 0.009;
thetaO2_mu  =  85;   thetaO2_sigma  = 3;
sigmaO2_mu  =  30;   sigmaO2_sigma  = 1;

% Build the params struct (all draws happen inside)
params = setup_params( ...
    N, ...
    Eleak_mu,    Eleak_sigma, ...
    gNaP_mu,     gNaP_sigma, ...
    gsyn_mu,     gsyn_sigma, ...
    phi_mu,      phi_sigma, ...
    thetaO2_mu,  thetaO2_sigma, ...
    sigmaO2_mu,  sigmaO2_sigma);

save( sprintf('saved_params_seed%d.mat',seed), 'params' );

%% Initial conditions (use the same inits for each cell)
% Single-cell inits from Diekman et al. 2017 (panel A):
% -41.7429 0.0313 0.3442 0.0025 2.4355 23.9533 23.3940
% -58.5754 0.0006 0.7252 0.0010 2.2665 103.3461 102.222
% -61.9818 0.0006 0.7252 0.0001157 2.02995 116.445 114.901 (eupneic, gsyn =
% 0.5, sigma = 0.1)
% -60.3782 0.0006 0.667918 0.00136 2.3223 99.1582 98.0934
%[V n h alpha vollung PO2lung PO2blood]
initsA = [-60.3782 0.0006 0.667918 0.00136 2.3223 99.1582 98.0934];

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
tf      = 500000;   % ms
opts    = odeset('RelTol',1e-9,'AbsTol',1e-9);
[t,U]   = ode15s(@(t,u) closedloop_population(t,u,params), [0 tf], u0, opts);
time_s  = t/1000;

%% 4) Extract all per‐neuron and loop signals
V_all    = U(:,    1:N);         % t×N
h_all    = U(:,2*N+1:3*N);       % t×N (for fig 2)
alpha    = U(:,4*N+1);           % t×1
vol_lung = U(:,4*N+2);           % t×1
PO2_lung = U(:,4*N+3);           % t×1
PO2_b    = U(:,4*N+4);           % t×1

[burstTimes, Nspike] = track_spikes_per_burst(time_s, PO2_b, V_all);

% compute g_tonic for all neurons
phi_vec = params.phi;
th_vec  = params.thetaO2;
sg_vec  = params.sigmaO2;
gtonic_all = zeros(length(t),N);
for i=1:N
  gtonic_all(:,i) = ...
    phi_vec(i).*(1 - tanh((PO2_b - th_vec(i))./sg_vec(i)));
end

%% 5) Figure 1: 2×3 closed‐loop grid, overlaying all N where appropriate
figure('Position',[100 100 900 600],'Color','w');
tl = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
colors = lines(N);

% (1,1) CPG: all V_i
nexttile(1)
hold on
for i=1:N
  plot(time_s, V_all(:,i),'Color',colors(i,:),'LineWidth',1.2);
end
hold off
ylabel('V_i (mV)'), title('CPG (V)'), xlim([0 tf/1000])
box off, grid on

% (1,2) Motor pool: α (single)
nexttile(2)
plot(time_s, alpha,'k','LineWidth',1.5)
ylabel('\alpha'), title('Motor pool'), xlim([0 tf/1000])
box off, grid on

% (1,3) Lung volume (single)
nexttile(3)
plot(time_s, vol_lung,'k','LineWidth',1.5)
ylabel('Vol_{lung}'), title('Lung Volume'), xlim([0 tf/1000])
box off, grid on

% (2,1) Chemosensation: all g_tonic,i
nexttile(4)
hold on
for i=1:N
  plot(time_s, gtonic_all(:,i),'Color',colors(i,:),'LineWidth',1.2);
end
hold off
ylabel('g_{tonic,i} (nS)'), title('Chemosensation'), xlim([0 tf/1000])
box off, grid on

% (2,2) Blood O2 (single)
nexttile(5)
plot(time_s, PO2_b,'k','LineWidth',1.5)
ylabel('P_{aO_2} (mmHg)'), title('Blood Oxygen'), xlim([0 tf/1000])
box off, grid on

% (2,3) Lung O2 (single)
nexttile(6)
plot(time_s, PO2_lung,'k','LineWidth',1.5)
ylabel('P_{O_2}^{lung} (mmHg)'), title('Lung Oxygen')
xlabel('Time (s)'), xlim([0 tf/1000])
box off, grid on

%% 6) Figure 2: CPG gating h and n for all neurons
figure('Position',[100 100 600 400],'Color','w');
tl2 = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

% Row 1: h_i for all i
nexttile
hold on
for i=1:N
  plot(time_s, h_all(:,i),'Color',colors(i,:),'LineWidth',1.2);
end
hold off
ylabel('h_i'), title('CPG inactivation (h)'), xlim([0 tf/1000])
box off, grid on

% Row 2: n_i for all i
n_all = U(:, N+1:2*N);
nexttile
hold on
for i=1:N
  plot(time_s, n_all(:,i),'Color',colors(i,:),'LineWidth',1.2);
end
hold off
ylabel('n_i'), title('CPG activation (n)')
xlabel('Time (s)'), xlim([0 tf/1000])
box off, grid on