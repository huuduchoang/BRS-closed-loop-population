%% parallel_find_gsyn_threshold.m
% Parallelized sweep of g_syn to find network period vs. coupling

%% 0) Launch parallel pool (once per MATLAB session)
if isempty(gcp('nocreate'))
    parpool('local');
end

%% 1) Network & simulation settings
N    = 5;                % number of neurons
tf   = 200000;            % simulation length (ms)
opts = odeset('RelTol',1e-6,'AbsTol',1e-9);

% Fixed heterogeneity (no scatter)
Eleak_mu      = -65;   Eleak_sigma    =  0;
gNaP_mu       =   2.8; gNaP_sigma     =  0;
phi_mu        =   0.3; phi_sigma      =  0;
thetaO2_mu    =    85; thetaO2_sigma  =  0;
sigmaO2_mu    =    30; sigmaO2_sigma  =  0;

% Sweep range for mean synaptic coupling
gsyn_vals = linspace(0,0.2,21);    % 21 points from 0 to 0.2 nS
periods   = nan(size(gsyn_vals));  % preallocate output

%% 2) Build explicit initial state vector u0
initsA    = [-56.8172, 9.5344e-04, 0.7454, 2.0026e-04, ...
              2.0525, 98.9638, 97.7927];
v0        = repmat(initsA(1), N, 1);
n0        = repmat(initsA(2), N, 1);
h0        = repmat(initsA(3), N, 1);
s0        = zeros(N,1);
alpha0    = initsA(4);
voll0     = initsA(5);
PO2lung0  = initsA(6);
PO2blood0 = initsA(7);
u0 = [ v0; n0; h0; s0; alpha0; voll0; PO2lung0; PO2blood0 ];

%% 3) Parallel sweep over gsyn_vals
parfor k = 1:numel(gsyn_vals)
    % 3a) Build params for this coupling
    params = setup_params( ...
      N, ...
      Eleak_mu,      Eleak_sigma, ...
      gNaP_mu,       gNaP_sigma, ...
      gsyn_vals(k),  0, ...        % gsyn mean + zero sigma
      phi_mu,        phi_sigma, ...
      thetaO2_mu,    thetaO2_sigma, ...
      sigmaO2_mu,    sigmaO2_sigma);
    % clamp any negatives
    params.gsyn(params.gsyn < 0) = 0;

    % 3b) Integrate the network
    [t,U] = ode15s(@(t,u) closedloop_population(t,u,params), [0 tf], u0, opts);

    % 3c) Detect lung‐volume peaks (min distance = 1 s)
    vol_ts = U(:,4*N+2);
    [~, locs] = findpeaks(vol_ts, t, 'MinPeakDistance',1000);

    % discard early transient (<2 s)
    locs = locs(locs > 2000);

    % 3d) Classify mean period or Inf if too few peaks
    if numel(locs) < 2
        periods(k) = Inf;
    else
        Tints      = diff(locs)/1000;   % to seconds
        periods(k) = mean(Tints);
    end
end

%% 4) Plot results
figure('Color','w');
plot(gsyn_vals, periods, '-o', 'LineWidth',2);
xlabel('g_{syn} (nS)');
ylabel('Mean cycle period (s)');
ylim([0 10]);
grid on;
title('Network period vs. synaptic coupling (N=20, parallel)');
