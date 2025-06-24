%% sweep_gsyn_eL_by_lungvolume.m
% Sweep σ_{E_L} vs σ_{g_syn}, classify breathing modes by lung‐volume peaks:
%   1 = eupneic   (3–5 s, CV ≤ 0.1)
%   2 = tachypneic (< 3 s)
%   3 = bradypneic (> 5 s)
%   4 = irregular (3–5 s & CV > 0.1)

%% 0) Parallel pool (once per session)
if isempty(gcp('nocreate'))
    parpool('local');
end

%% 1) Simulation & network settings
clearvars; close all;
N    = 5;             % pacemaker neurons
tf   = 200000;          % ms (20 s)
opts = odeset('RelTol',1e-6,'AbsTol',1e-9);

% fixed means/sigmas
Eleak_mu      = -65;   gNaP_mu       = 2.8;
gNaP_sigma    =  0;    phi_mu        = 0.3;
phi_sigma     =  0;    thetaO2_mu    = 85;
thetaO2_sigma =  0;    sigmaO2_mu    = 30;
sigmaO2_sigma =  0;

% sweep ranges
n1 = 12; n2 = 12;
sigmaEleak_vals = linspace(0,0.2,n1);    % mV
sigmaGsyn_vals  = linspace(0,0.015,n2);   % nS

mode = nan(n1,n2);  % classification matrix

%% 2) Build initial condition vector u0
initsA = [-56.8172, 9.5344e-04, 0.7454, 2.0026e-04, ...
           2.0525, 98.9638, 97.7927];
v0        = repmat(initsA(1),N,1);
n0        = repmat(initsA(2),N,1);
h0        = repmat(initsA(3),N,1);
s0        = zeros(N,1);
alpha0    = initsA(4);
voll0     = initsA(5);
PO2lung0  = initsA(6);
PO2blood0 = initsA(7);
u0 = [v0; n0; h0; s0; alpha0; voll0; PO2lung0; PO2blood0];

%% 3) Parallel sweep: detect lung‐volume peaks
parfor i = 1:n1
  local_row = nan(1,n2);
  for j = 1:n2
    % draw heterogeneities
    params = setup_params( ...
      N, ...
      Eleak_mu,     sigmaEleak_vals(i), ...
      gNaP_mu,      gNaP_sigma, ...
      0.015,         sigmaGsyn_vals(j), ...
      phi_mu,       phi_sigma, ...
      thetaO2_mu,   thetaO2_sigma, ...
      sigmaO2_mu,   sigmaO2_sigma);
    % clamp synapses ≥0
    params.gsyn(params.gsyn<0)=0;

    % integrate
    [t,U] = ode15s(@(t,u) closedloop_population(t,u,params), [0 tf], u0, opts);

    % extract lung volume and detect peaks
    vol_ts = U(:,4*N+2);   % lung volume
    % require peaks ≥1 s apart
    [~,locs] = findpeaks(vol_ts, t);
    % discard transients <2 s
    locs = locs(locs>2000);

    % classify
    if numel(locs)<2
      cls = 3;  % too few → bradypneic
    else
      Tints = diff(locs)/1000;   % sec
      Tmean = mean(Tints);
      Tcv   = std(Tints)/Tmean;
      if Tmean>5
        cls = 3;  % brady
      elseif Tmean<3
        cls = 2;  % tachy
      elseif Tcv>0.1
        cls = 4;  % irregular
      else
        cls = 1;  % eupneic
      end
    end

    local_row(j)=cls;
  end
  mode(i,:)=local_row;
end

%% 4) Plot 4-mode phase diagram
figure('Color','w');
imagesc(sigmaGsyn_vals, sigmaEleak_vals, mode);
set(gca,'YDir','normal','FontSize',12);
xlabel('\sigma_{g_{syn}} (nS)');
ylabel('\sigma_{E_L} (mV)');
title('Breathing mode (lung‐vol peaks)');

cmap = [ 0   1   0;    % eupneic (green)
         1   0   0;    % tachypneic (red)
         0   0   1;    % bradypneic (blue)
       0.7 0.7   0];   % irregular (olive)
colormap(cmap);

c = colorbar;
c.Ticks      = [1 2 3 4];
c.TickLabels = {'eupneic','tachypneic','bradypneic','irregular'};
