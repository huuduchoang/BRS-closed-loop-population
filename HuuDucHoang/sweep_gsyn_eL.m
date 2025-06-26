%% sweep_gsyn_eL_by_bloodO2_with_safety.m
% Sweep σ_{E_L} vs σ_{g_syn}, classify breathing modes by blood‐O₂ troughs,
% with safety override (PO₂ < 80 mmHg ⇒ tachypneic).

%% 0) Parallel pool (once per session)
if isempty(gcp('nocreate'))
    parpool('local');
end

%% 1) Simulation & network settings
clearvars; close all;
N    = 5;             % pacemaker neurons
tf   = 50000*4;        % ms (200 s)
opts = odeset('RelTol',1e-9,'AbsTol',1e-9);

% fixed means/sigmas
Eleak_mu      = -65;   gNaP_mu       = 2.8;
gNaP_sigma    =  0;    phi_mu        = 0.3;
phi_sigma     =  0;    thetaO2_mu    = 85;
thetaO2_sigma =  0;    sigmaO2_mu    = 30;
sigmaO2_sigma =  0;

% sweep ranges
n1 = 12; n2 = 12;
sigmaEleak_vals = linspace(0,19.5,n1);    % mV
sigmaGsyn_vals  = linspace(0,0.0219,n2);  % nS

mode = nan(n1,n2);  % classification matrix

%% 2) Build initial condition vector u0
initsA = [-58.5754 0.0006 0.7252 0.0010 2.2665 103.3461 102.2229];
v0        = repmat(initsA(1),N,1);
n0        = repmat(initsA(2),N,1);
h0        = repmat(initsA(3),N,1);
s0        = zeros(N,1);
alpha0    = initsA(4);
voll0     = initsA(5);
PO2lung0  = initsA(6);
PO2blood0 = initsA(7);
u0 = [v0; n0; h0; s0; alpha0; voll0; PO2lung0; PO2blood0];

%% 3) Parallel sweep: detect blood‐O₂ troughs + safety
parfor i = 1:n1
  local_row = nan(1,n2);
  for j = 1:n2
    % draw heterogeneities
    params = setup_params( ...
      N, ...
      Eleak_mu,     sigmaEleak_vals(i), ...
      gNaP_mu,      gNaP_sigma, ...
      0.1,        sigmaGsyn_vals(j), ...
      phi_mu,       phi_sigma, ...
      thetaO2_mu,   thetaO2_sigma, ...
      sigmaO2_mu,   sigmaO2_sigma);
    params.gsyn(params.gsyn<0)=0;

    % integrate
    [t,U] = ode15s(@(t,u) closedloop_population(t,u,params), [0 tf], u0, opts);

    % extract blood O₂ and detect troughs (peaks of –PO2)
    PO2b_ts = U(:,4*N+4);
    [~,locs] = findpeaks(-PO2b_ts, t, 'MinPeakDistance',2000);
    locs = locs(locs>2000);  % discard <2 s transient
    trough_vals = interp1(t, PO2b_ts, locs);

    % safety override: any hypoxic dip ⇒ tachypneic
    if any(trough_vals < 80)
      cls = 2;
    elseif numel(locs) < 2
      cls = 3;  % too few troughs ⇒ bradypneic
    else
      Tints = diff(locs)/1000;  % s
      Tmean = mean(Tints);
      Tcv   = std(Tints)/Tmean;
      if Tmean > 5
        cls = 3;  % brady
      elseif Tmean < 3
        cls = 2;  % tachy
      elseif Tcv > 0.1
        cls = 4;  % irregular
      else
        cls = 1;  % eupneic
      end
    end

    local_row(j) = cls;
  end
  mode(i,:) = local_row;
end

%% 4) Plot 4‐mode phase diagram
figure('Color','w');
imagesc(sigmaGsyn_vals, sigmaEleak_vals, mode);
set(gca,'YDir','normal','FontSize',12);
xlabel('\sigma_{g_{syn}} (nS)');
ylabel('\sigma_{E_L} (mV)');
title('Breathing mode (blood‐O_2 troughs + safety)');

cmap = [ 0   1   0;    % eupneic (green)
         1   0   0;    % tachypneic (red)
         0   0   1;    % bradypneic (blue)
       0.7 0.7   0];   % irregular (olive)
colormap(cmap);

c = colorbar;
c.Ticks      = [1 2 3 4];
c.TickLabels = {'eupneic','tachypneic','bradypneic','irregular'};
