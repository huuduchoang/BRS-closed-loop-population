clear all; clc;

%% 1) Common setup: single neuron with optional self-coupling (autapse)
N = 1;

Eleak_mu    = -65;   Eleak_sigma    = 0;
gNaP_mu     =   2.8; gNaP_sigma     = 0;
gsyn_mu     =   0.0; gsyn_sigma     = 0;   % base; we will override via params.W
phi_mu      =   0.3; phi_sigma      = 0;
thetaO2_mu  =  85;   thetaO2_sigma  = 0;
sigmaO2_mu  =  30;   sigmaO2_sigma  = 0;

params = setup_params( ...
    N, ...
    Eleak_mu, Eleak_sigma, ...
    gNaP_mu,  gNaP_sigma, ...
    gsyn_mu,  gsyn_sigma, ...
    phi_mu,   phi_sigma, ...
    thetaO2_mu, thetaO2_sigma, ...
    sigmaO2_mu, sigmaO2_sigma);

params.gtonic = 0.3;     % dynamic in closed loop; used as ref value elsewhere if needed

% --- self-coupling values to test (autapse strength) ---
gSelfVals = 0:0.02:0.20;    % 0, 0.02, ..., 0.20

%% 2) Initial conditions
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

tf   = 60000;                                   % 60 s blocks
opts = odeset('RelTol',1e-8,'AbsTol',1e-8);

%% 3) Sweep ranges
Mvals = 0.2e-5 : 0.1e-6 : 3.6e-5;   % metabolic demand range
nM    = numel(Mvals);
nG    = numel(gSelfVals);

avgClosed = nan(nG,nM);

%% 4) For each g_self: warm up at reference M, then parfor sweep of M (closed loop only)
warmRefM = 8e-6;

for gi = 1:nG
  gself = gSelfVals(gi);

  % local params with autapse:
  p0 = params;
  % Allow either matrix-based or scalar-based RHS implementations:
  p0.W      = gself;    % 1x1 adjacency (autapse)
  p0.softplus_r = 70000;

  % ---- warm-up closed loop at reference M ----
  p0.M = warmRefM;
  [~, U0c] = ode15s(@(t,u) closedloop_populationM_autapse_ReLU(t,u,p0), [0 tf],      u0, opts);
  [~, U1c] = ode15s(@(t,u) closedloop_populationM_autapse_ReLU(t,u,p0), [0 tf], U0c(end,:), opts);
  U1_closed_end = U1c(end,:);

  % ---- sweep M in parallel for this g_self (closed loop only) ----
  parfor ix = 1:nM
    p = p0;  p.M = Mvals(ix);

    % 5 consecutive 60 s to settle + average last (6th) 60 s
    inits = U1_closed_end;
    [~, U2] = ode15s(@(t,u) closedloop_populationM(t,u,p), [tf     2*tf], inits, opts);
    inits   = U2(end,:);
    [~, U3] = ode15s(@(t,u) closedloop_populationM(t,u,p), [2*tf  3*tf], inits, opts);
    inits   = U3(end,:);
    [~, U4] = ode15s(@(t,u) closedloop_populationM(t,u,p), [3*tf  4*tf], inits, opts);
    inits   = U4(end,:);
    [~, U5] = ode15s(@(t,u) closedloop_populationM(t,u,p), [4*tf  5*tf], inits, opts);
    inits   = U5(end,:);
    [t6, U6] = ode15s(@(t,u) closedloop_populationM(t,u,p), [5*tf  6*tf], inits, opts);

    pb = U6(:,4*N+4);  % PaO2 index (with N=1)
    avgClosed(gi,ix) = trapz(t6,pb) / (t6(end)-t6(1));
  end
end

%% 5) Plot: closed loop curves for each g_self
set(0,'DefaultAxesFontSize',16);
figure('Color','w'); hold on; box on; grid on;

% pick colors for different g_self
try
  cmap = distinguishable_colors(nG,[1 1 1]);  %#ok<*NASGU>
catch
  cmap = lines(max(nG,7));
end

for gi = 1:nG
  ccol = cmap(1+mod(gi-1,size(cmap,1)),:);
  plot(Mvals, avgClosed(gi,:), '-', 'LineWidth',2, 'Color', ccol);
end

xlabel('$M$','Interpreter','latex');
ylabel('$\langle P_{a}O_2\rangle$ (mmHg)','Interpreter','latex');
title('Closed loop: single Butera cell with autapse (self-coupling)');
xlim([min(Mvals) max(Mvals)]);

lg = legend(arrayfun(@(g) sprintf('$g_{\\rm self}=%.2f$',g), gSelfVals, 'uni',0), ...
            'Location','northeastoutside','Interpreter','latex');
set(lg,'Box','off');
