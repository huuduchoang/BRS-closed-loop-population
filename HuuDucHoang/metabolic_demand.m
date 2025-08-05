clear all


%% Parameter settings
N = 20;
Eleak_mu    = -65;   Eleak_sigma    = 0;
gNaP_mu     =   2.8; gNaP_sigma     = 0;
gsyn_mu     =   0.1;   gsyn_sigma     = 0;
phi_mu      =   0.3; phi_sigma      = 0;
thetaO2_mu  =  85;   thetaO2_sigma  = 0;
sigmaO2_mu  =  30;   sigmaO2_sigma  = 0;

params = setup_params( ...
    N, ...
    Eleak_mu,    Eleak_sigma, ...
    gNaP_mu,     gNaP_sigma, ...
    gsyn_mu,     gsyn_sigma, ...
    phi_mu,      phi_sigma, ...
    thetaO2_mu,  thetaO2_sigma, ...
    sigmaO2_mu,  sigmaO2_sigma);

params.gtonic = 0.3;
params.M=8e-6;

%% Initial conditions
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

%% ode15s run

tf = 60000;

% closed-loop
opts    = odeset('RelTol',1e-9,'AbsTol',1e-9);
[t0_closed, U0_closed]   = ode15s(@(t,u) closedloop_populationM(t,u,params), [0 tf], u0, opts);
inits1_closed = U0_closed(end, :);
[t1_closed, U1_closed]   = ode15s(@(t,u) closedloop_populationM(t,u,params), [0 tf], inits1_closed, opts);

% open-loop
[t0_open, U0_open]   = ode15s(@(t,u) openloop_populationM(t,u,params), [0 tf], u0, opts);
inits1_open = U0_open(end, :);
[t1_open, U1_open]   = ode15s(@(t,u) openloop_populationM(t,u,params), [0 tf], inits1_open, opts);
%% Sweep vector & pre‐alloc
Mvals = 2e-6:0.1e-6:18e-6;
nM    = numel(Mvals);
avgClosed = nan(1,nM);
avgOpen   = nan(1,nM);
%% Parallel loop
parfor ix = 1:nM
  % 1) Copy params and set this M
  p = params;
  p.M = Mvals(ix);

  % 2) CLOSED‐LOOP: run a single [0→6·tf] and then take last tf
  [t_all, U_all] = ode15s(@(t,u) closedloop_populationM(t,u,p), ...
                          [0 6*tf], U1_closed(end,:), opts);
  idx = t_all >= 5*tf;
  tb = t_all(idx);
  pb = U_all(idx,4*N+4);
  avgClosed(ix) = trapz(tb,pb)/(tb(end)-tb(1));

  % 3) OPEN‐LOOP: same idea
  [t_all, U_all] = ode15s(@(t,u) openloop_populationM(t,u,p), ...
                          [0 6*tf], U1_open(end,:), opts);
  idx = t_all >= 5*tf;
  tb = t_all(idx);
  pb = U_all(idx,4*N+4);
  avgOpen(ix) = trapz(tb,pb)/(tb(end)-tb(1));
end

%% Plot results
figure;
plot(Mvals, avgClosed, 'k-', 'LineWidth', 2);
hold on;
plot(Mvals, avgOpen,   'b-', 'LineWidth', 2);
xlabel('M'); ylabel('Mean P_{a}O_2');
legend('closed','open','Location','Best');