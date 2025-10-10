clear all

%% 1) Common setup
N = 20;
Eleak_mu    = -65;   Eleak_sigma    = 0;
gNaP_mu     =   2.8; gNaP_sigma     = 0;
gsyn_mu     =   0.1;   gsyn_sigma     = 0;
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

params.gtonic = 0.3;

%% 2) Initial conditions
initsA = [-58.5754 0.0006 0.7252 0.0010 2.2665 103.3461 102.2224];
v0       = repmat(initsA(1), N, 1);
n0       = repmat(initsA(2), N, 1);
h0       = repmat(initsA(3), N, 1);
s0       = zeros(N,1);
alpha0   = initsA(4);
voll0    = initsA(5);
PO2lung0 = initsA(6);
PO2blood0= initsA(7);
u0 = [v0; n0; h0; s0; alpha0; voll0; PO2lung0; PO2blood0];

tf   = 60000;
opts = odeset('RelTol',1e-8,'AbsTol',1e-8);

%% 3) Warm‐up at reference M to get a steady starting state
params.M = 8e-6;

[~, U0c] = ode15s(@(t,u) closedloop_populationM(t,u,params), [0 tf],      u0, opts);
[~, U1c] = ode15s(@(t,u) closedloop_populationM(t,u,params), [0 tf], U0c(end,:), opts);
U1_closed_end = U1c(end,:);

[~, U0o] = ode15s(@(t,u) openloop_populationM(t,u,params),   [0 tf],      u0, opts);
[~, U1o] = ode15s(@(t,u) openloop_populationM(t,u,params),   [0 tf], U0o(end,:), opts);
U1_open_end   = U1o(end,:);

%% 4) Sweep M in parallel
Mvals = 18e-6:0.1e-6:36e-6;
nM    = numel(Mvals);

avgClosed = nan(1,nM);
avgOpen   = nan(1,nM);

parfor ix = 1:nM
  % local copy of params & steady seed
  p = params;
  p.M = Mvals(ix);

  %— CLOSED LOOP 6×tf segments —%
  inits = U1_closed_end;
  [t2, U2] = ode15s(@(t,u) closedloop_populationM(t,u,p), [tf     2*tf], inits, opts);
  inits   = U2(end,:);
  [t3, U3] = ode15s(@(t,u) closedloop_populationM(t,u,p), [2*tf  3*tf], inits, opts);
  inits   = U3(end,:);
  [t4, U4] = ode15s(@(t,u) closedloop_populationM(t,u,p), [3*tf  4*tf], inits, opts);
  inits   = U4(end,:);
  [t5, U5] = ode15s(@(t,u) closedloop_populationM(t,u,p), [4*tf  5*tf], inits, opts);
  inits   = U5(end,:);
  [t6, U6] = ode15s(@(t,u) closedloop_populationM(t,u,p), [5*tf  6*tf], inits, opts);

  pb = U6(:,4*N+4);
  avgClosed(ix) = trapz(t6,pb) / (t6(end)-t6(1));

  %— OPEN LOOP 6×tf segments —%
  inits = U1_open_end;
  [t2o, U2o] = ode15s(@(t,u) openloop_populationM(t,u,p), [tf     2*tf], inits, opts);
  inits      = U2o(end,:);
  [t3o, U3o] = ode15s(@(t,u) openloop_populationM(t,u,p), [2*tf  3*tf], inits, opts);
  inits      = U3o(end,:);
  [t4o, U4o] = ode15s(@(t,u) openloop_populationM(t,u,p), [3*tf  4*tf], inits, opts);
  inits      = U4o(end,:);
  [t5o, U5o] = ode15s(@(t,u) openloop_populationM(t,u,p), [4*tf  5*tf], inits, opts);
  inits      = U5o(end,:);
  [t6o, U6o] = ode15s(@(t,u) openloop_populationM(t,u,p), [5*tf  6*tf], inits, opts);

  pb_o = U6o(:,4*N+4);
  avgOpen(ix) = trapz(t6o,pb_o) / (t6o(end)-t6o(1));
end

%% 5) Plot
set(0,'DefaultAxesFontSize',24)
figure; hold on
plot(Mvals,avgClosed,'k-','LineWidth',3)
plot(Mvals,avgOpen,  'b-','LineWidth',3)
xlim([18e-6 36e-6]), ylim([1 140])
xlabel('$M$','Interpreter','latex','FontSize',24)
ylabel('$P_\mathrm{a}O_2$','Interpreter','latex','FontSize',24)
h = legend('closed loop','open loop','Location','Northeast');
set(h,'Interpreter','latex','FontSize',20,'Box','off')
grid on
