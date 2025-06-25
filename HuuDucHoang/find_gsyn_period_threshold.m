%% find_gsyn_period_threshold.m
% Find g_syn^* where mean period = 3 s (tachypneic cutoff)
clear; close all;

%— 1) Settings —%
N        = 5;            % number of pacemaker neurons
tf       = 200e3;         % ms to simulate each test (100 s)
tol_g    = 1e-5;          % nS precision on g_syn
targetT  = 3;             % target period in s

function Tmean = computePeriod(gsyn)
    %— 1) Settings —%
    N        = 5;            % number of pacemaker neurons
    tf       = 300e3;         % ms to simulate each test (100 s)
    optsODE  = odeset('RelTol',1e-8,'AbsTol',1e-9);

% fixed heterogeneity = 0
    Eleak_mu = -65; gNaP_mu = 2.8;
    phi_mu   = 0.3; thetaO2_mu = 85; sigmaO2_mu = 30;
    % all sigmas zero
    Eleak_sd = 0; gNaP_sd = 0;
    phi_sd   = 0; thetaO2_sd = 0; sigmaO2_sd = 0;

% initial state vector (single‐cell inits × N)
    initsA   = [-56.8172, 9.5344e-04, 0.7454, 2.0026e-04, 2.0525, 98.9638, 94.7513];
    v0       = repmat(initsA(1), N,1);
    n0       = repmat(initsA(2), N,1);
    h0       = repmat(initsA(3), N,1);
    s0       = zeros(N,1);
    alpha0   = initsA(4);
    voll0    = initsA(5);
    po2l0    = initsA(6);
    po2b0    = initsA(7);
    u0_base  = [v0; n0; h0; s0; alpha0; voll0; po2l0; po2b0];
    % build params
    params = setup_params(N, Eleak_mu,Eleak_sd, gNaP_mu,gNaP_sd, gsyn,0, phi_mu,phi_sd, thetaO2_mu,thetaO2_sd, sigmaO2_mu,sigmaO2_sd);
    params.gsyn(params.gsyn<0)=0;

    % simulate
    [tU,U] = ode15s(@(t,u) closedloop_population(t,u,params), [0 tf], u0_base, optsODE);

    % detect lung‐volume bursts
    vol = U(:,4*N+2);
    [~,locs] = findpeaks(vol, tU, 'MinPeakDistance',2000);
    % drop first 10 s transient
    locs = locs(locs>10000);

    if numel(locs)<2
      Tmean = Inf;
    else
      Tints = diff(locs)/1000;    % s
      Tmean = mean(Tints);
    end
end

%— 3) Bracket gL, gH such that period(gL)>targetT and period(gH)<targetT —%
gL = 0;
TL = computePeriod(gL);
gH = 0.05;
TH = computePeriod(gH);
while TH > targetT
  gH = gH*2;
  TH = computePeriod(gH);
  if gH>1
    error('Cannot bracket period crossing in [0,1] nS');
  end
end

%— 4) Bisection —%
while gH - gL > tol_g
  gM = 0.5*(gL+gH);
  TM = computePeriod(gM);
  if TM > targetT
    gL = gM;  % still too slow
  else
    gH = gM;  % now in tachypneic
  end
end

gstar = 0.5*(gL+gH);
fprintf('g_{syn}^* ≈ %.4f nS (period ≈ %.3f s)\n', gstar, computePeriod(gstar));

%— 5) Optional: plot the period curve —%
gsyns = linspace(0,gH,50);
Tvals = arrayfun(@computePeriod, gsyns);
figure;
plot(gsyns, Tvals,'-o','LineWidth',1.5);
hold on;
plot([gstar gstar],[0 max(Tvals)],'--r');
xlabel('g_{syn} (nS)');
ylabel('Mean period (s)');
title('Mean breathing period vs. synaptic coupling');
legend('period','threshold');
grid on;