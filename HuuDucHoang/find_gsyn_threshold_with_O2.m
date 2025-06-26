function gstar = find_gsyn_threshold_with_O2()
% FIND_GSYN_THRESHOLD_WITH_O2  
%   Bisection search for g_syn* where mean breathing period crosses 3 s,
%   with safety check on blood‐O2 troughs.

  %% 1) Fixed settings
  N        = 5;            % number of pacemaker neurons
  tf       = 600e3;        % ms to simulate each test (600 s)
  tol_g    = 1e-5;         % nS precision on g_syn
  targetT  = 3;            % s, tachypneic cutoff
  safeO2   = 80;           % mmHg, safety threshold

  optsODE = odeset('RelTol',1e-9,'AbsTol',1e-9);

  % fixed heterogeneity = 0
  Eleak_mu = -65;   Eleak_sd = 0;
  gNaP_mu  =   2.8; gNaP_sd  = 0;
  phi_mu   =   0.3; phi_sd   = 0;
  th_mu    =  85;   th_sd    = 0;
  sg_mu    =  30;   sg_sd    = 0;

  % build base initial state (single‐cell inits × N)
  initsA   = [-58.5754 0.0006 0.7252 0.0010 2.2665 103.3461 102.2229];
  v0       = repmat(initsA(1), N,1);
  n0       = repmat(initsA(2), N,1);
  h0       = repmat(initsA(3), N,1);
  s0       = zeros(N,1);
  alpha0   = initsA(4);
  voll0    = initsA(5);
  po2l0    = initsA(6);
  po2b0    = initsA(7);
  u0_base  = [v0; n0; h0; s0; alpha0; voll0; po2l0; po2b0];

  %% 2) Bracket gL, gH so T(gL)>targetT, T(gH)<targetT
  gL = 0;
  [TL,~] = computeMetrics(gL);
  gH = 0.01;
  [TH,~] = computeMetrics(gH);
  while TH > targetT
    gH = 2*gH;
    [TH,~] = computeMetrics(gH);
    if gH > 1
      error('Cannot bracket period crossing in [0,1] nS');
    end 
  end

  %% 3) Bisection
  while gH - gL > tol_g
    gM = 0.5*(gL + gH);
    [TM,~] = computeMetrics(gM);
    if TM > targetT
      gL = gM;   % still too slow
    else
      gH = gM;   % now fast enough
    end
  end
  gstar = 0.5*(gL + gH);

  % final metrics
  [Tstar, minO2star] = computeMetrics(gstar);
  fprintf('g_{syn}^* ≈ %.5f nS\n', gstar);
  fprintf('  mean period = %.3f s, min PO₂ = %.1f mmHg\n', Tstar, minO2star);
  if minO2star < safeO2
    warning('Hypoxia: blood PO₂ dips below %g mmHg', safeO2);
  end

  %% 4) Optional plot
  gs = linspace(0,gstar*1.5,40);
  Tvals   = nan(size(gs));
  minO2   = nan(size(gs));
  for i=1:numel(gs)
    [Tvals(i), minO2(i)] = computeMetrics(gs(i));
  end
  figure;
  yyaxis left
    plot(gs, Tvals,'-o','LineWidth',1.5);
    ylabel('Mean period (s)');
  yyaxis right
    plot(gs, minO2,'-s','LineWidth',1.5);
    ylabel('Min blood PO_{2} (mmHg)');
  xlabel('g_{syn} (nS)');
  title('Period & blood‐O₂ trough vs. g_{syn}');
  grid on;
  legend('Period','Min PO₂','Location','Best');

  %––––––––––––––––––––––––––––––––––––––––––––––––––––––––
  function [Tmean, minPO2] = computeMetrics(gsyn)
    % build params
    params = setup_params(N, ...
      Eleak_mu,Eleak_sd, gNaP_mu,gNaP_sd, ...
      gsyn,0,            phi_mu,phi_sd, ...
      th_mu,th_sd,       sg_mu,sg_sd);
    params.gsyn(params.gsyn<0)=0;

    % simulate
    [tU,U] = ode15s(@(t,u) closedloop_population(t,u,params), ...
                   [0 tf], u0_base, optsODE);

    PO2b = U(:,4*N+4);
    % troughs = peaks of -PO2
    [~,locs] = findpeaks(-PO2b, tU, 'MinPeakDistance',2000);
    locs     = locs(locs>10000);  % drop first 10 s

    if numel(locs)<2
      Tmean = Inf;
    else
      Tints  = diff(locs)/1000;
      Tmean  = mean(Tints);
    end
    if isempty(locs)
      minPO2 = min(PO2b);
    else
      minPO2 = min(interp1(tU,PO2b,locs));
    end
  end

end
