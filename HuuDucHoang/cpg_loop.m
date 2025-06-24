function dydt = cpg_loop(t,y,par)
% CPG_LOOP  Closed‐loop respiratory CPG population for MATCONT
%
%   y   : state vector of length 4*N + 4
%     v(1:N)      – membrane potentials
%     n(1:N)      – K⁺ gating
%     h(1:N)      – NaP inactivation
%     s(1:N)      – synaptic gating
%     alpha       – motor‐pool drive
%     vollung     – lung volume
%     PO2lung     – lung O₂
%     PO2blood    – blood O₂
%
%   par : parameter vector
%     par(1) = N           (number of pacemaker neurons)
%     par(2) = Eleak       (E_L in mV)
%     par(3) = gNaP        (g_NaP in nS)
%     par(4) = gsyn        (mean synaptic conductance in nS)
%     par(5) = phi         (chemo‐gain φ in nS)
%     par(6) = thetaO2     (chemo half‐activ in mmHg)
%     par(7) = sigmaO2     (chemo slope in mmHg)
%
% This uses homogeneous parameters (no randn calls) so it’s
% amenable to continuation in gsyn.

  % unpack parameters
  N         = par(1);
  Eleak     = par(2);
  gNaP      = par(3);
  gsyn      = par(4);
  phi       = par(5);
  thetaO2   = par(6);
  sigmaO2   = par(7);

  % unpack state
  v        = y(1:N);
  n        = y(N+1:2*N);
  h        = y(2*N+1:3*N);
  s        = y(3*N+1:4*N);
  alpha    = y(4*N+1);
  vollung  = y(4*N+2);
  PO2lung  = y(4*N+3);
  PO2blood = y(4*N+4);

  %% CPG biophysics
  C  = 21;
  gna= 28; gk = 11.2; gl = 2.8;
  Ena= 50; Ek = -85; El_const = Eleak*ones(N,1);
  Esyn = 0;

  % Persistent Na
  theta_mp = -40; sigma_mp = -6;
  theta_h  = -48; sigma_h  =  6; taumax_h = 10000;
  mp_inf = 1./(1+exp((v - theta_mp)/sigma_mp));
  h_inf  = 1./(1+exp((v - theta_h )/sigma_h ));
  tau_h  = taumax_h ./ cosh((v - theta_h)/(2*sigma_h));
  Inap   = gNaP .* mp_inf .* h .* (v - Ena);

  % Transient Na
  theta_m = -34; sigma_m = -5;
  m_inf   = 1./(1+exp((v - theta_m)/sigma_m));
  Ina     = gna .* (m_inf.^3) .* (1-n) .* (v - Ena);

  % K current
  theta_n = -29; sigma_n = -4; taumax_n = 10;
  Ik      = gk .* (n.^4) .* (v - Ek);
  n_inf   = 1./(1+exp((v - theta_n)/sigma_n));
  tau_n   = taumax_n ./ cosh((v - theta_n)/(2*sigma_n));

  % Leak
  Il = gl .* (v - El_const);

  % Synaptic coupling (homogeneous all‐to‐all)
  Gsyn = gsyn*(ones(N)-eye(N));
  tau_s = 5; k_r = 1;
  theta_s = -10; sigma_s = -5;
  s_inf   = 1./(1+exp(-(v - theta_s)/sigma_s));
  g_in    = Gsyn.' * s;             % sum_j g_ji s_j
  Isyn    = g_in .* (v - Esyn);

  %% Chemosensory feedback
  % single, homogeneous phi, thetaO2, sigmaO2
  gtonic = phi * (1 - tanh((PO2blood - thetaO2)/sigmaO2));
  Itonic = gtonic .* (v - Esyn);

  %% Motor pool
  r    = 0.001; Tmax = 1; VT = 2; Kp = 5;
  NT   = 1./(1+exp(-(v - VT)/Kp));
  Tpop = Tmax * mean(NT);

  %% Lung mechanics & gas exchange
  E1    = 0.0025; E2 = 0.4; Vol0 = 2;
  dvol  = max(0, -E1*(vollung-Vol0) + E2*alpha);
  PO2ext= (760-47)*0.21;  R = 62.364; Temp = 310;
  taulb = 500;

  M      = 8e-6; Hb = 150; volblood = 5;
  eta    = Hb*1.36; gamma = volblood/22400; betaO2 = 0.03;
  c      = 2.5; K = 26;
  SaO2   = (PO2blood^c)/(PO2blood^c + K^c);
  CaO2   = eta*SaO2 + betaO2*PO2blood;
  partial= (c*PO2blood^(c-1))*(1/(PO2blood^c+K^c) - PO2blood^c/(PO2blood^c+K^c)^2);
  Jlb    = (1/taulb)*(PO2lung - PO2blood)*(vollung/(R*Temp));
  Jbt    = M*CaO2*gamma;

  %% Differential equations
  dydt            = zeros(size(y));
  dydt(1:N)       = (-Inap - Ina - Ik - Il - Itonic - Isyn)./C;
  dydt(N+1:2*N)   = (n_inf - n)./tau_n;
  dydt(2*N+1:3*N) = (h_inf - h)./tau_h;
  dydt(3*N+1:4*N) = ((1 - s).*s_inf - k_r.*s)./tau_s;
  dydt(4*N+1)     = r*Tpop*(1-alpha) - r*alpha;
  dydt(4*N+2)     = -E1*(vollung - Vol0) + E2*alpha;
  dydt(4*N+3)     = (1/vollung)*(PO2ext - PO2lung)*dvol - Jlb*(R*Temp/vollung);
  dydt(4*N+4)     = (Jlb - Jbt)/(gamma*(betaO2 + eta*partial));
end
