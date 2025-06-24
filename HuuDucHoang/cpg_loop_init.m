function [sys,x0,par0] = cpg_loop_init
% CPG_LOOP_INIT   Initialization for MATCONT continuation of cpg_loop
%
%   [sys,x0,par0] = cpg_loop_init() returns:
%     sys.ode  = @cpg_loop           % handle to the ODE file
%     sys.npar = 7                   % total number of parameters
%     sys.pnames = { ...             % names (just for display)
%       'N','Eleak','gNaP','gsyn', ...
%       'phi','thetaO2','sigmaO2' }
%     x0    = initial state (4N+4 × 1)
%     par0  = initial parameter vector [N,Eleak,gNaP,gsyn,phi,thetaO2,sigmaO2]
%
%   MATCONT will treat par0(4) = gsyn as the free continuation parameter.

  %% 1) System definition
  sys.ode     = @cpg_loop;
  sys.npar    = 7;
  sys.pnames  = { ...
    'N'       , ...
    'Eleak'   , ...
    'gNaP'    , ...
    'gsyn'    , ... 
    'phi'     , ...
    'thetaO2' , ...
    'sigmaO2' ...
  };

  %% 2) Initial parameter values
  N        = 5;     % number of pacemaker neurons
  Eleak    = -65;    % leak reversal (mV)
  gNaP     =  2.8;   % persistent Na conductance (nS)
  gsyn     =  0.0;   % start with zero coupling
  phi      =  0.3;   % chemosensory gain
  thetaO2  = 85;     % chemosensory half‐activation (mmHg)
  sigmaO2  = 30;     % chemosensory slope (mmHg)

  par0 = [N, Eleak, gNaP, gsyn, phi, thetaO2, sigmaO2];

  %% 3) Initial state vector x0
  % Diekman panel‐A single‐cell inits, replicated for N cells
  initsA    = [-56.8172, 9.5344e-04, 0.7454, 2.0026e-04, ...
                2.0525,  98.9638,  97.7927];
  v0        = repmat(initsA(1), N, 1);
  n0        = repmat(initsA(2), N, 1);
  h0        = repmat(initsA(3), N, 1);
  s0        = zeros(N,1);
  alpha0    = initsA(4);
  voll0     = initsA(5);
  PO2lung0  = initsA(6);
  PO2blood0 = initsA(7);

  x0 = [ v0; n0; h0; s0; alpha0; voll0; PO2lung0; PO2blood0 ];
end
