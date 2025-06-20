function params = setup_params( ...
    N, ...
    Eleak_mu,    Eleak_sigma, ...
    gNaP_mu,     gNaP_sigma, ...
    gsyn_mu,     gsyn_sigma, ...
    phi_mu,      phi_sigma, ...
    thetaO2_mu,  thetaO2_sigma, ...
    sigmaO2_mu,  sigmaO2_sigma)
% SETUP_PARAMS  Build only the heterogeneity fields needed by closedloop_population
%
%   params = SETUP_PARAMS(N, Eleak_mu, Eleak_sigma, ...
%                        gNaP_mu, gNaP_sigma, ...
%                        gsyn_mu, gsyn_sigma, ...
%                        phi_mu, phi_sigma, ...
%                        thetaO2_mu, thetaO2_sigma, ...
%                        sigmaO2_mu, sigmaO2_sigma)
%
%   Output fields in params:
%     .N           - number of neurons
%     .Eleak       - N×1 vector of E_L draws
%     .gNaP        - N×1 vector of g_NaP draws
%     .Gsyn        - N×N matrix of g_syn draws (zero diagonal)
%     .phi         - N×1 vector of chemosensory φ draws
%     .thetaO2     - N×1 vector of θ_O2 draws
%     .sigmaO2     - N×1 vector of σ_O2 draws

  params.N = N;

  % leak reversal heterogeneity
  params.Eleak = Eleak_mu    + Eleak_sigma   * randn(N,1);

  % persistent Na+ conductance heterogeneity
  params.gnap  = gNaP_mu     + gNaP_sigma    * randn(N,1);

  % synaptic coupling heterogeneity (all-to-all, zero self)
  G = gsyn_mu + gsyn_sigma .* randn(N);
  G(1:N+1:end) = 0;
  params.gsyn = G;

  % chemosensory drive heterogeneity
  params.phi      = phi_mu      + phi_sigma      * randn(N,1);
  params.thetaO2  = thetaO2_mu  + thetaO2_sigma  * randn(N,1);
  params.sigmaO2  = sigmaO2_mu  + sigmaO2_sigma  * randn(N,1);
end
