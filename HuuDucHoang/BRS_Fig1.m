% reproduce_BRS_Fig1_openloop_parfor.m
% Parallel version: each panel runs independently, discards transients,
% and plots the final 5 s of steady‐state bursting or silence.

clear; close all; clc;

%% 1) Parameter grid
gsyn_vals   = [0, 1, 2, 10, 12, 15];    % rows
gtonic_vals = [0.2, 0.3, 0.4, 0.7, 1.0]; % columns

nRows = numel(gsyn_vals);
nCols = numel(gtonic_vals);
nJobs = nRows * nCols;

%% 2) Simulation settings
N      = 2;
tmax   = 65000;   % total simulation time (ms)
tDiscard = 60000; % ms to discard as transient
tspan  = [0 tmax];
opts   = odeset('RelTol',1e-6,'AbsTol',1e-6);

% fixed heterogeneity means
Eleak_mu = -65; Eleak_sigma = 0;
gNaP_mu  =  2.8; gNaP_sigma  = 0;
phi_mu   =  0.3; phi_sigma   = 0;
th_mu    = 85;   th_sigma    = 0;
sg_mu    = 30;   sg_sigma    = 0;

% initial state (from Diekman et al. Panel A)
inits = [-58.5754 0.0006 0.7252 0.0010 2.2665 103.3461 102.2229];
v0  = repmat(inits(1), N, 1);
n0  = repmat(inits(2), N, 1);
h0  = repmat(inits(3), N, 1);
s0  = zeros(N,1);
u0  = [v0; n0; h0; s0; inits(4); inits(5); inits(6); inits(7)];

%% 3) Preallocate results struct
results(nJobs) = struct('row',[], 'col',[], 't',[], 'V',[]);

%% 4) Launch pool if needed
if isempty(gcp('nocreate'))
    parpool('local');
end

%% 5) Parallel sweep
parfor job = 1:nJobs
  % map linear job -> grid (i = row, j = col)
  [i, j] = ind2sub([nRows, nCols], job);

  % build params
  params = setup_params( ...
    N, Eleak_mu,Eleak_sigma, ...
    gNaP_mu, gNaP_sigma, ...
    gsyn_vals(i), 0, ...
    phi_mu,   phi_sigma, ...
    th_mu,    th_sigma, ...
    sg_mu,    sg_sigma);

  % inject tonic drive into params (openloop_population reads params.gtonic)
  params.gtonic = gtonic_vals(j);

  % integrate
  [T, U] = ode15s(@(t,u) openloop_population(t,u,params), tspan, u0, opts);

  % discard transient and record final 5 s
  sel = T >= tDiscard;
  t_ss = (T(sel) - tDiscard) / 1000;  % in seconds
  V1   = U(sel, 1);                   % plot cell 1 only

  % store
  results(job).row = i;
  results(job).col = j;
  results(job).t   = t_ss;
  results(job).V   = V1;
end

%% 6) Plotting
figure('Color','w','Position',[100 100 1000 800]);
for job = 1:nJobs
  i = results(job).row;
  j = results(job).col;
  % invert so high gsyn at top
  rowPlot = nRows - i + 1;
  colPlot = j;

  idx = (rowPlot-1)*nCols + colPlot;
  subplot(nRows, nCols, idx);
    plot(results(job).t, results(job).V, 'k', 'LineWidth',1);
    xlim([0, 5]);    % show 5 s
    ylim([-65, 5]);
    box off; set(gca,'TickDir','out','FontSize',10,'XTick',[],'YTick',[]);

    % annotate left edge
    if colPlot == 1
      ylabel(sprintf('g_{syn} = %d', gsyn_vals(i)));
    end
    % annotate bottom row
    if rowPlot == 1
      xlabel(sprintf('g_{tonic} = %.1f', gtonic_vals(j)));
    end
end

% 5 s scalebar top-right
axes('Position',[0.92 0.92 0.05 0.05],'Visible','off');
text(0,0.5,'5 s','FontSize',12,'HorizontalAlignment','center');
