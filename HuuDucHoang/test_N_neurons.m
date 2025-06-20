% test_N_neurons.m
clear all, close all

%% Number of neurons and heterogeneity settings
N = 2; 

% Heterogeneity (you can set any sigma to zero if you want no variability)
Eleak_mu    = -65;   Eleak_sigma    = 0.05;
gNaP_mu     =   2.8; gNaP_sigma     = 0;
gsyn_mu     =   0.10;gsyn_sigma     = 0.025;
phi_mu      =   0.3; phi_sigma      = 0;  
thetaO2_mu  =  85;   thetaO2_sigma  = 0;
sigmaO2_mu  =  30;   sigmaO2_sigma  = 0;

% Build the params struct (all draws happen inside)
params = setup_params( ...
    N, ...
    Eleak_mu,    Eleak_sigma, ...
    gNaP_mu,     gNaP_sigma, ...
    gsyn_mu,     gsyn_sigma, ...
    phi_mu,      phi_sigma, ...
    thetaO2_mu,  thetaO2_sigma, ...
    sigmaO2_mu,  sigmaO2_sigma);

%% Initial conditions (use the same inits for each cell)
% Single-cell inits from Diekman et al. 2017 (panel A):
initsA = [-54.9027, 0.0015, 0.7519, 0.0001, 2.0152, 95.9033, 94.7513];

v0       = repmat(initsA(1), N, 1);
n0       = repmat(initsA(2), N, 1);
h0       = repmat(initsA(3), N, 1);
s0       = zeros(N,1);            % synaptic gates
alpha0   = initsA(4);
voll0    = initsA(5);
PO2lung0 = initsA(6);
PO2blood0= initsA(7);

u0 = [v0; n0; h0; s0; alpha0; voll0; PO2lung0; PO2blood0];

%% Integrate with ode15s
tf      = 15000;   % ms
opts    = odeset('RelTol',1e-9,'AbsTol',1e-9);
[t,U]   = ode15s(@(t,u) closedloop_population(t,u,params), [0 tf], u0, opts);

%% Extract V, h, and compute gtonic for all N
V_all      = U(:, 1:N);                % t × N
h_all      = U(:, 2*N+1 : 3*N);        % t × N
PO2blood   = U(:,4*N+4);               % t × 1

% per‐neuron tonic conductance:
phi_vec    = params.phi;
th_vec     = params.thetaO2;
sg_vec     = params.sigmaO2;

gtonic_all = zeros(length(t), N);
for i = 1:N
  gtonic_all(:,i) = phi_vec(i).* ...
    (1 - tanh((PO2blood - th_vec(i))./sg_vec(i)));
end

%% Plot (panel A) for all N
colors = lines(N);   % pick N distinct colors

figure('Position',[100 100 600 800])
set(gcf,'Color','w')

%--- Voltages ---
subplot(3,1,1)
hold on
for i = 1:N
  plot(t/1000, V_all(:,i),'LineWidth',2,'Color',colors(i,:));
end
hold off
xlim([0 tf/1000]), ylim([-70 20])
ylabel('V (mV)')
title(sprintf('Membrane potential, %d uncoupled neurons', N))
legend(arrayfun(@(i) sprintf('Neuron %d',i),1:N,'UniformOutput',false),...
       'Location','best')
box off, grid on

%--- Inactivation h ---
subplot(3,1,2)
hold on
for i = 1:N
  plot(t/1000, h_all(:,i),'LineWidth',2,'Color',colors(i,:));
end
hold off
xlim([0 tf/1000]), ylim([0.5 0.8])
ylabel('h')
box off, grid on

%--- Tonic conductance ---
subplot(3,1,3)
hold on
for i = 1:N
  plot(t/1000, gtonic_all(:,i),'LineWidth',2,'Color',colors(i,:));
end
hold off
xlim([0 tf/1000]), ylim([0 max(gtonic_all(:))*1.1])
ylabel('g_{tonic} (nS)')
xlabel('Time (s)')
box off, grid on
