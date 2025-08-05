% Figure Setup
%% Number of neurons and heterogeneity settings
N = 5; 
% seed  = 2;
% rng(seed, 'twister');

% Heterogeneity settings
Eleak_mu    = -65;   Eleak_sigma    = 2;
gNaP_mu     =   2.8; gNaP_sigma     = 0.1;
gsyn_mu     =   0.1;   gsyn_sigma     = 0.003;
phi_mu      =   0.3; phi_sigma      = 0.003;
thetaO2_mu  =  85;   thetaO2_sigma  = 3;
sigmaO2_mu  =  30;   sigmaO2_sigma  = 1;

% Build and save parameters
params = setup_params(N, Eleak_mu, Eleak_sigma, gNaP_mu, gNaP_sigma, ...
                      gsyn_mu, gsyn_sigma, phi_mu, phi_sigma, ...
                      thetaO2_mu, thetaO2_sigma, sigmaO2_mu, sigmaO2_sigma);

%% Initial conditions
initsA = [-60.3782 0.0006 0.667918 0.000964558873889936 2.19871978350039 90.1180326599434 89.1092102145405];
v0       = repmat(initsA(1), N, 1);
n0       = repmat(initsA(2), N, 1);
h0       = repmat(initsA(3), N, 1);
s0       = zeros(N,1);
alpha0   = initsA(4);
voll0    = initsA(5);
PO2lung0 = initsA(6);
PO2blood0= initsA(7);
u0 = [v0; n0; h0; s0; alpha0; voll0; PO2lung0; PO2blood0];

%% Integrate with ode15s for 60000 ms
tf      = 200000;   % ms
opts    = odeset('RelTol',1e-9,'AbsTol',1e-9);
[t,U]   = ode15s(@(t,u) closedloop_population(t,u,params), [0 tf], u0, opts);
time_s  = t/1000;

%% Determine display window: last 15000 ms
win_start = (tf - 20000)/1000;  % seconds
win_end   = tf/1000;            % seconds

%% Extract signals
V_all    = U(:,1:N);
h_all    = U(:,2*N+1:3*N);
alpha    = U(:,4*N+1);
vol_lung = U(:,4*N+2);
PO2_lung = U(:,4*N+3);
PO2_b    = U(:,4*N+4);

%% Track spikes per burst
[burstTimes, Nspike] = track_spikes_per_burst(time_s, PO2_b, V_all);

%% Compute g_tonic
phi_vec = params.phi;
th_vec  = params.thetaO2;
sg_vec  = params.sigmaO2;
gtonic_all = zeros(length(t),N);
for i=1:N
  gtonic_all(:,i) = phi_vec(i).*(1 - tanh((PO2_b - th_vec(i))./sg_vec(i)));
end

%% Figure 1: 2×3 subplots, big fonts, thick lines, LaTeX labels
colors = distinguishable_colors(N);
set(0,'DefaultAxesFontSize',24)
lw = 3;

figure(1)
clf

% Panel 1: CPG voltage traces
subplot(2,3,1)
hold on
for i=1:N
  plot(time_s, V_all(:,i), 'Color', colors(i,:), 'LineWidth', lw);
end
hold off
xlim([win_start win_end])
set(gca,'box','off','TickDir','out')
ylabel('$V_i$ (mV)','Interpreter','latex')
xlabel('$t$ (s)','Interpreter','latex')

% Panel 2: motor‐pool α
subplot(2,3,2)
plot(time_s, alpha, 'k', 'LineWidth', lw)
xlim([win_start win_end])
set(gca,'box','off','TickDir','out')
ylabel('$\alpha$','Interpreter','latex')
xlabel('$t$ (s)','Interpreter','latex')

% Panel 3: lung volume
subplot(2,3,3)
plot(time_s, vol_lung, 'k', 'LineWidth', lw)
xlim([win_start win_end])
set(gca,'box','off','TickDir','out')
ylabel('vol$_\mathrm{L}$','Interpreter','latex')
xlabel('$t$ (s)','Interpreter','latex')

% Panel 4: chemosensation (g_tonic)
subplot(2,3,4)
hold on
for i=1:N
  plot(time_s, gtonic_all(:,i), 'Color', colors(i,:), 'LineWidth', lw);
end
hold off
xlim([win_start win_end])
set(gca,'box','off','TickDir','out')
ylabel('$g_\mathrm{tonic}$','Interpreter','latex')
xlabel('$t$ (s)','Interpreter','latex')

% Panel 5: blood O2
subplot(2,3,5)
plot(time_s, PO2_b, 'k', 'LineWidth', lw)
xlim([win_start win_end])
set(gca,'box','off','TickDir','out')
ylabel('$P_{a}\mathrm{O}_2$','Interpreter','latex')
xlabel('$t$ (s)','Interpreter','latex')

% Panel 6: lung O2
subplot(2,3,6)
plot(time_s, PO2_lung, 'k', 'LineWidth', lw)
xlim([win_start win_end])
set(gca,'box','off','TickDir','out')
ylabel('$P_{A}\mathrm{O}_2$','Interpreter','latex')
xlabel('$t$ (s)','Interpreter','latex')

% optionally maximize the window
set(gcf,'Units','normalized','Position',[0.05 0.05 0.9 0.9]);
