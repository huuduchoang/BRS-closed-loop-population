% Duration of chemosensory feedback interruption
breakDur = 80000; % Continue from here
% seed = 6; % somewhere between 136s and 137s)
seed = 2;

load( sprintf('saved_params_seed%d.mat',seed), 'params' );

% %% Number of neurons and heterogeneity settings
% N =20; 
% 
% % Heterogeneity (you can set any sigma to zero if you want no variability)
% Eleak_mu    = -65;   Eleak_sigma    = 2;
% gNaP_mu     =   2.8; gNaP_sigma     = 0.1;
% gsyn_mu     =   0.1;   gsyn_sigma     = 0.005;
% phi_mu      =   0.3; phi_sigma      = 0.003;
% thetaO2_mu  =  85;   thetaO2_sigma  = 8.5;
% sigmaO2_mu  =  30;   sigmaO2_sigma  = 3;
% 
% % Build the params struct (all draws happen inside)
% params = setup_params( ...
%     N, ...
%     Eleak_mu,    Eleak_sigma, ...
%     gNaP_mu,     gNaP_sigma, ...
%     gsyn_mu,     gsyn_sigma, ...
%     phi_mu,      phi_sigma, ...
%     thetaO2_mu,  thetaO2_sigma, ...
%     sigmaO2_mu,  sigmaO2_sigma);

params.gClamp = 0.01;

%% Initial conditions
initsA = [-60.3782 0.0006 0.667918 0.00136 2.3223 99.1582 98.0934];

v0       = repmat(initsA(1), N, 1);
n0       = repmat(initsA(2), N, 1);
h0       = repmat(initsA(3), N, 1);
s0       = zeros(N,1);            
alpha0   = initsA(4);
voll0    = initsA(5);
PO2lung0 = initsA(6);
PO2blood0= initsA(7);

u0 = [v0; n0; h0; s0; alpha0; voll0; PO2lung0; PO2blood0];

%% Simulation run
tf = 6e4;
opts = odeset('RelTol',1e-10,'AbsTol',1e-10);

params.clamp = false;
[t0,U0]   = ode15s(@(t,u) closedloop_population_clampgtonic(t,u,params), [0 tf], u0, opts);

inits1 = U0(end, :);
[t1,U1]   = ode15s(@(t,u) closedloop_population_clampgtonic(t,u,params), [0 tf], inits1, opts);

inits2 = U1(end, :);
params.clamp = true;

[t2,U2]   = ode15s(@(t,u) closedloop_population_clampgtonic(t,u,params), [tf tf+breakDur], inits2, opts);

inits3 = U2(end, :);
params.clamp = false;

[t3,U3]   = ode15s(@(t,u) closedloop_population_clampgtonic(t,u,params), [tf+breakDur tf+breakDur+tf*3], inits3, opts);

t = [t1; t2; t3];
u = [U1; U2; U3];

V_all    = u(:,1:N);
h_all    = u(:,2*N+1:3*N);
h_mean   = mean(h_all,2);
alpha    = u(:,4*N+1);
vol_lung = u(:,4*N+2);
PO2_lung = u(:,4*N+3);
PO2_b    = u(:,4*N+4);

BCinds=find(t>=0 & t<=tf);
DCinds=find(t>=tf & t<=(tf+breakDur));
ACinds=find(t>=(tf+breakDur)); 

%% Draw figures

set(0,'DefaultAxesFontSize',16)

tsec=t/1000; tfsec=max(t)/1000;

xlo=0;
xhi=tfsec;

figure(1)
plot(tsec(BCinds),PO2_b(BCinds),'k','Linewidth',2)
hold on
plot(tsec(DCinds),PO2_b(DCinds),'b','Linewidth',2)
if PO2_b(end)<80
    plot(tsec(ACinds),PO2_b(ACinds),'r','Linewidth',2)
else
    plot(tsec(ACinds),PO2_b(ACinds),'g','Linewidth',2)
end
xlabel('t (s)')
ylabel('P_aO_2')
xlim([xlo xhi])
set(gca,'box','off','XTick',0:50:350)
grid on

figure(2)
plot3(h_mean(BCinds),vol_lung(BCinds),PO2_b(BCinds),'k','Linewidth',2)
hold on
plot3(h_mean(DCinds),vol_lung(DCinds),PO2_b(DCinds),'b','Linewidth',2)
if PO2_b(end)<80
    plot3(h_mean(ACinds),vol_lung(ACinds),PO2_b(ACinds),'r','Linewidth',2)
else
    plot3(h_mean(ACinds),vol_lung(ACinds),PO2_b(ACinds),'g','Linewidth',2)
end
plot3(h_mean(end),vol_lung(end),PO2_b(end),'ko','MarkerSize',8,'MarkerFaceColor','k')
xlabel('$h$','Interpreter','Latex')
ylabel('vol$_\mathrm{L}$','Interpreter','Latex')
zlabel('$P_\mathrm{a}\mathrm{O}_2$','Interpreter','Latex')
grid on