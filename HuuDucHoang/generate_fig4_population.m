% generate_figure4_population_chunks_with_original_axes.m
% Same as before, simulate in 10×60 000 ms chunks and use the original xlim/ylim settings.

clearvars; close all; rng(0);

TruncTol = 1; 

%% 1. Build parameter struct
N = 2;
Eleak_mu    = -65;    Eleak_sigma = 0;
gNaP_mu     = 2.8;    gNaP_sigma  = 0;
gsyn_mu     = 0.1;    gsyn_sigma  = 0;
phi_mu      = 0.3;    phi_sigma   = 0;
thetaO2_mu  = 85;     thetaO2_sigma = 0;
sigmaO2_mu  = 30;     sigmaO2_sigma = 0;

params = setup_params( ...
    N, ...
    Eleak_mu,    Eleak_sigma, ...
    gNaP_mu,     gNaP_sigma, ...
    gsyn_mu,     gsyn_sigma, ...
    phi_mu,      phi_sigma, ...
    thetaO2_mu,  thetaO2_sigma, ...
    sigmaO2_mu,  sigmaO2_sigma);

opts = odeset('RelTol',1e-8,'AbsTol',1e-8);
nChunks = 10;
Tchunk  = 60000;  % ms

%% 2. Simulate open-loop in 10 chunks (panel A)
params.gtonic = 0.30;
inits1 = [-55, 0.001, 0.74, 0.0001, 2, 100, 100];
v0       = repmat(inits1(1),N,1);
n0       = repmat(inits1(2),N,1);
h0       = repmat(inits1(3),N,1);
s0       = zeros(N,1);
alpha0   = inits1(4);
voll0    = inits1(5);
PO2lung0 = inits1(6);
PO2blood0= inits1(7);

u0_open = [v0; n0; h0; s0; alpha0; voll0; PO2lung0; PO2blood0];

topen = [];
uopen = [];
t0 = 0; uprev = u0_open;
for i=1:nChunks
    tspan = [t0, t0+Tchunk];
    [t_seg, u_seg] = ode15s(@(t,u) openloop_population(t,u,params), tspan, uprev, opts);
    if i>1
        t_seg = t_seg(2:end); u_seg = u_seg(2:end,:);
    end
    topen = [topen; t_seg];
    uopen = [uopen; u_seg];
    t0    = topen(end);
    uprev = uopen(end,:);
end

vopen = uopen(:,1);
hopen = uopen(:,2*N+1);

%% 3. Simulate closed-loop in 10 chunks (panels B–D)
inits2 = [-50.0010863, 0.0052147493, 0.5152414947, 0.0007852332, 2.1785621040, 77.39545093, 76.9];
v0       = repmat(inits2(1),N,1);
n0       = repmat(inits2(2),N,1);
h0       = repmat(inits2(3),N,1);
s0       = zeros(N,1);
alpha0   = inits2(4);
voll0    = inits2(5);
PO2lung0 = inits2(6);
PO2blood0= inits2(7);

u0_closed = [v0; n0; h0; s0; alpha0; voll0; PO2lung0; PO2blood0];

tclosed = [];
uclosed = [];
t0 = 0; uprev = u0_closed;
for i=1:nChunks
    tspan = [t0, t0+Tchunk];
    [t_seg, u_seg] = ode15s(@(t,u) closedloop_population(t,u,params), tspan, uprev, opts);
    if i>1
        t_seg = t_seg(2:end); u_seg = u_seg(2:end,:);
    end
    tclosed = [tclosed; t_seg];
    uclosed = [uclosed; u_seg];
    t0       = tclosed(end);
    uprev    = uclosed(end,:);
end

vclosed = uclosed(:,1);
hclosed = uclosed(:,2*N+1);

% Identify one burst window for panel overlays
regburstmin = 177270-4900;
regburstmax = 178570;
regburstInds = find(tclosed>=regburstmin & tclosed<=regburstmax);

%% 5. Compute h-nullcline
vrange = -70:0.1:0;
hnull  = 1./(1+exp((vrange - (-48))/6));

%% 6. Plot panels A–D using split‐and‐style logic
figure('Units','normalized','Position',[0.1 0.1 0.8 0.6]);
set(gcf,'DefaultAxesFontSize',14);

% Panel A: open‐loop gtonic=0.30
subplot(2,4,1); hold on;
D      = dlmread('2_cell_data/gtonic_0.3_gsyn_0.1.dat');
hA     = D(:,1);    vA     = D(:,2);
styleA = D(:,4);    branchA = D(:,6);
plot_segments(hA, vA, styleA, branchA, TruncTol, 'k',3,'k--',1);
plot(hopen, vopen, 'b','LineWidth',1);
xlim([0.53 .93]);  ylim([-70 10]);
xlabel('$h$','Interpreter','Latex');
ylabel('$V$','Interpreter','Latex');
set(gca,'XTick',[.53 .73 .93]);
title('Open loop ($g_\mathrm{tonic}=0.3$)','FontSize',18,'Interpreter','Latex');

subplot(2,4,5); hold on;
plot_segments(hA, vA, styleA, branchA, TruncTol, 'k',3,'k--',1);
plot(hopen, vopen, 'b','LineWidth',1);
plot(hnull, vrange, '--','Color',[.5 .5 .5],'LineWidth',3);
plot(.6090,-50.6576,'o','Color',[.5 .5 .5],'MarkerSize',10,'LineWidth',3);
xlim([0.56 .62]);  ylim([-55 -47]);
xlabel('$h$','Interpreter','Latex');
ylabel('$V$','Interpreter','Latex');
set(gca,'XTick',[.56 .59 .62]);

% Panel B: closed‐loop gtonic=0.12
subplot(2,4,2); hold on;
D      = dlmread('2_cell_data/gtonic_0.12_gsyn_0.1.dat');
hB     = D(:,1);    vB     = D(:,2);
styleB = D(:,4);    branchB = D(:,6);
plot_segments(hB, vB, styleB, branchB, TruncTol, 'c',3,'k',1);
plot(hclosed(regburstInds), vclosed(regburstInds), 'k');
plot(hclosed(regburstInds(9679)), vclosed(regburstInds(9679)), ...
     'co','MarkerFaceColor','c','MarkerSize',10);
xlim([0.6 1]); ylim([-70 10]);
xlabel('$h$','Interpreter','Latex');
ylabel('$V$','Interpreter','Latex');
set(gca,'box','off');
title('Closed loop ($g_\mathrm{tonic}=0.12$)','FontSize',18,'Interpreter','Latex');

subplot(2,4,6); hold on;
plot_segments(hB, vB, styleB, branchB, TruncTol, 'c',3,'k',1);
plot(hclosed(regburstInds), vclosed(regburstInds), 'k');
plot(hclosed(regburstInds(9679)), vclosed(regburstInds(9679)), ...
     'co','MarkerFaceColor','c','MarkerSize',10);
plot(hnull, vrange, '--','Color',[.5 .5 .5],'LineWidth',3);
xlim([0.65 0.9]); ylim([-62 -45]);
xlabel('$h$','Interpreter','Latex');
ylabel('$V$','Interpreter','Latex');
set(gca,'box','off','XTick',[.65 .75 .85]);

% Panel C: closed‐loop gtonic=0.22
subplot(2,4,3); hold on;
D      = dlmread('2_cell_data/gtonic_0.22_gsyn_0.1.dat');
hC     = D(:,1);    vC     = D(:,2);
styleC = D(:,4);    branchC = D(:,6);
plot_segments(hC, vC, styleC, branchC, TruncTol, 'g',3,'k',1);
plot(hclosed(regburstInds), vclosed(regburstInds), 'k');
plot(hclosed(292457), vclosed(292457), 'go','MarkerFaceColor','g','MarkerSize',10);
xlim([0.6 1]); ylim([-70 10]);
xlabel('$h$','Interpreter','Latex');
ylabel('$V$','Interpreter','Latex');
set(gca,'box','off','XTick',[.6 .8 1]);
title('Closed loop ($g_\mathrm{tonic}=0.22$)','FontSize',18,'Interpreter','Latex');

subplot(2,4,7); hold on;
plot_segments(hC, vC, styleC, branchC, TruncTol, 'g',3,'k',1);
plot(hclosed(regburstInds), vclosed(regburstInds), 'k');
plot(hclosed(292457), vclosed(292457), 'go','MarkerFaceColor','g','MarkerSize',10);
plot(hnull, vrange, '--','Color',[.5 .5 .5],'LineWidth',3);
xlim([0.65 0.9]); ylim([-62 -45]);
xlabel('$h$','Interpreter','Latex');
ylabel('$V$','Interpreter','Latex');
set(gca,'box','off','XTick',[.65 .75 .85]);

% Panel D: closed‐loop gtonic=0.18
subplot(2,4,4); hold on;
D      = dlmread('2_cell_data/gtonic_0.18_gsyn_0.1.dat');
hD     = D(:,1);    vD     = D(:,2);
styleD = D(:,4);    branchD = D(:,6);
plot_segments(hD, vD, styleD, branchD, TruncTol, 'm',3,'k',1);
plot(hclosed(regburstInds), vclosed(regburstInds), 'k');
plot(hclosed(301655),   vclosed(301655),   'mo','MarkerFaceColor','m','MarkerSize',10);
xlim([0.6 1]); ylim([-70 10]);
xlabel('$h$','Interpreter','Latex');
ylabel('$V$','Interpreter','Latex');
set(gca,'box','off');

subplot(2,4,8); hold on;
plot_segments(hD, vD, styleD, branchD, TruncTol, 'm',3,'k',1);
plot(hclosed(regburstInds), vclosed(regburstInds), 'k');
plot(hclosed(301655),   vclosed(301655),   'mo','MarkerFaceColor','m','MarkerSize',10);
plot(hnull, vrange, '--','Color',[.5 .5 .5],'LineWidth',3);
xlim([0.65 0.9]); ylim([-62 -45]);
xlabel('$h$','Interpreter','Latex');
ylabel('$V$','Interpreter','Latex');
set(gca,'box','off','XTick',[.65 .75 .85]);


%% ─── Helper: split & plot each segment ───────────────────────────────────────
function plot_segments(x,y,style,branch,tol,cS,lwS,cU,lwU)
    dx = abs(diff(x));   dy = abs(diff(y));
    breaks = [ find(dx>tol | dy>tol) ; numel(x) ];
    i0 = 1;
    for b = breaks'
        % within [i0:b], further split whenever style/branch change
        idx = i0:b;
        ds = diff(style(idx));  db = diff(branch(idx));
        splits = find(ds~=0 | db~=0) + i0;
        edges  = [ splits(:)'  b ];
        i1 = i0;
        for e = edges
            if e - i1 + 1 >= 2
                if style(e)==1 && branch(e)==0
                    plot(x(i1:e), y(i1:e), cS, 'LineWidth', lwS);
                elseif style(e)==2 && branch(e)==0
                    plot(x(i1:e), y(i1:e), cU, 'LineWidth', lwU, 'LineStyle', '--');
                end
            end
            i1 = e + 1;
        end
        i0 = b + 1;
    end
end