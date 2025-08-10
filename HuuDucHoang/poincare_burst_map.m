function poincare_burst_map(seed, discardSec, Vth)
% poincare_burst_map(seed, discardSec, Vth)
%
%   seed        : which saved_params_seed#.mat to load (1…20)
%   discardSec  : how many seconds of the simulation to discard
%   Vth         : voltage threshold (mV) that the mean-V must cross
%
% Example: poincare_burst_map(7,20,-35)

% ========= 1) load heterogeneity for this seed ===========================
S       = load(sprintf('saved_params_seed%d.mat',seed),'params');
params  = S.params;
N = params.N;

% ========= 2) build initial state (same for all cells) ===================
initsA = [-60.3782 0.0006 0.667918 0.00136 2.3223 99.1582 98.0934];
v0 = repmat(initsA(1),N,1);
n0 = repmat(initsA(2),N,1);
h0 = repmat(initsA(3),N,1);
s0 = zeros(N,1);
alpha0   = initsA(4);
voll0    = initsA(5);
PO2lung0 = initsA(6);
PO2blood0= initsA(7);
u0 = [v0; n0; h0; s0; alpha0; voll0; PO2lung0; PO2blood0];

% ========= 3) integrate 200 s of the closed loop =========================
odefun = @closedloop_population;
opts   = odeset('RelTol',1e-9,'AbsTol',1e-9);
sim_ms = 300000;                  % 200 s
[t,U]  = ode15s(@(t,u) odefun(t,u,params), [0 sim_ms], u0, opts);

% ========= 4) discard the first discardSec seconds =======================
keep    = t >= discardSec*1000;      % logical index
t       = t(keep);
U       = U(keep,:);

% ========= 5) population-average voltage & burst detection ===============
Vmat    = U(:,1:N);                  % all neuron voltages
Vmean   = mean(Vmat,2);              % mean across neurons

% find upward crossings (start) and downward crossings (end) of Vth
above  = Vmean > Vth;
crossU = find( diff(above) == +1 ) + 1;   % index immediately after crossing up
crossD = find( diff(above) == -1 ) + 1;   % index immediately after crossing down

% ensure pairs start < stop
if isempty(crossU) || isempty(crossD)
    warning('No bursts detected at the chosen threshold.'); return
end
if crossD(1) < crossU(1)      % drop leading downslope if needed
    crossD(1) = [];
end
minLen = min(numel(crossU),numel(crossD));
crossU = crossU(1:minLen);
crossD = crossD(1:minLen);

% ========= 6) area under each burst (trapezoid) ==========================
nBursts = numel(crossU);
areas   = zeros(1,nBursts);
for k = 1:nBursts
    idx = crossU(k):crossD(k);
    areas(k) = trapz(t(idx), Vmean(idx));
end

% ========= 7) Poincaré map ==============================================
Ak   = areas(1:end-1);
Ak1  = areas(2:end);

figure('Color','w'); hold on
plot(Ak,Ak1,'ko','MarkerFaceColor','k')
plot([min(areas) max(areas)],[min(areas) max(areas)],'k:') % 45° diagonal
xlabel('A_k  (mV·ms)','FontSize',12)
ylabel('A_{k+1}  (mV·ms)','FontSize',12)
title(sprintf('Poincaré map, seed %d, V_{th}=%.1f mV',seed,Vth))
axis equal tight; grid on

figure;
plot(t, Vmean)
end
