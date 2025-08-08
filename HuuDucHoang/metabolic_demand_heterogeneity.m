function metabolic_demand_heterogeneity

odefun = @closedloop_populationM;             

% ---------- sweep settings ----------------------------------------------
seeds     = 1:20;
Mvals     = 2e-6 : 0.1e-6 : 18e-6;
tf_ms     = 60000;             % 60 s
warm_ms   = 180000;            % warm-up
nSeeds    = numel(seeds);  nM = numel(Mvals);

% ---------- solver options (looser) -------------------------------------
opts = odeset('RelTol',1e-6,'AbsTol',1e-6,'Refine',1);

% ---------- cache warm-up per seed --------------------------------------
warmState = cell(nSeeds,1);
u0_proto  = buildInitialState();        % IC template for all seeds
parfor si = 1:nSeeds
    p          = setupParams(seeds(si));
    p.M        = 8e-6;                  % dummy value during warm-up
    [~,U]      = ode15s(odefun,[0 warm_ms],u0_proto,opts,p);
    warmState{si} = U(end,:);           % store final state
end

% ---------- job list (seed, M) ------------------------------------------
[jSeed,jM] = ndgrid(1:nSeeds, 1:nM);    % jSeed is 20x181 matrix with columns of 1:nSeeds, jM is 20x181 matrix with rows of 1:nM
jobs       = [jSeed(:), jM(:)];         % jobs = 3220x2 matrix
allAvgLin  = nan(size(jobs,1),1);

% ---------- spin up pool if needed --------------------------------------
if isempty(gcp('nocreate')), parpool; end

% ---------- main parallel loop ------------------------------------------
nJobs = size(jobs,1);
parfor idx = 1:nJobs
    si = jobs(idx,1);                        % seed index
    mi = jobs(idx,2);                        % M-value index

    p   = setupParams(seeds(si));            % load params
    p.M = Mvals(mi);

    u_blk   = warmState{si};                 % IC from warm-up
    lastAvg = NaN;                           % PaO2 average over final block

    % -------- six consecutive 60-s windows with CUMULATIVE time --------
    for blk = 1:6
        t0 = (blk-1)*tf_ms;
        t1 =  blk    *tf_ms;
        [t_blk,U_blk] = ode15s(odefun,[t0 t1],u_blk,opts,p);
        u_blk         = U_blk(end,:);

        if blk == 6 
            po2 = U_blk(:,4*p.N + 4);
            lastAvg = trapz(t_blk,po2) / (t1 - t0);
        end
    end

    allAvgLin(idx) = lastAvg;                % store (seed, M) result
end   % parfor

% ---------- reshape & analyse -------------------------------------------
allAvg = reshape(allAvgLin,nSeeds,nM);

% ---------- plot --------------------------------------------------------
x = Mvals * 1e5;   

% grab nSeeds distinct colors
colors = distinguishable_colors(nSeeds);

figure('Color','w');
hold on;

% 1) plot each seed’s curve
for si = 1:nSeeds
    plot(x, allAvg(si,:), 'Color', colors(si,:), 'LineWidth', 1.2);
end

% 2) compute mean and SEM
meanPO2 = mean(allAvg, 1);
semPO2  = std(allAvg, 0, 1) / sqrt(nSeeds);

% 3) shaded ± SEM envelope
fill( [x, fliplr(x)], ...
      [meanPO2 + semPO2, fliplr(meanPO2 - semPO2)], ...
      [0.9 0.9 0.9], 'EdgeColor','none');

% 4) mean curve
plot(x, meanPO2, 'k-', 'LineWidth', 2);

xlabel('M \times 10^{-5}', 'FontSize',14);
ylabel('\langle P_{a}O_2\rangle \pm SEM (mmHg)', ...
       'Interpreter','latex','FontSize',14);
grid on; box on;
xlim([min(x) max(x)]);

% if you want a legend:
legend_entries = [arrayfun(@(s) sprintf('seed %d',s),1:nSeeds,'uni',0), ...
                  {'\pm1SEM envelope','Mean'}];
legend(legend_entries, 'Location','SouthWest','Interpreter','none');
end
%%
function p = setupParams(seed)
S = load(sprintf('saved_params_seed%d.mat',seed),'params');
p = S.params;
end
% ------------------------------------------------------------------------
function u0 = buildInitialState
initsA   = [-60.3782 0.0006 0.667918 0.00136 2.3223 99.1582 98.0934];
N        = 20;     % change if you ever vary neuron count per seed
v0 = repmat(initsA(1),N,1);  n0 = repmat(initsA(2),N,1);
h0 = repmat(initsA(3),N,1);  s0 = zeros(N,1);
u0 = [v0; n0; h0; s0; initsA(4:7).'];
end
