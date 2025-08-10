function metabolic_demand_heterogeneity

odefun = @closedloop_populationM;

% ---------- sweep settings ----------------------------------------------
seeds     = 1:20;
Mvals     = 2e-6 : 0.1e-6 : 18e-6;
tf_ms     = 60000;              % 60 s
warm_ms   = 180000;             % warm-up
nSeeds    = numel(seeds); 
nM        = numel(Mvals);

% ---------- output/checkpoint folder ------------------------------------
outdir = 'mdh_results';
if ~exist(outdir,'dir'), mkdir(outdir); end

% ---------- solver options (looser) -------------------------------------
opts = odeset('RelTol',1e-6,'AbsTol',1e-6);

% ---------- cache warm-up per seed (resume if available) ----------------
warmFile = fullfile(outdir,'warmState.mat');
if exist(warmFile,'file')
    S = load(warmFile,'warmState');
    warmState = S.warmState;
else
    warmState = cell(nSeeds,1);
    u0_proto  = buildInitialState();  % IC template for all seeds
    parfor si = 1:nSeeds
        p     = setupParams(seeds(si));
        p.M   = 8e-6;                 % dummy value during warm-up
        [~,U] = ode15s(odefun,[0 warm_ms],u0_proto,opts,p);
        warmState{si} = U(end,:);     % store final state
    end
    save(warmFile,'warmState','-v7.3');  % so future runs skip warm-up
end

% ---------- job list (seed, M) ------------------------------------------
[jSeed,jM] = ndgrid(1:nSeeds, 1:nM);
jobs       = [jSeed(:), jM(:)];       % 3220 x 2
nJobs      = size(jobs,1);

% ---------- prefill from any existing per-job files ---------------------
allAvgLin = nan(nJobs,1);
for idx = 1:nJobs
    si = jobs(idx,1);  mi = jobs(idx,2);
    fname = fullfile(outdir, sprintf('seed%03d_M%03d.mat', si, mi));
    if exist(fname,'file')
        S = load(fname,'lastAvg');
        if isfield(S,'lastAvg') && ~isempty(S.lastAvg)
            allAvgLin(idx) = S.lastAvg;
        end
    end
end

% ---------- spin up pool if needed --------------------------------------
if isempty(gcp('nocreate')), parpool; end

% ---------- main parallel loop (skips completed jobs) -------------------
parfor idx = 1:nJobs
    % quick skip if already done
    if ~isnan(allAvgLin(idx)), continue; end

    si = jobs(idx,1);
    mi = jobs(idx,2);

    % checkpoint filename for this (seed, Mindex)
    fname = fullfile(outdir, sprintf('seed%03d_M%03d.mat', si, mi));
    if exist(fname,'file')
        % double-check inside parfor (race-safe)
        S = load(fname,'lastAvg');
        if isfield(S,'lastAvg') && ~isempty(S.lastAvg)
            allAvgLin(idx) = S.lastAvg;
            continue
        end
    end

    % load params and run six consecutive 60-s windows with cumulative time
    p   = setupParams(seeds(si));
    p.M = Mvals(mi);

    u_blk   = warmState{si};   % start from cached warm-up
    lastAvg = NaN;

    for blk = 1:6
        t0 = (blk-1)*tf_ms; 
        t1 =  blk   *tf_ms;
        [t_blk,U_blk] = ode15s(odefun,[t0 t1],u_blk,opts,p);
        u_blk         = U_blk(end,:);
        if blk == 6
            po2     = U_blk(:,4*p.N + 4);
            lastAvg = trapz(t_blk,po2) / (t1 - t0);
        end
    end

    allAvgLin(idx) = lastAvg;

    % atomic-ish save: write temp then move
    tmpname = [fname '.tmp'];

    % pack variables into a *scalar* struct
    Ssave = struct('lastAvg', lastAvg, 'si', si, 'mi', mi);

    try
        % parfor-safe save from struct (optionally keep -v7.3)
        save(tmpname, "-fromstruct", Ssave, "-v7");  % or "-v7.3" if needed
        movefile(tmpname, fname, 'f');               % rename atomically
    catch ME
        warning('Save failed for seed=%d, Mindex=%d: %s', si, mi, ME.message);
        if exist(tmpname,'file'), delete(tmpname); end
    end

end

% ---------- reshape & analyse (allow missing as NaN) --------------------
allAvg = reshape(allAvgLin,nSeeds,nM);

% compute mean/SEM with NaNs omitted, SEM uses per-M effective n
meanPO2 = mean(allAvg,1,'omitnan');
stdPO2  = std(allAvg,0,1,'omitnan');
nEff    = sum(~isnan(allAvg),1);
semPO2  = stdPO2 ./ max(nEff,1).^.5;

% ---------- plot --------------------------------------------------------
x = Mvals * 1e5;

% pick colors (fallback if distinguishable_colors is unavailable)
if exist('distinguishable_colors','file')
    colors = distinguishable_colors(nSeeds);
else
    colors = lines(nSeeds);
end

figure('Color','w'); hold on;

% 1) plot each seed’s curve (NaNs break the line where missing)
for si = 1:nSeeds
    plot(x, allAvg(si,:), 'Color', colors(si,:), 'LineWidth', 1.2);
end

% 2) shaded ± SEM envelope (only where nEff>0)
valid = nEff>0;
fill([x(valid) fliplr(x(valid))], ...
     [meanPO2(valid)+semPO2(valid) fliplr(meanPO2(valid)-semPO2(valid))], ...
     [0.9 0.9 0.9], 'EdgeColor','none');

% 3) mean curve
plot(x, meanPO2, 'k-', 'LineWidth', 2);

xlabel('M \times 10^{-5}', 'FontSize',14);
ylabel('\langle P_{a}O_2\rangle \pm SEM (mmHg)', 'Interpreter','latex','FontSize',14);
grid on; box on; xlim([min(x) max(x)]);
legend_entries = [arrayfun(@(s) sprintf('seed %d',s),1:nSeeds,'uni',0), ...
                  {'\pm1SEM envelope','Mean'}];
legend(legend_entries, 'Location','SouthWest','Interpreter','none');

% summary of progress
pctDone = 100*sum(~isnan(allAvgLin))/numel(allAvgLin);
fprintf('Completed %.1f%% (%d/%d jobs). Results in %s\n', pctDone, sum(~isnan(allAvgLin)), numel(allAvgLin), outdir);
end

%% helpers
function p = setupParams(seed)
S = load(sprintf('saved_params_seed%d.mat',seed),'params');
p = S.params;
end

function u0 = buildInitialState
initsA   = [-60.3782 0.0006 0.667918 0.00136 2.3223 99.1582 98.0934];
N        = 20;     % change if you ever vary neuron count per seed
v0 = repmat(initsA(1),N,1);  n0 = repmat(initsA(2),N,1);
h0 = repmat(initsA(3),N,1);  s0 = zeros(N,1);
u0 = [v0; n0; h0; s0; initsA(4:7).'];
end
