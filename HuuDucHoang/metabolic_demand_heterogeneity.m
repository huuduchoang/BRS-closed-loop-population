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

% ---------- solver options ----------------------------------------------
opts = odeset('RelTol',1e-6,'AbsTol',1e-6,'Refine',1);

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
jobs       = [jSeed(:), jM(:)];           % 3220 x 2
nJobs      = size(jobs,1);

% ---------- prefill from any existing per-job files ---------------------
allAvgLin = nan(nJobs,1);
doneMask  = false(nJobs,1);
for idx = 1:nJobs
    si = jobs(idx,1);  mi = jobs(idx,2);
    fname = fullfile(outdir, sprintf('seed%03d_M%03d.mat', si, mi));
    if exist(fname,'file')
        S = load(fname,'lastAvg');
        if isfield(S,'lastAvg') && ~isempty(S.lastAvg)
            allAvgLin(idx) = S.lastAvg;
            doneMask(idx)  = true;
        end
    end
end

% ---------- start pool & a DataQueue to save on the client --------------
if isempty(gcp('nocreate')), parpool; end
dq = parallel.pool.DataQueue;

% client-side state for progress
nDone = sum(doneMask);
fprintf('Already done: %d/%d (%.1f%%)\n', nDone, nJobs, 100*nDone/nJobs);

% callback runs on CLIENT (not in the workers) → OK to use save here
afterEach(dq, @(R) clientSaveOne(R));

% ---------- main parallel loop (workers compute, client saves) ----------
parfor idx = 1:nJobs
    if doneMask(idx), continue; end  % skip completed

    si = jobs(idx,1);
    mi = jobs(idx,2);

    % load params and run six consecutive 60-s windows with cumulative time
    p   = setupParams(seeds(si));
    p.M = Mvals(mi);

    u_blk   = warmState{si};   % start from cached warm-up
    lastAvg = NaN;

    for blk = 1:6
        t0 = (blk-1)*tf_ms;  t1 = blk*tf_ms;
        [t_blk,U_blk] = ode15s(odefun,[t0 t1],u_blk,opts,p);
        u_blk         = U_blk(end,:);
        if blk == 6
            po2     = U_blk(:,4*p.N + 4);
            lastAvg = trapz(t_blk,po2) / (t1 - t0);
        end
    end

    % send result to client for saving (no save in parfor!)
    send(dq, struct('idx',idx,'si',si,'mi',mi,'lastAvg',lastAvg,'outdir',outdir));
end

% ---------- reshape & analyse (allow missing as NaN) --------------------
allAvg = reshape(allAvgLin, nSeeds, nM);
meanPO2 = mean(allAvg,1,'omitnan');
stdPO2  = std(allAvg,0,1,'omitnan');
nEff    = sum(~isnan(allAvg),1);
semPO2  = stdPO2 ./ max(nEff,1).^.5;

% ---------- plot --------------------------------------------------------
x = Mvals * 1e5;
if exist('distinguishable_colors','file'), colors = distinguishable_colors(nSeeds);
else, colors = lines(nSeeds); end

figure('Color','w'); hold on;
for si = 1:nSeeds
    plot(x, allAvg(si,:), 'Color', colors(si,:), 'LineWidth', 1.2);
end
valid = nEff>0;
fill([x(valid) fliplr(x(valid))], ...
     [meanPO2(valid)+semPO2(valid) fliplr(meanPO2(valid)-semPO2(valid))], ...
     [0.9 0.9 0.9], 'EdgeColor','none');
plot(x, meanPO2, 'k-', 'LineWidth', 2);
xlabel('M \times 10^{-5}','FontSize',14);
ylabel('\langle P_{a}O_2\rangle \pm SEM (mmHg)','Interpreter','latex','FontSize',14);
grid on; box on; xlim([min(x) max(x)]);
legend([arrayfun(@(s) sprintf('seed %d',s),1:nSeeds,'uni',0), ...
        {'\pm1SEM envelope','Mean'}], 'Location','SouthWest','Interpreter','none');

fprintf('Finished sweep (some entries may still be NaN if wall-time hit).\n');

% ========= nested client callback (runs on client, can use save) ========
    function clientSaveOne(R)
        % R has fields: idx, si, mi, lastAvg, outdir
        allAvgLin(R.idx) = R.lastAvg;  % update in-client buffer

        fname   = fullfile(R.outdir, sprintf('seed%03d_M%03d.mat', R.si, R.mi));
        tmpname = [fname '.tmp'];

        % save on client; R2021b supports -struct here (not in parfor)
        Ssave = struct('lastAvg', R.lastAvg, 'si', R.si, 'mi', R.mi);
        try
            save(tmpname, '-struct', 'Ssave', '-v7.3');  % or '-v7.3'
            movefile(tmpname, fname, 'f');
        catch ME
            warning('Client save failed for seed=%d, Mindex=%d: %s', R.si, R.mi, ME.message);
            if exist(tmpname,'file'), delete(tmpname); end
        end

        % simple progress ping
        nDoneLocal = sum(~isnan(allAvgLin));
        if mod(nDoneLocal,50)==0
            fprintf('Progress: %d/%d (%.1f%%)\n', nDoneLocal, nJobs, 100*nDoneLocal/nJobs);
        end
    end
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
