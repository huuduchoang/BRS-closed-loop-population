%% setup
seeds  = 1:20;                           % your 20 seed files
nSeeds = numel(seeds);
Mvals  = 2e-6 : 0.1e-6 : 18e-6;          % metabolic‐demand sweep
nM     = numel(Mvals);
tf_sec = 60000;                            % seconds per block
warmup = 200000;                         % ms warm‐up
  
%% pre‐alloc
allAvg = nan(nSeeds, nM);
opts    = odeset('RelTol',1e-8,'AbsTol',1e-8);

%% loop over seeds
for si = 1:nSeeds
  seed = seeds(si);
  % load that seed's heterogeneous params
  S = load(sprintf('saved_params_seed%d.mat',seed),'params');
  params = S.params;
  
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
  [~, Uwu] = ode15s(@(t,u) closedloop_populationM(t,u,params), [0 warmup], u0, opts);
  u_prev = Uwu(end, :);
  
  % now sweep M in parallel
  avgPO2 = nan(1,nM);
  parfor ix = 1:nM
    p = params;
    p.M = Mvals(ix);
    u_blk = u_prev;
    % six consecutive 60 s runs
    for blk = 1:6
        t0 = (blk-1)*tf;
        t1 =  blk   *tf;
        [t_blk, U_blk] = ode15s(@(t,u) closedloop_populationM(t,u,p), [t0 t1], u_blk, opts);
        u_blk = U_blk(end, :);  
    end
    % compute time‐avg PaO2 over last block
    po2b = U_blk(:, 4*params.N + 4);
    avgPO2(ix) = trapz(t_blk,po2b) / (t_blk(end)-t_blk(1));
  end
  
  allAvg(si,:) = avgPO2;
end
  
%% compute mean ± SD across seeds
meanPO2 = mean(allAvg,1);
stdPO2  = std (allAvg,[],1);
  
%% plot
x = Mvals*1e5;   % scale M by 1e5 for the axis
y = meanPO2;
e = stdPO2;
  
figure('Color','w');
hold on;
% shaded ±1SD
fill( [x, fliplr(x)], ...
      [y+e, fliplr(y-e)], ...
      [0.8 1 0.8], 'EdgeColor','none');
% mean line
plot(x,y,'k-','LineWidth',2);
  
xlabel('M \times 10^{-5}','FontSize',14);
ylabel('\langle P_{a}O_2\rangle\pm1\,SD (mmHg)','Interpreter','latex','FontSize',14);
grid on; box on;
xlim([min(x) max(x)]);
  
legend('1\,SD envelope','Mean','Location','SouthWest');
