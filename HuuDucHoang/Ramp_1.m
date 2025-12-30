%% ----- Plots: superimposed traces for both W ---------------------------
FS  = 18;   % axes tick + label base font
LFS = 16;   % legend font
TFS = 20;   % title font

figure('Color','w','Position',[100 100 1250 800]);

% Panel 1 — mean V (for N=1, it's just V)
ax1 = subplot(3,1,1); hold on; box on; grid on;
hV = gobjects(1,numel(Wvals));
for k = 1:numel(Wvals)
    ts    = T_all{k};
    Vmean = mean(U_all{k}.V,2);
    hV(k) = plot(ts, Vmean, '-', 'LineWidth', 1.6, 'Color', colors(k,:));
end
ylabel('\langle V\rangle (mV)','FontSize',FS);
title('Burst structure during a slow ramp of M','FontSize',TFS);
xlim([0, T_all{end}(end)]);

% --- key fix: legend OUTSIDE so it never covers the trace ---
leg1 = legend(ax1, hV, arrayfun(@(w) sprintf('W=%.2f',w), Wvals, 'uni', 0), ...
              'Location','eastoutside');
set(leg1,'Box','off','FontSize',LFS);

% Panel 2 — motor pool (left) and lung volume (right), both runs
ax2 = subplot(3,1,2); hold on; box on; grid on;
yyaxis left
hA = gobjects(1,numel(Wvals));
for k = 1:numel(Wvals)
    ts = T_all{k};
    hA(k) = plot(ts, U_all{k}.alpha, '-', 'LineWidth', 1.4, 'Color', colors(k,:));
end
ylabel('\alpha','FontSize',FS);

yyaxis right
hVol = gobjects(1,numel(Wvals));
for k = 1:numel(Wvals)
    ts = T_all{k};
    hVol(k) = plot(ts, U_all{k}.vol, '--', 'LineWidth', 1.4, 'Color', colors(k,:));
end
ylabel('Vol_{lung}','FontSize',FS);
xlim([0, T_all{end}(end)]);

% (keep your legend choices; just enlarge fonts)
yyaxis left
leg2 = legend(hA, arrayfun(@(w) sprintf('\\alpha, W=%.2f',w), Wvals, 'uni', 0), ...
              'Location','northwest');
set(leg2,'Box','off','FontSize',LFS);

yyaxis right
leg3 = legend(hVol, arrayfun(@(w) sprintf('vol_{lung}, W=%.2f',w), Wvals, 'uni', 0), ...
              'Location','northeast');
set(leg3,'Box','off','FontSize',LFS);

% Panel 3 — PaO2 (left) with M(t) (right axis, single ramp shown once)
ax3 = subplot(3,1,3); hold on; box on; grid on;
hP = gobjects(1,numel(Wvals));
for k = 1:numel(Wvals)
    ts = T_all{k};
    hP(k) = plot(ts, U_all{k}.PO2b, '-', 'LineWidth', 1.6, 'Color', colors(k,:));
end
ylabel('P_{aO_2} (mmHg)','FontSize',FS);
xlabel('Time (s)','FontSize',FS);
xlim([0, T_all{end}(end)]);

yyaxis right
plot(T_all{1}, Mt_all{1}*1e5, 'k:', 'LineWidth', 1.6);
ylabel('M \times 10^{-5}','FontSize',FS);

yyaxis left
leg4 = legend(hP, arrayfun(@(w) sprintf('P_{aO_2}, W=%.2f',w), Wvals, 'uni', 0), ...
              'Location','best');
set(leg4,'Box','off','FontSize',LFS);

% Make tick labels bigger everywhere
set([ax1 ax2 ax3], 'FontSize', FS);

linkaxes([ax1 ax2 ax3],'x');
