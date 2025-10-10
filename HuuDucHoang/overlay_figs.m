function hFig = overlay_figs(baseFigFile, overlayFigFile, dottedStyle)
% OVERLAY_FIGS  Overlay lines from overlayFigFile onto baseFigFile.
% - Copied lines are dotted and tagged so re-runs replace (not pile up).
% - Legend entries are suffixed once: base -> (uncoupled), overlay -> (coupled).
% - Legend is deduplicated and AutoUpdate is off.

if nargin < 3, dottedStyle = ':'; end

% --- Open base and get its first axes
f1  = openfig(baseFigFile,'invisible'); set(f1,'Visible','on');
ax1 = findobj(f1,'Type','axes','-not','Tag','legend','-not','Tag','Colorbar');
ax1 = ax1(1); hold(ax1,'on');

% --- Remove any prior overlay lines we added on earlier runs
oldOver = findobj(ax1,'Type','line','Tag','overlay_dotted');
if ~isempty(oldOver), delete(oldOver); end

% --- Get base lines and label them "(uncoupled)" once
baseLines = findobj(ax1,'Type','line','-not','Tag','legend','-not','Tag','Colorbar');
for k = 1:numel(baseLines)
    nm = strtrim(char(string(get(baseLines(k),'DisplayName'))));
    if isempty(nm), nm = sprintf('Curve %d', k); end
    low = lower(nm);
    if ~contains(low,'(uncoupled)') && ~contains(low,'(coupled)')
        set(baseLines(k),'DisplayName',[nm ' (uncoupled)'],'HandleVisibility','on');
    end
end

% --- Open overlay, copy its lines, make dotted, tag them
f2  = openfig(overlayFigFile,'invisible');
ax2 = findobj(f2,'Type','axes','-not','Tag','legend','-not','Tag','Colorbar');
ax2 = ax2(1);
ln2 = findobj(ax2,'Type','line');
newLines = copyobj(ln2, ax1);
set(newLines, 'LineStyle', dottedStyle, 'Tag','overlay_dotted');

% --- Label copied lines "(coupled)" once
for k = 1:numel(newLines)
    nm = strtrim(char(string(get(newLines(k),'DisplayName'))));
    if isempty(nm), nm = sprintf('Overlay %d', k); end
    low = lower(nm);
    if ~contains(low,'(coupled)') && ~contains(low,'(uncoupled)')
        set(newLines(k),'DisplayName',[nm ' (coupled)'],'HandleVisibility','on');
    end
end

% --- Fit limits to include both figs
xl1 = xlim(ax1); yl1 = ylim(ax1);
xl2 = xlim(ax2); yl2 = ylim(ax2);
xlim(ax1,[min(xl1(1),xl2(1)) max(xl1(2),xl2(2))]);
ylim(ax1,[min(yl1(1),yl2(1)) max(yl1(2),yl2(2))]);
close(f2);

% --- Deduplicate legend entries and freeze it
allLines = findobj(ax1,'Type','line','-not','Tag','Colorbar','-not','Tag','legend');
allLines = flipud(allLines); % keep visual stacking order
names = get(allLines,'DisplayName'); if ~iscell(names), names = {names}; end
names = cellfun(@(s) strtrim(char(string(s))), names, 'UniformOutput', false);
[~, ia] = unique(lower(names),'stable');
uniqLines = allLines(ia); uniqNames = names(ia);

lgd = legend(ax1, uniqLines, uniqNames, 'AutoUpdate','off'); %#ok<NASGU>
drawnow;

if nargout, hFig = f1; end
end
