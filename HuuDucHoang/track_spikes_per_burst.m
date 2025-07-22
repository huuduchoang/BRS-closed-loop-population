% track_spikes_per_burst.m
% Given:
%   time_s   = t/1000;           % time vector in seconds (t from ode15s)
%   PO2_b    = U(:,4*N+4);       % arterial O2 (mmHg)
%   V_all    = U(:,1:N);         % membrane voltages, size [length(t) × N]
%
% Usage: 
%   [burstTimes, Nspike] = track_spikes_per_burst(time_s, PO2_b, V_all);

function [burstTimes, Nspike] = track_spikes_per_burst(time_s, PO2_b, V_all)

  %% 1) Find all local minima (troughs) and maxima (peaks) in PO2_b
  %   We use findpeaks on -PO2_b to get troughs
  [~, troughLocs] = findpeaks(-PO2_b, time_s, ...
                              'MinPeakProminence',1, ...   % adjust if needed
                              'MinPeakDistance',1);         % at least 1 s apart
  [peakVals, peakLocs] = findpeaks( PO2_b, time_s, ...
                              'MinPeakProminence',1, ...
                              'MinPeakDistance',1);

  % Make sure every peak has a preceding trough
  peakLocs = peakLocs(peakLocs > min(troughLocs));
  
  %% 2) For each peak, find the most‐recent trough before it
  burstTimes = [];      % will store the time of each detected inspiration (peak time)
  Nspike     = [];      % number of neurons spiking in that burst

  thresh = -20;         % mV threshold for “spike”

  for k = 1:length(peakLocs)
    t_peak = peakLocs(k);
    % find the last trough before this peak
    prevTroughs = troughLocs(troughLocs < t_peak);
    if isempty(prevTroughs)
      continue
    end
    t_trough = prevTroughs(end);

    % indices for this burst window
    idx = time_s >= t_trough & time_s <= t_peak;
    if sum(idx)<2
      continue
    end

    % count how many neurons ever exceed the threshold
    Vcycle = V_all(idx, :);           % [window × N]
    didSpike = any(Vcycle > thresh);  % 1×N logical
    nSpk = sum(didSpike);

    % record
    burstTimes(end+1) = t_peak;      %#ok<AGROW>
    Nspike(end+1)     = nSpk;         %#ok<AGROW>
  end

  %% 3) Plot results
  figure('Color','w');
  plot(burstTimes, Nspike, 'o-', 'LineWidth',1.5);
  xlabel('Burst time (s)');
  ylabel('Number of neurons spiking');
  title('Spikes per burst');
  grid on;

end
