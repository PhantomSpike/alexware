function spikeTimes = find_spikes(data,fs,spikeThreshold)
% spikeFilter = makeSpikeFilter(fs)
%
% Makes a filter appropriate for filtering waveforms for spike detection
% 
% fs: Sampling frequency
% 
% The deadTime output is just a guess at how many samples are invalid
% when this filter is used to filter a signal.

Wp = [300 3000];
n = 2;
[spikeFilter.B,spikeFilter.A] = ellip(n, 0.01, 40, Wp/(fs/2));

%Filtert the data as in filterData from Benware
fprintf('== Applying eliptic filter to the data ==\n');tic;
data = double(data);
data = data';
data = filtfilt(spikeFilter.B, spikeFilter.A, data);
fprintf('== Done! Filtering took %.1f sec ==\n',toc);
%From appendSpikeTimes

% Alex
% for ii = 6:(size(data, 2)-6)
%   data(:,ii) = data(:,ii) - median(data(:,(ii-5):(ii+5)),2);
% end
%Alex

mn = mean(data, 1); %mean across all time points for each channel
sd = std(data, [], 1);

for ii = 1:size(data, 2) %for each channel
  data(:,ii) = (data(:,ii)-mn(ii)) / sd(ii) - spikeThreshold;
end



% find threshold crossings
data = sign(data);
data = diff(data)<0;

% append new spike times to spikeTimes
num_chans = size(data, 2);
spikeTimes = cell(1,num_chans);

for chan = 1:num_chans
  spikeSamples = find(data(:,chan));
  spikeTimesMs = (spikeSamples) / fs * 1000;
  spikeTimes{chan} = spikeTimesMs;
end
end