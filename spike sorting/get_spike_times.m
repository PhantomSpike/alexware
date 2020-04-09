function Y = get_spike_times(dataPath,qualityOpt)
%function [Y, expInfo] = loadSpikeTimes(dataPath,qualityOpt)
%>>>INPUT>>>
%dataPath = The absolute path to the folder containing the manually sorted
%Kilosort files
%qualityOpt = Quality option that allows the selection of specific data to
%be included in the final output Y file
%The quality options are:
%nonoise = INclude all the units which were not classified as noise
%Good = Include olny the units that were labelled as single units during
%the manual spike sorting
%MUA = Include olny the units that were labelled as multi units units during
%the manual spike sorting
%Both = Include both the single and multi units 
%All = Include all units, even the ones labelled as noise
%<<<OUTPUT<<<
% Y - [spike absolute time - spike relative times - unit - stimulus # - repeat # - sweep #]
%expInfo - Info file about the stimuli and the recording from the Benware
%software
% Is the channel still relevant ? Or do we just work with units ?

if ~exist('dataPath','var') || isempty(dataPath)
    dataPath = cd;
end
if ~exist('qualityOpt','var') || isempty(qualityOpt)
    qualityOpt = 'Both';
end
load(fullfile(dataPath,'rez.mat'));
% Saving the spike times in mat files
spike_times = rez.st3(:,1);

fs = rez.ops.fs;

spike_units = readNPY(fullfile(dataPath,'spike_clusters.npy'));

if ~isempty(dir(fullfile(dataPath, 'cluster_groups.csv')))
    fid = fopen(fullfile(dataPath,'cluster_groups.csv'),'r');
else
    fid = fopen(fullfile(dataPath,'cluster_group.tsv'),'r');
end

dataArray = textscan(fid,'%f%s%[^\n\r]' , 'Delimiter', '\t', 'HeaderLines' ,1, 'ReturnOnError', false);
units = dataArray{1};
qualia = dataArray{2};

% Filter the units according the result of the manual clustering
switch qualityOpt
    case 'Good'
        idx = strcmp(qualia,'good');
    case 'MUA'
        idx = strcmp(qualia,'mua');
    case 'nonoise'
        idx = ~(strcmp(qualia,'noise'));
    case 'Both'
        idx = strcmp(qualia,'good') | strcmp(qualia,'mua');
    case 'All'
        idx =true(length(units),1);
    otherwise
        error('Unknown quality type.');
end

unitsList = units(idx);
idxSpikes = ismember(double(spike_units),unitsList);
spike_times = spike_times(idxSpikes);
spike_units = spike_units(idxSpikes);

% Build Y matrix
Y(:,1) = spike_times;
Y(:,2) = spike_units;
Y(:,1) = Y(:,1) ./ fs;

end