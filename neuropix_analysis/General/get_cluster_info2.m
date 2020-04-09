function clust_info = get_cluster_info2(sorted_dir)
%This function is a wrapper function for the function
%templatePositionsAmplitudes by Nick Steinmetz
%>>INPUT>>
%sorted_dir - The directory containing the spike sorted data
%<<OUTPUT<<
%cluster_info - struct with fields corresponding to the output of templatePositionsAmplitudes




fprintf('== Extracting cluster information ==\n');
%Load the Spike GLX Meta file
meta_dir = dir(fullfile(sorted_dir, '*ap*.meta')); % meta file from spikeGLX specifically
meta_info = readSpikeGLXmeta(fullfile(meta_dir.folder, meta_dir.name));

fs = meta_info.sRateHz; % The sampling rate in Hz
nCh = meta_info.nSavedChans; % The number of saved channels
clust_info.fs = fs;
option = meta_info.imProbeOpt;
switch option
    case {1,2,3}
        chan_map_path = '/home/alex/Desktop/Code/neuropixel-utils/map_files/neuropixPhase3A_kilosortChanMap.mat'; %aboslute path to the Pixels 3A chan map file from KS2
    case 4
        chan_map_path = '/home/alex/Desktop/Code/neuropixel-utils/map_files/neuropixPhase3A_kilosortChanMap_option4.mat';
end
%Record the quality of each cluster and extract the good and mua ones
%leaving the noise out
fid = fopen(fullfile(sorted_dir,'cluster_group.tsv'),'r'); %Original
dataArray = textscan(fid,'%f%s%[^\n\r]' , 'Delimiter', '\t', 'HeaderLines' ,1, 'ReturnOnError', false);
units = dataArray{1};
units_quality = dataArray{2};

ix_all = ~(strcmp(units_quality,'noise'));
all_units_id = units(ix_all);
clust_info.clust_id = all_units_id;
clust_info.clust_quality = units_quality(ix_all);
no_clusters = numel(clust_info.clust_id);

ix_good = strcmp(units_quality,'good');
good_units_id = units(ix_good);
clust_info.good_units(:,1) = good_units_id;
ix = ismember(clust_info.clust_id,clust_info.good_units(:,1));
clust_info.good_units(:,2) = find(ix == 1);

ix_mua = strcmp(units_quality,'mua');
mua_units_id = units(ix_mua);
clust_info.mua_units(:,1) = mua_units_id;
ix = ismember(clust_info.clust_id,clust_info.mua_units(:,1));
clust_info.mua_units(:,2) = find(ix == 1);

fclose(fid);

%Load various .npy files containing information and remove the noise clusters from
%them
spike_clusters = readNPY(fullfile(sorted_dir,'spike_clusters.npy')); 
ix_all_sp = ismember(spike_clusters,clust_info.clust_id);


spike_clusters = spike_clusters(ix_all_sp);
temps = readNPY(fullfile(sorted_dir,'templates.npy'));
winv = readNPY(fullfile(sorted_dir,'whitening_mat_inv.npy'));
ycoords = double(readNPY(fullfile(sorted_dir,'channel_positions.npy')));
ycoords(:,1) = [];
spikeTimes = readNPY(fullfile(sorted_dir,'spike_times.npy'));
spikeTimes = spikeTimes(ix_all_sp);
clust_info.spikeTimes = double(spikeTimes)./fs; %Convert spike times to s
spike_clusters = double(spike_clusters); %Convert spike clusters to double so we can combine with the spike times variable
clust_info.spikeTimes = [clust_info.spikeTimes,spike_clusters];
spikeTemplates = readNPY(fullfile(sorted_dir,'spike_templates.npy'));
spikeTemplates = spikeTemplates(ix_all_sp);
tempScalingAmps = readNPY(fullfile(sorted_dir,'amplitudes.npy'));
tempScalingAmps = tempScalingAmps(ix_all_sp);


%Run Nick's function for extracting the relevant information about the
%clusters
[clust_info.spikeAmps, clust_info.spikeDepths, clust_info.templateDepths, clust_info.tempAmps, clust_info.tempsUnW, clust_info.templateDuration, clust_info.waveforms] = templatePositionsAmplitudes(temps, winv, ycoords, spikeTemplates, tempScalingAmps);

%Get spike Depths which are better suited for the Drift map
[~, ~, clust_info.spikeDepthsDrift, ~] = ksDriftmap(sorted_dir);
%Get the channel site for each spike
% one could potentially use the unwhitened templates, but that shouldn't really change the results
[~,max_site] = max(max(abs(clust_info.tempsUnW),[],2),[],3); % the maximal site for each template
clust_info.templateChannels = max_site;

clust_info.templateDuration = (clust_info.templateDuration./fs).*1000 ; %Convert to ms
clust_info.spikeDepths = clust_info.spikeDepths; %Round the depth
clust_info.templateDepths = round(clust_info.templateDepths); %Round the depth


%Find which templates contributed to which clusters, the depth of
%the clusters, channels they came from, spike waveforms, spike amplitudes
%and duration
for ii = 1:no_clusters
    cluster =  clust_info.clust_id(ii);
    ix = find(spike_clusters == cluster);
    clust_info.no_spikes(ii,1) = numel(ix);
    clust_info.template_id{ii,1} = unique(spikeTemplates(ix)) + 1;
    clust_info.clust_depth(ii,1) = mean(clust_info.templateDepths(clust_info.template_id{ii,1}));
    clust_info.clust_spikedur(ii,1) = mean(clust_info.templateDuration(clust_info.template_id{ii,1}));
    clust_info.clust_amps(ii,1) = mean(clust_info.tempAmps(clust_info.template_id{ii,1}));
    clust_info.clust_waveforms(ii,:) = mean(clust_info.waveforms(clust_info.template_id{ii,1},:),1);
    clust_info.clust_ch(ii,:) = floor(mean(max_site(clust_info.template_id{ii,1},:))) - 1;
end



%% Loading raw waveforms
fprintf('== Extracting raw waveforms ==\n');tic;
num_wfs = 500; %How many raw wfs to extract for each cluster
setenv('NEUROPIXEL_MAP_FILE',chan_map_path); % Set the path to the chanmap for the O'Shea code
ks = Neuropixel.KiloSortDataset(sorted_dir); %Load the results for O'Shea code
ks.load();

num_units = length(all_units_id);
raw_wfs = cell(num_units,1);
raw_wfs_ix = cell(num_units,1);
mean_wfs = [];
std_wfs = [];

parfor cl = 1:num_units
    fprintf('== Extracting cluster %0.f/%0.f ==\n',cl,num_units);
    cluster_id = all_units_id(cl);
    snippet_wfs = ks.getWaveformsFromRawData('cluster_ids', cluster_id, 'num_waveforms', num_wfs, 'best_n_channels', 1, 'car', true,'subtractOtherClusters', false); %Get 500 raw wfs for each good cluster. Perform CRA on the data and remove wfs of other neurons spiking at the same time
    raw_wfs{cl} = double(squeeze(snippet_wfs.data));
    mean_wfs(cl,:) = nanmean(raw_wfs{cl},2);
    std_wfs(cl,:) = nanstd(raw_wfs{cl},0,2);
    raw_wfs_ix{cl} = snippet_wfs.sample_idx;
end

clust_info.raw_wfs = raw_wfs;
clust_info.mean_wfs = mean_wfs;
clust_info.std_wfs = std_wfs;
clust_info.raw_wfs_ix = raw_wfs_ix;

snip = ks.getWaveformsFromRawData('cluster_ids', all_units_id(1), 'num_waveforms', 10, 'best_n_channels', 1, 'car', true,'subtractOtherClusters', true); %Run again to get some basic parameters about the wfs
clust_info.scaleToUv = snip.scaleToUv;
clust_info.wf_time = snip.time_ms;

save([sorted_dir,'/clust_info.mat'],'clust_info');
fprintf('== Done! Processing took %0.fs ==\n',toc);