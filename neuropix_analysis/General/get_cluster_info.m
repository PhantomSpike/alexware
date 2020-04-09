function clust_info = get_cluster_info(sorted_dir,batch_size)
%This function is a wrapper function for the function
%templatePositionsAmplitudes by Nick Steinmetz
%>>INPUT>>
%sorted_dir - The directory containing the spike sorted data
%batch_size - How many clusters to use per batch during extraction. If you run out of memory decrease batch size
%<<OUTPUT<<
%cluster_info - struct with fields corresponding to the output of templatePositionsAmplitudes

switch nargin
    case 1
        batch_size = 30; 
end

fprintf('== Extracting cluster information ==\n');tic;
%Load the Spike GLX Meta file
meta_dir = dir(fullfile(sorted_dir, '*ap*.meta')); % meta file from spikeGLX specifically
meta_info = readSpikeGLXmeta(fullfile(meta_dir.folder, meta_dir.name));

fs = meta_info.sRateHz; % The sampling rate in Hz
nCh = meta_info.nSavedChans; % The number of saved channels
clust_info.fs = fs;
%Record the quality of each cluster and extract the good and mua ones
%leaving the noise out
fid = fopen(fullfile(sorted_dir,'cluster_group.tsv'),'r'); %Original
dataArray = textscan(fid,'%f%s%[^\n\r]' , 'Delimiter', '\t', 'HeaderLines' ,1, 'ReturnOnError', false);
units = dataArray{1};
units_quality = dataArray{2};

ix_all = ~(strcmp(units_quality,'noise'));
clust_info.clust_id = units(ix_all);
clust_info.clust_quality = units_quality(ix_all);
no_clusters = numel(clust_info.clust_id);

ix_good = strcmp(units_quality,'good');
clust_info.good_units(:,1) = units(ix_good);
ix = ismember(clust_info.clust_id,clust_info.good_units(:,1));
clust_info.good_units(:,2) = find(ix == 1);

ix_mua = strcmp(units_quality,'mua');
clust_info.mua_units(:,1) = units(ix_mua);
ix = ismember(clust_info.clust_id,clust_info.mua_units(:,1));
clust_info.mua_units(:,2) = find(ix == 1);

fclose(fid);

%Load various .npy files containing information and remove the noise clusters from
%them
spike_clusters = readNPY([sorted_dir,'spike_clusters.npy']); 
ix_all_sp = ismember(spike_clusters,clust_info.clust_id);


spike_clusters = spike_clusters(ix_all_sp);
temps = readNPY([sorted_dir,'templates.npy']);
winv = readNPY([sorted_dir,'whitening_mat_inv.npy']);
ycoords = double(readNPY([sorted_dir,'channel_positions.npy']));
ycoords(:,1) = [];
spikeTimes = readNPY([sorted_dir,'spike_times.npy']);
spikeTimes = spikeTimes(ix_all_sp);
clust_info.spikeTimes = double(spikeTimes)./fs; %Convert spike times to s
spike_clusters = double(spike_clusters); %Convert spike clusters to double so we can combine with the spike times variable
clust_info.spikeTimes = [clust_info.spikeTimes,spike_clusters];
spikeTemplates = readNPY([sorted_dir,'spike_templates.npy']);
spikeTemplates = spikeTemplates(ix_all_sp);
tempScalingAmps = readNPY([sorted_dir,'amplitudes.npy']);
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

% To get the true waveforms of the spikes (not just kilosort's template
% shapes), use the getWaveForms function. This takes 2000 random waveforms
% examples throught the recording from a given cluster

win_ms = 1.4; %The window around each spike time to take the spike waveforms
win_samp = ceil((win_ms/1000)*fs);
clust_info.wf_time = [-win_ms:1000/fs:win_ms];

gwfparams.nWf = 2000;                    % Number of waveforms per unit to pull out
gwfparams.dataDir = sorted_dir;    % KiloSort/Phy output folder
% apD = dir(fullfile(sorted_dir, '*ap*.bin')); % AP band file from spikeGLX specifically
apD = dir(fullfile(sorted_dir, '*.bin')); % AP band file from spikeGLX specifically %Alex 14/03/19
gwfparams.fileName = apD(1).name;         % .dat file containing the raw 
gwfparams.dataType = 'int16';            % Data type of .dat file (this should be BP filtered)
gwfparams.nCh = nCh;                      % Number of channels that were streamed to disk in .dat file
gwfparams.wfWin = [-win_samp win_samp];              % Number of samples before and after spiketime to include in waveform

%We need to run the extraction algorithm in batches as otherwise PC will
%run out of memory
num_batch = ceil(no_clusters/batch_size); %Find number of batches to run
batch_ix = [1:batch_size:no_clusters]; 
clust_info.spikeTimeKeeps = [];
clust_info.waveFormsMean = [];
clust_info.raw_wfs = [];
counter = 0;

for batch = 1:num_batch
    fprintf('== Running batch %0.f/%0.f ==\n',batch,num_batch);
    if batch == num_batch
        current_clusters = clust_info.clust_id(batch_ix(batch):no_clusters);
    else
        current_clusters = clust_info.clust_id(batch_ix(batch):batch*batch_size);
    end
    no_curr_clust = numel(current_clusters);
    gwfparams.spikeClusters = spike_clusters(ismember(spike_clusters,current_clusters)); %Alex 22/08/18
    gwfparams.spikeTimes = spikeTimes(ismember(spike_clusters,current_clusters)); % Vector of cluster spike times (in samples) same length as .spikeClusters %Alex 22/08/18
    wf = getWaveForms(gwfparams); %Run the function to extract the waveforms
    [~,ix_max_ch] = max(max(abs(wf.waveFormsMean),[],3),[],2); %Find the channels with the maximum value of the waveform for each cluster
    
    for i = 1:no_curr_clust
        counter = counter + 1;
        clust_info.raw_wfs(counter,:,:) = squeeze(wf.waveForms(i,:,ix_max_ch(i),:));
    end
    clust_info.waveFormsMean = cat(1,clust_info.waveFormsMean,wf.waveFormsMean);
    clust_info.spikeTimeKeeps = [clust_info.spikeTimeKeeps;wf.spikeTimeKeeps];
    clear wf;
    fprintf('== Done! ==\n');
end

clust_info.mean_wfs = squeeze(nanmean(clust_info.raw_wfs,2));
clust_info.std_wfs = squeeze(nanstd(clust_info.raw_wfs,0,2));
clust_info.plot_mean_wfs = squeeze(nanmean(clust_info.waveFormsMean,1));
save([sorted_dir,'/clust_info.mat'],'clust_info');
fprintf('== Done! Processing took %0.fs ==\n',toc);
