function plot_cl_info(sorted_dir)

%Load the clust info containing all the data
load(fullfile(sorted_dir,'clust_info'));
%Load the Spike GLX Meta file
meta_dir = dir(fullfile(sorted_dir, '*ap*.meta')); % meta file from spikeGLX specifically
meta_info = readSpikeGLXmeta(fullfile(meta_dir.folder, meta_dir.name));
%Load ks object for O'Shea's code
option = meta_info.imProbeOpt;
switch option
    case {1,2,3}
        chan_map_path = '/home/alex/Desktop/Code/neuropixel-utils/map_files/neuropixPhase3A_kilosortChanMap.mat'; %aboslute path to the Pixels 3A chan map file from KS2
    case 4
        chan_map_path = '/home/alex/Desktop/Code/neuropixel-utils/map_files/neuropixPhase3A_kilosortChanMap_option4.mat';
end
setenv('NEUROPIXEL_MAP_FILE',chan_map_path); % Set the path to the chanmap for the O'Shea code
ks = Neuropixel.KiloSortDataset(sorted_dir); %Load the results for O'Shea code
ks.load();
metrics = ks.computeMetrics();
%% Plot the individual clsuters with their std for both good and mua units separately 
save_dir = [sorted_dir,'/Plots/'];

if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
good_units_dir = fullfile(save_dir,'Good_units');
mua_units_dir = fullfile(save_dir,'MUA_units');

if ~exist(good_units_dir,'dir')
    mkdir(good_units_dir);
end

if ~exist(mua_units_dir,'dir')
    mkdir(mua_units_dir);
end

plot_clusters(clust_info,clust_info.good_units(:,1),good_units_dir);
plot_clusters(clust_info,clust_info.mua_units(:,1),mua_units_dir);

%% Plot all template waveforms for both good and mua units separately
figure('units','normalized','outerposition',[0 0 1 1]);
imagesc(clust_info.clust_waveforms(clust_info.good_units(:,2),:));
set(gca, 'YDir', 'normal'); xlabel('Time (samples)'); ylabel('Cluster #');
title('Good units waveforms')
colormap(colormap_BlueWhiteRed); caxis([-1 1]*max(abs(caxis()))/2); box off;

file_name = ['Spike_wfs_good_clusters'];
save_name = fullfile(good_units_dir,[file_name,'.png']);
export_fig(save_name);
close;

figure('units','normalized','outerposition',[0 0 1 1]);
imagesc(clust_info.clust_waveforms(clust_info.mua_units(:,2),:));
set(gca, 'YDir', 'normal'); xlabel('Time (samples)'); ylabel('Cluster #');
title('MUA units waveforms')
colormap(colormap_BlueWhiteRed); caxis([-1 1]*max(abs(caxis()))/2); box off;

file_name = ['Spike_wfs_mua_clusters'];
save_name = fullfile(mua_units_dir,[file_name,'.png']);
export_fig(save_name);
close

%% Plot a drift map of the channels
figure('units','normalized','outerposition',[0 0 1 1]);
metrics.plotDriftmap();

file_name = ['Spike_Drift_map_total'];
save_name = fullfile([save_dir,file_name,'.png']);
export_fig(save_name);
close;

%% Plot a drift map of individual clusters
figure('units','normalized','outerposition',[0 0 1 1]);
metrics.plotClusterDriftmap('cluster_ids', clust_info.good_units(:,1));
file_name = ['Spike_Drift_map_good'];
save_name = fullfile(good_units_dir,[file_name,'.png']);
export_fig(save_name);
close;

%The same but without the dots and just the trend lines
figure('units','normalized','outerposition',[0 0 1 1]);
metrics.plotClusterDriftmap('cluster_ids', clust_info.good_units(:,1), 'showSmooth', true, 'showIndividual', false);
file_name = ['Spike_Drift_map_good_line'];
save_name = fullfile(good_units_dir,[file_name,'.png']);
export_fig(save_name);
close;

%% Plot clusters center of mass across the probe
figure('units','normalized','outerposition',[0 0 1 1]);
metrics.plotClusterWaveformAtCenterOfMass();
file_name = ['Clusters_COM_all'];
save_name = [save_dir,file_name,'.png'];
export_fig(save_name);
close;
%% Basic quantification of spiking plot (pdf and cdf along the electrode)
max_depth = max(clust_info.spikeDepths) + 20;
depthBins = 0:30:max_depth;
ampBins = 0:20:min(max(clust_info.spikeAmps),800);
recordingDur = clust_info.spikeTimes(end);

% [pdfs, cdfs] = computeWFampsOverDepth(clust_info.spikeAmps, clust_info.spikeDepthsDrift, ampBins, depthBins, recordingDur);
[pdfs, cdfs] = computeWFampsOverDepth(clust_info.spikeAmps, clust_info.spikeDepthsDrift(1:numel(clust_info.spikeAmps)), ampBins, depthBins, recordingDur); %Alex 23/10/18
plotWFampCDFs(pdfs, cdfs, ampBins, depthBins);

file_name = ['pdf_cdf_spike_map'];
save_name = [save_dir,file_name,'.png'];
export_fig(save_name);
close;