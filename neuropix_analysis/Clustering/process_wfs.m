function [clust_info] = process_wfs(sorted_dir,features)
%function [clust_info] = process_wfs(sorted_dir)
%This functions processes the raw wfs of the 'good' clusters so as to
%extract key features for clustering later on. Before feature extraction
%interpolation is done to improve resoltuion and reduce noise in estimates
%>>INPUT>>
%sorted_dir - The directory with the sorted data which contains the
%clust_info file
%<<OUTPUT<<
%clust_info - The updated clust_info variable

%% Params
upsamp_factor = 100; %How many times to upsample the wfs for feature extraction
method = 'makima'; %Which method to use for itnerpolation
plot_on = 1;
if ~exist('features','var')
    features = [1 2];
end
%1 - Spike width (peak-to-trough width), this is the distance form peak (depolarization) to
%trough (hyperpolarization) in the wf, var: spike_width_ms
%2 - Width at 1/2 height trough, this is the width of the trough part at 1/2
%max value, var: width_hh_trough_ms
%3 - End slope, this is the slope of the decaying part of the trough, var: end_slope
%4 - Peak-to-trough ratio, this is the ratio of the amplitudes of the peak
%to the trough PA/TA, var: peak_trough_ratio
%5 - Width at 1/2 height peak, this is the width of the peak part at 1/2
%max value, var: width_hh_peak_ms
%% Load and get name
file_path = fullfile(sorted_dir,'clust_info.mat');
slash_ix1 = find(file_path == '/', 1, 'last');
var_name = file_path(slash_ix1+1:end);
var_name = strrep(var_name,'_',' ');
load(file_path);
%% Interpolation
ix_good_wfs = clust_info.good_units(:,2);
clust_info.good_wfs = clust_info.mean_wfs(ix_good_wfs,:); %Use only waveforms that come from single units
wf_time = clust_info.wf_time;
[clust_info.good_int_wf,clust_info.int_time_ms] = interpolate_wf(clust_info.good_wfs,wf_time,upsamp_factor,method);
%% Extract features
clust_info.good_int_wf = clust_info.good_int_wf.*clust_info.scaleToUv; %Convert the y-scale to uV
[clust_info] = feature_extract(clust_info);
%% Run k-means clustering and save the results
[clust_info] = run_kmeans(clust_info,features,plot_on);
save_dir = fullfile(sorted_dir,'/Plots/Clustering/');
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end
fig_handles = findobj('Type', 'figure');
save_name = fullfile(save_dir,[var_name,'_IE_Clustering_norm.png']);
export_fig(save_name,fig_handles(1));
save_name = fullfile(save_dir,[var_name,'_IE_Clustering.png']);
export_fig(save_name,fig_handles(2));
close all;
save_name = fullfile(sorted_dir,'clust_info.mat');
save(save_name,'clust_info');
