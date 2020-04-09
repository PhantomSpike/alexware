function sahani_slow(sorted_dir,t_bin_ms)

min_trig_length_s = 40; %The minmum trigger length in seconds. This will depend on the stimulus. It is best practice to check the synch channel always
min_inter_trig_length_s = 0.1; %The minmum inter trigger length in seconds. This will depend on the stimulus. It is best practice to check the synch channel always

%Give the correct Benware Grid folder
grid_list = dir(fullfile(sorted_dir,'Meta'));
grid_list = grid_list(~ismember({grid_list.name},{'.','..'}));
grid_root = fullfile(grid_list(1).folder,grid_list(1).name); 
%Give the synch_ch file
synch_path =   fullfile(sorted_dir,'synch_ch.mat');
%% Get the clust_info
load(fullfile(sorted_dir,'clust_info'),'clust_info');
%% Get the synch channel
load(synch_path,'synch_ch');
%% Get the triggers
fs = clust_info.fs;
[start_time_ms] = get_triggers(synch_ch,min_trig_length_s,min_inter_trig_length_s,fs);
%% Generate the PSTHs 
psth_end_time_s = 36; %This is the beginning of the noise burst that we want to exclude
Y = clust_info.spikeTimes; %Get all the spike times and cluster ids 
psth_sahani = psth_func(Y,start_time_ms,psth_end_time_s,t_bin_ms,grid_root);
%% Run sahani_quick to get NPSP measures
no_clusters = size(psth_sahani,2);
psth_sahani(no_clusters).NPSP = [];
psth_sahani(no_clusters).SP = [];
psth_sahani(no_clusters).NP = [];
psth_sahani(no_clusters).TP = [];
psth_sahani(no_clusters).SP_std_error = [];
psth_sahani(1).t_bin_ms = t_bin_ms;

fprintf('== Calculating NPSPs ==\n');tic;
parfor cluster = 1:no_clusters
    fprintf('== Processing cluster %.0f/%.0f ==\n',cluster,no_clusters);
    [SP, NP, TP, SP_std_error] = sahani_quick(psth_sahani(cluster).data); 
    psth_sahani(cluster).NPSP = NP/SP;
    psth_sahani(cluster).SP = SP;
    psth_sahani(cluster).NP= NP;
    psth_sahani(cluster).TP = TP;
    psth_sahani(cluster).SP_std_error = SP_std_error;
end
fprintf('== Done! Processing took %.0fs ==\n',toc);
%% Save the results
clust_info.psth_sahani = psth_sahani;
save(fullfile(sorted_dir,'clust_info'),'clust_info');