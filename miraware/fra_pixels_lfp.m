function  [fra_psth] = fra_pixels_lfp(synch_dir,grid_dir,sorted_dir,qualityOpt)
%This function is very similar to fra_pixels2 with some modifications that
%keeps the individual trials and is used for FRA analysis simialr to the
%one used in the imaging 

dir_inf = dir([synch_dir '/*synch_ch.mat']); %Get the names of all files with this extension
synch_filename = fullfile(synch_dir, dir_inf.name); %Form the name of the file to be opened
load(synch_filename);

%% Define a bunch of parameters
psth_end_ms = 200; %The end of the psth window in ms
psth_start_ms = -100; %The beginning of the psth window
t_bin_ms = 5; %Time bin in ms
t_ms = [psth_start_ms:t_bin_ms:psth_end_ms]; %The edges of the histogram
fs = 30000;

%% Load the grid info and extract relevant information
dir_info = dir([grid_dir '/*Info.mat']); %Get the names of all files with this extension
grid_filename = fullfile(grid_dir, dir_info.name); %Form the name of the file to be opened
grid_load = load(grid_filename);
grid = grid_load.grid;
num_reps = grid.repeatsPerCondition; %Find the number of reps
min_trig_length_ms = grid.stimGrid(1,2); %The minimum length of one trigger in samples
min_trig_length_s = (min_trig_length_ms/1000); %The minimum length of one trigger in s
min_inter_trig_length_s = grid.postStimSilence(1); %The time between two sweeps
% min_inter_trig_length_s = 0.1; %Alex hack !!! 19/07/19
dB_lvls = unique(grid.stimGrid(:,3));
dB_lvls = sort(dB_lvls,'descend');
num_dB_lvls = numel(dB_lvls); %Find the number of different dB levels used
freqs = unique(grid.stimGrid(:,1)); %Find the specific frequencies used
num_freqs = numel(freqs); %Find the number of different frequencies used
num_stim = grid.nStimConditions; %Total number of stimuli

%Define the two vectors that represent the freqs and levels for all the
%repeats that will be used for the anovan test
% all_f = grid.stimGrid(:,1);
% all_lev = grid.stimGrid(:,3);
% all_f = repmat(all_f,1,num_reps);
% all_lev = repmat(all_lev,1,num_reps);
% all_f = all_f';
% all_lev = all_lev';
% all_f = all_f(:);
% all_lev = all_lev(:);
% g_f = num2str(all_f./1000,'%.1f');
% g_lev = num2str(all_lev);
%% Find the triggers
[start_time_ms] = get_triggers_new(synch_ch,min_trig_length_s,min_inter_trig_length_s,fs);

num_triggers = numel(start_time_ms);

%% Get the spike times for all MUA and Good clusters
load(fullfile(sorted_dir,'clust_info.mat'));
spike_times_ms = clust_info.spikeTimes(:,1).*1000; % Get the spike times in ms
clusters = clust_info.spikeTimes(:,2);
switch qualityOpt
    case 'Good'
        cluster_id = clust_info.good_units(:,1);
    case 'MUA'
        cluster_id = clust_info.mua_units(:,1);
    case 'Both'
        cluster_id = clust_info.clust_id;
end

total_no_clusters = length(cluster_id); %The total number of unqiue clusters

fra_psth.params.t_bin_ms = t_bin_ms;
fra_psth.params.dB_levels = dB_lvls;
fra_psth.params.freqs = ceil(freqs)/1000; 
fra_psth.params.t_ms = t_ms;
fra_psth.params.psth_start_ms = psth_start_ms;
fra_psth.params.psth_end_ms = psth_end_ms;
fra_psth(total_no_clusters).X_rt = cell(num_stim);
fra_psth(1).cluster_id = cluster_id;


%% Find the psths for every cluster, stimulus and repetition and store them in a cell array

for cluster = 1:total_no_clusters
    current_cluster_id = cluster_id(cluster);
    fprintf('== Processing cluster %.0f/%.0f ==\n',cluster,total_no_clusters);
    for stim = 1:num_stim
        ix_rep = find(grid.randomisedGridSetIdx(1:num_triggers,1)==stim);
        ix_rep_ms = start_time_ms(ix_rep);
        num_actual_reps = length(ix_rep_ms); %Find out the number of repeats that were played
        for rep = 1:num_actual_reps
            fra_psth(cluster).X_rt{stim}(rep,:) = histc(spike_times_ms(clusters == current_cluster_id),ix_rep_ms(rep) + t_ms); 
        end
    end
    fra_psth(cluster).X_rt = reshape(fra_psth(cluster).X_rt,num_dB_lvls,num_freqs,[]);
    fra_psth(cluster).X_rt = flipud(fra_psth(cluster).X_rt);
end
fra_psth(end).X_rt(:,:,2:end) = []; %Delete extra stuff from the last cell
clust_info.fra_psth = fra_psth;
save_name = fullfile(sorted_dir,'clust_info');
save(save_name,'clust_info');
