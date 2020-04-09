function  fra_psth = fra_pixels2(synch_dir,grid_dir,sorted_dir,qualityOpt)

%% Define a bunch of parameters
sum_window_start_ms = 0; %The start of the summation window for generating FRAs in ms
sum_window_end_ms = 100; %The end of the summation window for generating FRAs in ms
t_bin_ms = 5; %Time bin in ms
psth_end_ms = 250; %The end of the psth window in ms
psth_start_ms = -250; %The beginning of the psth window
base_window_start_ms = -200; %The start of the summation window for calculating baseline rate
base_window_end_ms = -5; %The end of the summation window for calculating the baseline rate
t_ms = [psth_start_ms:t_bin_ms:psth_end_ms]; %The edges of the histogram
sum_window_start_bin = round((sum_window_start_ms - psth_start_ms)/t_bin_ms)+1; %Specifies which bins to extract from the psth for generating the fra
sum_window_end_bin = round((sum_window_end_ms - psth_start_ms)/t_bin_ms)+1;
base_window_start_bin = round((base_window_start_ms - psth_start_ms)/t_bin_ms)+1; %Specifies which bins to extract from the psth for the baseline period
base_window_end_bin = round((base_window_end_ms - psth_start_ms)/t_bin_ms)+1;
fs = 30000;

%% Make names
dir_inf = dir([synch_dir '/*synch_ch.mat']); %Get the names of all files with this extension
synch_filename = fullfile(synch_dir, dir_inf.name); %Form the name of the file to be opened
load(synch_filename);
save_folder = [fullfile(sorted_dir,'Plots')];

if ~exist(save_folder, 'dir')
    mkdir(save_folder)
end
%% Load the grid info and extract relevant information
dir_info = dir([grid_dir '/*Info.mat']); %Get the names of all files with this extension
grid_filename = fullfile(grid_dir, dir_info.name); %Form the name of the file to be opened
grid_load = load(grid_filename);
grid = grid_load.grid;
num_reps = grid.repeatsPerCondition; %Find the number of reps
min_trig_length_ms = grid.stimGrid(1,2); %The minimum length of one trigger in samples
min_trig_length_s = (min_trig_length_ms/1000); %The minimum length of one trigger in s
min_inter_trig_length_s = grid.postStimSilence(1); %The time between two sweeps
dB_lvls = unique(grid.stimGrid(:,3));
dB_lvls = sort(dB_lvls,'descend');
num_dB_lvls = numel(dB_lvls); %Find the number of different dB levels used
freqs = unique(grid.stimGrid(:,1)); %Find the specific frequencies used
num_freqs = numel(freqs); %Find the number of different frequencies used
num_stim = grid.nStimConditions; %Total number of stimuli

%Define the two vectors that represent the freqs and levels for all the
%repeats that will be used for the anovan test
all_f = grid.stimGrid(:,1);
all_lev = grid.stimGrid(:,3);
%% Find the triggers
[start_time_ms] = get_triggers(synch_ch,min_trig_length_s,min_inter_trig_length_s,fs);

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
fra_psth.params.t_ms = t_ms(1:end-1);
min_reps = floor(num_triggers/num_stim);
for jj = 1:total_no_clusters
    fra_psth(jj).X_dbft = cell(num_stim,1);
    for ii = 1:num_stim
        fra_psth(jj).X_dbft{ii} = zeros(min_reps,length(t_ms));
    end
end
fra_psth(1).cluster_id = cluster_id;
fra_psth(total_no_clusters).base_rate = [];


%% Find the psths for every cluster, stimulus and repetition and store them in a cell array

parfor cluster = 1:total_no_clusters
    current_cluster_id = cluster_id(cluster);
    fprintf('== Processing cluster %.0f/%.0f ==\n',cluster,total_no_clusters);
    g_f = [];
    g_lev = [];
    for stim = 1:num_stim
        ix_rep = find(grid.randomisedGridSetIdx(1:num_triggers,1)==stim);
        ix_rep_ms = start_time_ms(ix_rep);
        num_actual_reps = length(ix_rep_ms); %Find out the number of repeats that were played
        temp_f = all_f(stim)*ones(num_actual_reps,1);
        temp_lev = all_lev(stim)*ones(num_actual_reps,1);
        g_f = [g_f;temp_f];
        g_lev = [g_lev;temp_lev];
        for rep = 1:num_actual_reps
            fra_psth(cluster).X_dbft{stim}(rep,:) = histc(spike_times_ms(clusters == current_cluster_id),ix_rep_ms(rep) + t_ms); 
        end
    end
    psth_temp = cellfun(@(x) mean(x(:,sum_window_start_bin:sum_window_end_bin),2),fra_psth(cluster).X_dbft,'UniformOutput',false);
    base_psth_temp = cellfun(@(x) mean(x(:,base_window_start_bin:base_window_end_bin),2),fra_psth(cluster).X_dbft,'UniformOutput',false);
    psth_temp = cell2mat(psth_temp);
    base_psth_temp = cell2mat(base_psth_temp);
    fra_psth(cluster).base_rate = mean(base_psth_temp);
    pval(cluster,:) = anovan(psth_temp,{g_f,g_lev},'display','off','model','interaction')';
    fra_psth(cluster).X_dbft = cellfun(@(x) mean(x),fra_psth(cluster).X_dbft,'UniformOutput',false);
    fra_psth(cluster).X_dbft = cellfun(@(x) mean(x(sum_window_start_bin:sum_window_end_bin)),fra_psth(cluster).X_dbft,'UniformOutput',false);
    fra_psth(cluster).X_dbft = cell2mat(fra_psth(cluster).X_dbft);
    fra_psth(cluster).X_dbft = flipud(reshape(fra_psth(cluster).X_dbft,num_dB_lvls,num_freqs));
end
fra_psth(1).pval = pval;
fra_psth(1).pval(:,4) = [1:length(pval)];
fra_psth(1).pval(:,5) = cluster_id;
fra_psth(1).col1 = 'pval freq anova';
fra_psth(1).col2 = 'pval level anova';
fra_psth(1).col3 = 'pval freq_level interaction anova';
fra_psth(1).col4 = 'cluster ix';
fra_psth(1).col5 = 'cluster id';

save_name = fullfile(save_folder,'fra_psth');
save(save_name,'fra_psth')
