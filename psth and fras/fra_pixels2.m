function  fra_psth = fra_pixels2(synch_path,grid_root,sorted_root)
load(synch_path);

qualityOpt = 'Both';

%% Define a bunch of parameters
psth_end_ms = 200; %The end of the psth window in ms
psth_start_ms = -100; %The beginning of the psth window
t_bin_ms = 5; %Time bin in ms
sum_window_start_ms = 5; %The summation window for generating FRAs in ms
sum_window_end_ms = 65;
t_ms = [psth_start_ms:t_bin_ms:psth_end_ms]; %The edges of the histogram
sum_window_start_bin = round((sum_window_start_ms - psth_start_ms)/t_bin_ms); %Specifies which bins to extract from the psth for generating the fra
sum_window_end_bin = round((sum_window_end_ms - psth_start_ms)/t_bin_ms);
fs = 30000;

%% Load the grid info and extract relevant information
dir_info = dir([grid_root '/*Info.mat']); %Get the names of all files with this extension
grid_filename = fullfile(dir_info.folder, dir_info.name); %Form the name of the file to be opened
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
all_f = repmat(all_f,1,num_reps);
all_lev = repmat(all_lev,1,num_reps);
all_f = all_f';
all_lev = all_lev';
all_f = all_f(:);
all_lev = all_lev(:);
g_f = num2str(all_f./1000,'%.1f');
g_lev = num2str(all_lev);
%% Find the triggers
[start_time_ms] = get_triggers_new(synch_ch,min_trig_length_s,min_inter_trig_length_s,fs);

num_triggers = numel(start_time_ms);

%% Get the spike times for all MUA and Good clusters
Y = get_spike_times(sorted_root,qualityOpt);
spike_times_ms = Y(:,1).*1000; % Get the spike times in ms
clusters = Y(:,2);
cluster_id = unique(clusters); %Sort the clusters which are good in ascending order
total_no_clusters = length(cluster_id); %The total number of unqiue clusters

fra_psth.params.t_bin_ms = t_bin_ms;
fra_psth.params.dB_levels = dB_lvls;
fra_psth.params.freqs = ceil(freqs)/1000; 
fra_psth.params.t_ms = t_ms(1:end-1);
fra_psth(total_no_clusters).X_dbft = [];
fra_psth(1).cluster_id = cluster_id;



%% Find the psths for every cluster, stimulus and repetition and store them in a cell array

parfor cluster = 1:total_no_clusters
    current_cluster_id = cluster_id(cluster);
    fprintf('== Processing cluster %.0f/%.0f ==\n',cluster,total_no_clusters);
    for stim = 1:num_stim
        ix_rep = find(grid.randomisedGridSetIdx(1:num_triggers,1)==stim);
        ix_rep_ms = start_time_ms(ix_rep);
        for rep = 1:num_reps
            fra_psth(cluster).X_dbft(stim,rep,:) = histc(spike_times_ms(clusters == current_cluster_id),ix_rep_ms(rep) + t_ms); 
        end
    end
    fra_psth(cluster).X_dbft = fra_psth(cluster).X_dbft(:,:,1:end-1); %Delete the last bin which is weird
    psth_temp = squeeze(sum(fra_psth(cluster).X_dbft(:,:,sum_window_start_bin:sum_window_end_bin),3));
    psth_temp = psth_temp';
    psth_temp = psth_temp(:);
    pval(cluster,:) = anovan(psth_temp,{g_f,g_lev},'display','off')';
    fra_psth(cluster).X_dbft = mean(fra_psth(cluster).X_dbft,2);
    fra_psth(cluster).X_dbft = reshape(fra_psth(cluster).X_dbft,num_dB_lvls,num_freqs,numel(t_ms)-1);
    fra_psth(cluster).X_dbft = sum(fra_psth(cluster).X_dbft(:,:,sum_window_start_bin:sum_window_end_bin),3);
    fra_psth(cluster).X_dbft = flipud(fra_psth(cluster).X_dbft);
end
fra_psth(1).pval = pval;
fra_psth(1).pval(:,3) = [1:length(pval)];