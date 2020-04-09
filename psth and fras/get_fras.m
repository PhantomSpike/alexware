function [psth_final,fras] = get_fras(data_root,grid,sorted_root)

start_time_s = 0;
end_time_s = 2100;
chunk_s = 50;
synch_ch = load_synch(data_root,start_time_s,end_time_s,chunk_s);

trig_min_length = 3000; %The minimum length of one trigger in samples

qualityOpt = 'nonoise';
Y = get_spike_times(sorted_root,qualityOpt);


sweep_duration_ms = 200; %The duration of one time interval of interest in ms. Actual Benware sweep duration is around 105.2ms
t_bin_ms = 5; %Time bin in ms
fs_s = 30000; %Sampling rate in samples/s
num_dB_levels = 5; %The number of different dB levels that were presented
num_freq = 15; %The number of different freq
num_stim = num_dB_levels*num_freq; %Totla number of stimuli
% t_ms = [0:t_bin_ms:sweep_duration_ms]; %The edges of the histogram
t_ms = [-50:t_bin_ms:sweep_duration_ms]; %The edges of the histogram

diff_sig = diff(synch_ch); %Find the difference between every n+1 sample - n sample. This tells us the the beginning/end of each sweep

%Find the sample no for the beginning of each sweep
start_ix = find(diff_sig==1); %Note that this diff will differ depending on which sync channel we use so always check it before analyzing
end_ix = find(diff_sig==-1); %Note that this diff will differ depending on which sync channel we use so always check it before analyzing
start_ix = start_ix(1:length(end_ix));
diff_ix = end_ix - start_ix; %Find the length of each sweep in samples
start_ix = start_ix(diff_ix >= trig_min_length); %Keep only the triggers which have length >= minimum triger length
start_ix_ms = (start_ix/fs_s).*1000; %Convert the starting sample numbers to times in ms

num_triggers = length(start_ix);

spike_times_ms = Y(:,1).*1000; % Get the spike times in ms
clusters = Y(:,2);
cluster_id = unique(clusters); %Sort the clusters which are good in ascending order
total_no_clusters = length(cluster_id); %The total number of unqiue clusters

psth = cell(num_stim,total_no_clusters); %Initialize the psth variable

%Find the psths for every cluster, stimulus and repetition and store
%them in a cell array

parfor cluster = 1:total_no_clusters
    current_cluster_id = cluster_id(cluster);
    fprintf('== Processing cluster %.0f/%.0f ==\n',cluster,total_no_clusters);
    for stim = 1:num_stim
        ix_rep = find(grid.randomisedGridSetIdx(1:num_triggers,1)==stim);
        ix_rep_ms = start_ix_ms(ix_rep);
        for rep = 1:length(ix_rep)
            psth{stim,cluster}(rep,:) = histc(spike_times_ms(clusters == current_cluster_id),ix_rep_ms(rep) + t_ms);
        end
        psth{stim,cluster} = psth{stim,cluster}(:,1:end-1); %Delete the last bin which is weird
    end
end
avg_psth = cellfun(@(x)(mean(x,1)),psth,'UniformOutput',false); %Find the average across repetitions
%Convert the cell array into 3D array with dimesnions stimuli x
%clusters x time bins 
int_mat = cellfun(@(x)reshape(x,1,1,[]),avg_psth,'un',0);
psth_final(:,:,:) = cell2mat(int_mat);
fras = reshape(psth_final,num_dB_levels,num_freq,total_no_clusters,length(t_ms)-1); %Reshape the mean_psth into an fra
end

