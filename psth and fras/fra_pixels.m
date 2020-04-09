load('/home/phant0msp1ke/Desktop/Code/alexware/synch_channels/P02_quning_synch_ch');
grid_root = '/media/phant0msp1ke/4TB SSD/DATA/Ferret_exp_17_10_2017/Success_ferret_grids/P02-quning/';
sorted_root = '/media/phant0msp1ke/4TB SSD/DATA/Ferret_exp_17_10_2017/Penetration2/P02-quning/CRA/P02-quning_cra_allch_3000Hz/';
qualityOpt = 'nonoise';

% synch_ch = get_synch(data_root); %Load the synch channel

sweep_duration_ms = 200; %The duration of one time interval of interest in ms. Actual Benware sweep duration is around 105.2ms
t_bin_ms = 5; %Time bin in ms
sum_window_ms = 100; %The summation window for generating FRAs in ms
fs = 30000;

% if ~exist('qualityOpt','var')
%     qualityOpt = 'MUA';
% end
Y = get_spike_times(sorted_root,qualityOpt);

dir_info = dir([grid_root '/*Info.mat']); %Get the names of all files with this extension
grid_filename = fullfile(dir_info.folder, dir_info.name); %Form the name of the file to be opened
load(grid_filename);


trig_length_ms = grid.stimGrid(1,2); %The minimum length of one trigger in samples
trig_length_samples = (trig_length_ms/1000)*fs;
dB_lvls = unique(grid.stimGrid(:,3));
dB_lvls = sort(dB_lvls,'descend');
num_dB_lvls = numel(dB_lvls); %Find the number of different dB levels used
freqs = unique(grid.stimGrid(:,1)); %Find the specific frequencies used
num_freqs = numel(freqs); %Find the number of different frequencies used
num_stim = grid.nStimConditions; %Totla number of stimuli
t_ms = [0:t_bin_ms:sweep_duration_ms]; %The edges of the histogram
sum_window_bin = sum_window_ms/t_bin_ms; %Specifies which bins to extract from the psth for generating the fra

diff_sig = diff(synch_ch); %Find the difference between every n+1 sample - n sample. This tells us the the beginning/end of each sweep

%Find the sample no for the beginning of each sweep
start_ix = find(diff_sig>0); %Note that this diff will differ depending on which sync channel we use so always check it before analyzing
end_ix = find(diff_sig<0); %Note that this diff will differ depending on which sync channel we use so always check it before analyzing
start_ix = start_ix(1:length(end_ix));
diff_ix = end_ix - start_ix; %Find the length of each sweep in samples
start_ix = start_ix(diff_ix >= trig_length_samples); %Keep only the triggers which have length >= minimum triger length
start_ix_ms = (start_ix/fs).*1000; %Convert the starting sample numbers to times in ms

num_triggers = numel(start_ix_ms);

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
psth_final = reshape(psth_final,num_dB_lvls,num_freqs,total_no_clusters,length(t_ms)-1); %Reshape the mean_psth into an fra; 
% fra = sum(psth_final(:,:,:,1:sum_window_bin),4);
fra = sum(psth_final(:,:,:,20:34 ),4); %Alex 18/12/17
fra = flipud(fra); %Flip the fra so the y axis will be from loud to more quiet levels goung up -> down