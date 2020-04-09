%% Define parameters
%Give the correct directory with Kilosorted files
sorted_dir = '/mnt/40086D4C086D41D0/Reverb_data/Sorted/P08_Derry/P08-alex_reverb_wtih_noise_same.3';
%Give the correct Benware Grid folder
grid_root =     '/mnt/40086D4C086D41D0/Reverb_data/Meta_data/benware/Derry_P08'; 
%Give the synch_ch file
synch_path =   '/mnt/40086D4C086D41D0/Reverb_data/Meta_data/synch_ch/Derry_P08_synch_ch.mat';

min_trig_length_s = 40; %The minmum trigger length in seconds. This will depend on the stimulus. It is best practice to check the synch channel always
min_inter_trig_length_s = 0.1; %The minmum inter trigger length in seconds. This will depend on the stimulus. It is best practice to check the synch channel always
actual_stimlength_s = 36; %The length of each presentation without a noise burst in s
num_stim = 6; 
num_repeats = 10; %Number of different stimuli. They are 6 of them (3 reverb conditions and 2 sound stim sets) with 10 repeat each. Therefore the indices are like this
%1 – anechoic stim 1 -> [1:10]
%2 – anechoic stim 2 -> [11:20]
%3 – small reverb stim 1 -> [21:30]
%4 – small reverb stim 2 -> [31:40]
%5 – big reverb stim 1 -> [41:50]
%6 – big reverb stim 2 -> [51:60]
%% Get the clust_info
load(fullfile(sorted_dir,'clust_info'));
%% Get the synch channel
load(synch_path);
%% Get the triggers
fs = clust_info.fs;
[start_time_ms] = get_triggers(synch_ch,min_trig_length_s,min_inter_trig_length_s,fs);
start_time_s = start_time_ms/1000; % Convert to s
%% Load the benware file
dir_info = dir([grid_root '/*Info.mat']); %Get the names of all files with this extension
grid_filename = fullfile(dir_info.folder, dir_info.name); %Form the name of the file to be opened
ben_grid = load(grid_filename);
%% Loop through every cluster and save the spiketimes for the particular stim and repeat
num_clust = length(clust_info.clust_id); %Find the total number of clusters
max_clust = max(clust_info.clust_id);
data(max_clust).stim(num_stim).repeat(num_repeats).spiketimes = []; %Initialize the struct
fprintf('== Getting clusters ==\n');tic;
for c = 1:num_clust
    fprintf('== Processing cluster%0.f/%0.f ==\n',c,num_clust);
    stim_ix = 0; %Counter var to keep track of which index we are using
    clust_id = clust_info.clust_id(c)+1; %Get the clust id for that neuron. We add one to be able to index;
    spiketimes = clust_info.spikeTimes(clust_info.spikeTimes(:,2)==clust_id,1); %Get the spiketimes that belong to it
    
    for s = 1:num_stim
        
        for r = 1:num_repeats
            stim_ix = stim_ix + 1;
            ix_order = find(ben_grid.grid.randomisedGrid==stim_ix); %Find in which order was this particular stim and rep played
            %Get only the spiketimes that correspond to this stimulus and
            %this repeat
            spiketimes_rep = spiketimes(spiketimes>=start_time_s(ix_order) & spiketimes<(start_time_s(ix_order)+actual_stimlength_s));
            %Normalize relative to the beginning
            data(clust_id).stim(s).repeat(r).spiketimes = spiketimes_rep - start_time_s(ix_order);
        end
        
    end
    
end
fprintf('== Done! This took %.1fs ==\n',toc);
save(fullfile(sorted_dir,'Derry_P08_data.mat'),'data');