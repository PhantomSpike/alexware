%Give the correct directory with Kilosorted files
data_root = '/media/alex/4TB SSD/DATA/Ronnie_23_01_2018/bin_files/Preprocessed/CRA_done/P13-reverb_with_noise_same_cra/P13-reverb_with_noise_same_cra_sorted';
%Give the correct Benware Grid folder
grid_root = '/media/alex/4TB SSD/DATA/Ronnie_23_01_2018/benware_gridinfo/P13-reverb_with_noise_same'; 
%Give the metadata file
meta_path = '/media/alex/4TB SSD/DATA/Ronnie_23_01_2018/metadata/reverb_ordered_metadata_same';
%Give the synch_ch file
synch_path = '/media/alex/4TB SSD/DATA/Ronnie_23_01_2018/metadata/synch_ch_reverbwithnoisesame_P13';

t_bin_ms = 10; %The time bin in ms that we use for the psths
first_bin_ms = -100; %The first time bin in ms relative to the stim onset
last_bin_ms = 600; %The last time bin in ms relative to the stim onset

qualityOpt = 'nonoise'; %Which untis you want to analyse

min_trig_length_s = 20; %The minmum trigger length in seconds. This will depend on the stimulus. It is best practice to check the synch channel always

fs = 30000; %Sampling rate in seconds

%% Get the noise position in the reverb stim
load(meta_path);
for stim_no = 1:numel(reverb_metadata) 
    noisepos(stim_no) = reverb_metadata(stim_no).noisepos;
    noisepos = noisepos';
end


%% Get the synch channel
load(synch_path);
% synch_ch = get_synch(data_root);

%% Get the triggers
[start_time_ms] = get_triggers(synch_ch,min_trig_length_s,fs);

%% Get the spike times
Y = get_spike_times(data_root,qualityOpt);

%% Get the noise positions and generate the PSTHs 

psth_sahani = psth_func_noise(Y,start_time_ms,first_bin_ms,last_bin_ms,t_bin_ms,noisepos,grid_root);

%% Run sahani_quick to get NPSP measures
no_clusters = size(psth_sahani,2);
psth_sahani(no_clusters).NPSP = [];
psth_sahani(no_clusters).SP = [];
psth_sahani(no_clusters).NP = [];
psth_sahani(no_clusters).TP = [];
psth_sahani(no_clusters).SP_std_error = [];

fprintf('== Calculating NPSPs ==\n');tic;
for cluster = 1:no_clusters
    fprintf('== Processing cluster %.0f/%.0f ==\n',cluster,no_clusters);
    [SP, NP, TP, SP_std_error] = sahani_quick(psth_sahani(cluster).data); 
    psth_sahani(cluster).NPSP = NP/SP;
    psth_sahani(cluster).SP = SP;
    psth_sahani(cluster).NP= NP;
    psth_sahani(cluster).TP = TP;
    psth_sahani(cluster).SP_std_error = SP_std_error;
end
fprintf('== Done! Processing took %.0fs ==\n',toc);