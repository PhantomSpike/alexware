clear all; close all;
data_root = input('Data root: ');
grid_root = input('Benware grid root: ');
channel_no = input('Choose channel [1 - 384]: ');
start_time_s = input('Start time [s]: ');
end_time_s = input('End time [s]: ');
method = input('Spike extraction method (default is benware): ');
if isempty(method)
    method = 'b';
end

[spike_times_ms,synch_ch,fs] = extract_spikes_pixels_ind(data_root,channel_no,start_time_s,end_time_s,method);


sweep_duration_ms = 200; %The duration of one time interval of interest in ms. Actual Benware sweep duration is around 105.2ms
t_bin_ms = 5; %Time bin in ms
sum_window_ms = 80; %The summation window for generating FRAs in ms

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
start_ix = find(diff_sig==1); %Note that this diff will differ depending on which sync channel we use so always check it before analyzing
end_ix = find(diff_sig==-1); %Note that this diff will differ depending on which sync channel we use so always check it before analyzing
start_ix = start_ix(1:numel(end_ix)); %Remove triggers that did not finish
diff_ix = end_ix - start_ix; %Find the length of each sweep in samples
start_ix = start_ix(diff_ix >= trig_length_samples); %Keep only the triggers which have length >= minimum triger length
start_ix_ms = (start_ix/fs).*1000; %Convert the starting sample numbers to times in ms

num_triggers = numel(start_ix_ms);

psth = cell(num_stim,1); %Initialize the psth variable

%Find the psths for every cluster, stimulus and repetition and store
%them in a cell array

    for stim = 1:num_stim
        ix_rep = find(grid.randomisedGridSetIdx(1:num_triggers,1)==stim);
        ix_rep_ms = start_ix_ms(ix_rep);
        for rep = 1:length(ix_rep)
            psth{stim}(rep,:) = histc(spike_times_ms,ix_rep_ms(rep) + t_ms); %The PSTH is a tensor of dimensions stimuli(dB_lvls*freq)*time bins*repetitions
        end
        psth{stim} = psth{stim,1}(:,1:end-1); %Delete the last bin which is weird
    end
    
    %Find the average across repetitions
    avg_psth = cellfun(@(x)(mean(x,1)),psth,'UniformOutput',false);
    %Convert the cell array into 2D array with dimesnions stimuli(dB_lvls*freq)*time bins
    psth_final = cell2mat(avg_psth);
    psth_final = reshape(psth_final,num_dB_lvls,num_freqs,numel(t_ms)-1); %Reshape to convert to dimensions dB_lvls*freq*time bins
    fra = sum(psth_final(:,:,1:sum_window_bin),3); %Sum the spikes in the first n bins
    fra = flipud(fra); %Flip the fra so the y axis will be from loud to more quiet levels goung up -> down
    
 %% Plotting the FRA
 figure;
 imagesc(fra);
 colormap('jet');
 colorbar;
 xlabel('Frequency [kHz]');
 ylabel('Sound Level [dB]');
 title(['FRA for channel #',num2str(channel_no)]);
 xticks([1:num_freqs]);
 freqs = ceil(freqs)/1000; %Convert the frequencies into kHz
 x_labels = string(freqs(1:num_freqs));
 xticklabels(x_labels);
 yticks(1:num_dB_lvls);
 y_labels = string(dB_lvls);
 yticklabels(y_labels);
 set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold');