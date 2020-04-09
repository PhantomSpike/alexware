function [fra,psth_f_t,t_ms] = plot_fra(bin_fullname,grid_fullname,start_ix_ms,channel_no,start_time_s,end_time_s,method,probe,start_window_ms,end_window_ms,t_bin_ms)


start_hist_ms = -25; %The start of the histogram in ms from the trigger
extra_time_ms = 70; %How many ms to add after the sweep is done


[spike_times_ms] = extract_spikes_pixels(bin_fullname,channel_no,start_time_s,end_time_s,method,probe);

load(grid_fullname);
% total_time_s = end_time_s - start_time_s;
sweep_duration_ms = grid.stimGrid(1,2); %The duration of one time interval of interest in ms. Actual Benware sweep duration is around 105.2ms


dB_lvls = unique(grid.stimGrid(:,3));
dB_lvls = sort(dB_lvls,'descend');
num_dB_lvls = numel(dB_lvls); %Find the number of different dB levels used
freqs = unique(grid.stimGrid(:,1)); %Find the specific frequencies used
num_freqs = numel(freqs); %Find the number of different frequencies used
num_stim = grid.nStimConditions; %Total number of stimuli

t_ms = [start_hist_ms:t_bin_ms:(sweep_duration_ms + extra_time_ms)]; %The edges of the histogram
end_window_bin = round(end_window_ms/t_bin_ms + abs(start_hist_ms/t_bin_ms) + 1); %Specifies which bins to extract from the psth for generating the fra
start_window_bin = round(start_window_ms/t_bin_ms + abs(start_hist_ms/t_bin_ms) + 1); 
% start_pre_window_bin = round(start_hist_ms/t_bin_ms + abs(start_hist_ms/t_bin_ms) + 1);
% end_pre_window_ms = -5;
% end_pre_window_bin = round(end_pre_window_ms/t_bin_ms + abs(start_hist_ms/t_bin_ms) + 1);


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
psth_temp = cell2mat(avg_psth);
% NPSP_all = cellfun(@(x)(sahani_quick_online(x)),psth,'UniformOutput',false);
% NPSP = nanmin(cell2mat(NPSP_all));
psth_l_f_t = reshape(psth_temp,num_dB_lvls,num_freqs,numel(t_ms)-1); %Reshape to convert to dimensions dB_lvls*freq*time bins 
psth_f_t = squeeze(mean(psth_l_f_t,1)); %Find the mean psth for every freqeuncy across levels
fra = sum(psth_l_f_t(:,:,start_window_bin:end_window_bin),3); %Find the mean number of spikes  in the first n ms and subtract them from the pre-stimulus activity 
fra = flipud(fra); %Flip the fra so the y axis will be from loud to more quiet levels goung up -> down
end
