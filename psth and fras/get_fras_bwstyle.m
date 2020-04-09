function [mean_psth,fras] = get_fras_bwstyle(data_root,grid,start_time_s,end_time_s,chunk_size_s)

sweep_duration_ms = 200; %The duration of one time interval of interest in ms. Actual Benware sweep duration is around 105.2ms
t_bin_ms = 5; %Time bin in ms
fs_s = 30000; %Sampling rate in samples/s
total_chan = 384; %Total number of channels
num_dB_levels = 5; %The number of different dB levels that were presented
num_freq = 24; %The number of different freq
num_stim = num_dB_levels*num_freq; %Totla number of stimuli  
t_ms = [0:t_bin_ms:sweep_duration_ms]; %The edges of the histogram

total_time = end_time_s - start_time_s; %Find the total time in sec that will be extracted
num_chunks = total_time./chunk_size_s; %Find the number of chunks necessary to load the data

if ~(floor(num_chunks)==num_chunks)
    error('The total time in seconds should be a multiple of the chunk size in seconds');
end

%Define the start and end times of all chunks
chunk_start_times_s = [start_time_s:chunk_size_s:end_time_s - chunk_size_s];
chunk_end_times_s = chunk_start_times_s + chunk_size_s;

chunk_psth = zeros(num_stim,total_chan,length(t_ms)-1,num_chunks); %Prealocate variable 

for chunk = 1:num_chunks
    fprintf('== Processing chunk %.0f/%.0f ==\n',chunk,num_chunks);tic;
    [spikeTimes,synch_ch] = get_spikes(data_root,chunk_start_times_s(chunk),chunk_end_times_s(chunk));
    
    diff_sig = diff(synch_ch); %Find the difference between every n+1 sample - n sample. This tells us the the beginning/end of each sweep
    
    %Find the sample no for the beginning of each sweep
    start_ix = find(diff_sig==1);
    
    %Check whether the beginning of the chunk coincides with the beginning
    %of a sweep. If so, add the first sample to the start ix
    if synch_ch(1) == -1
        warning('Chunk no %.0f begins with a sweep start. Investigate fruther',chunk);
        start_ix = [1,start_ix];
    end
    
    start_ix_ms = (start_ix/fs_s).*1000; %Convert the starting sample numbers to times in ms
    
    %Check if the last sweep is interupted before it finishes. If so, delete
    %the last starting time
    if (start_ix_ms(end) + sweep_duration_ms)>chunk_end_times_s(chunk)*1000
        warning('Chunk no %.0f ends before the sweep has finished. Investigate further',chunk);
        start_ix_ms(end) = [];
    end
    
    %Note the last grid set index which was used and shift by one for the
    %next chunk so that the correct indices are used
    if chunk == 1
        next_grid_set_ix(1) = 1;
        last_grid_set_ix(1) =  length(start_ix_ms);
    else
        next_grid_set_ix(chunk) = last_grid_set_ix(chunk-1) + 1;
        last_grid_set_ix(chunk) = next_grid_set_ix(chunk) + length(start_ix_ms) - 1;
    end
    
    
    psth = cell(num_stim,total_chan); %Initialize the psth variable
    
    %Find the psths for every channel, stimulus and repetition and store
    %them in a cell array
    for ch = 1:total_chan
        for stim = 1:num_stim
            ix_rep = find(grid.randomisedGridSetIdx(next_grid_set_ix(chunk):last_grid_set_ix(chunk),1)==stim);
            ix_rep_ms = start_ix_ms(ix_rep);
            for rep = 1:length(ix_rep)
                psth{stim,ch}(rep,:) = histc(spikeTimes{ch},ix_rep_ms(rep) + t_ms);
            end
            psth{stim,ch} = psth{stim,ch}(:,1:end-1); %Delete the last bin which is weird
        end
    end
    avg_psth = cellfun(@(x)(mean(x,1)),psth,'UniformOutput',false); %Find the average across repetitions
    %Convert the cell array into 4D array with dimesnions stimuli x
    %channels x time bins x chunks
    int_mat = cellfun(@(x)reshape(x,1,1,[]),avg_psth,'un',0);
    chunk_psth(:,:,:,chunk) = cell2mat(int_mat);
    fprintf('== Processing chunk %.0f/%.0f took %.1f sec ==\n',chunk,num_chunks,toc);
end
mean_psth = mean(chunk_psth,4); %Average across all chunks
fras = reshape(mean_psth,num_dB_levels,num_freq,total_chan,length(t_ms)-1); %Reshape the mean_psth into an fra
end

