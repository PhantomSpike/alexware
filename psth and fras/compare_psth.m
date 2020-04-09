% function ch_sp_rate = compare_psth(data_root,start_time,end_time)
data_root = '/data/alex/Data/PSTH_test/P01-repeated-noise-6';
start_time = 10;
end_time = 133;
fs = 30000;
total_chan = 384;
[spikeTimes,synch_ch] = get_spikes(data_root,start_time,end_time);
diff_sig = diff(synch_ch);
time_bin = 62.5; %Time bin in ms which is 125ms/2 i.e. the duration of the noise burst. This way we have the the onset and offset response being 1/2 of the total time
%Find the sample no for the beginning and end of each sweep
start_ix = find(diff_sig==1);
end_ix = find(diff_sig==-1);
%Convert the start and end samples into time in ms
start_time = (start_ix./fs).*1000; 
end_time = (end_ix./fs).*1000;
no_pres = length(start_ix);
onset_bins = [1:4:557]; %Assuming a time bin of 62.5ms there will be 540 bins in a 35s interval. We want the beignning of every 125ms sweep and the end.
offset_bins = [2:4:558];

for sweep = 1:no_pres
    edges(sweep,:) = [start_time(sweep):time_bin:end_time(sweep)]; %Make the edges of the histogram
end

for chan = 1:total_chan
    spike_times = spikeTimes{chan};
    for sweep_no = 1:no_pres
    [spike_count{chan}(sweep_no,:)] = histcounts(spike_times,edges(sweep_no,:));
    end
    mean_count_sweep = mean(spike_count{chan});
    mean_count_onset = mean(mean_count_sweep(onset_bins));
    spike_rate_onset = (mean_count_onset./time_bin).*1000; %Find the spike rate per second for the onset response
    mean_count_offset = mean(mean_count_sweep(offset_bins));
    spike_rate_offset = (mean_count_offset./time_bin).*1000; %Find the spike rate per second for the offset response
    ch_sp_rate(chan,1) = spike_rate_onset;
    ch_sp_rate(chan,2) = spike_rate_offset;
    ch_sp_rate(chan,3) = spike_rate_onset/spike_rate_offset;
end