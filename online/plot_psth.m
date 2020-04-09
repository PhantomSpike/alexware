function [psth_final,raster,t_ms,t_ms_r,NPSP,total_spikes,edges_ms,psth_check] = plot_psth(bin_fullname,start_ix_ms,sweep_duration_ms,extra_time_ms,channel_no,start_time_s,end_time_s,method,probe,t_bin_ms)

[spike_times_ms] = extract_spikes_pixels(bin_fullname,channel_no,start_time_s,end_time_s,method,probe);


t_bin_ms_r = 0.1; %Time bin for the raster plot
edges_ms = [t_bin_ms:t_bin_ms:(end_time_s - start_time_s)*1000];
psth_check = histc(spike_times_ms,edges_ms);
total_spikes = numel(spike_times_ms);

t_ms = [-extra_time_ms:t_bin_ms:sweep_duration_ms]; %The edges of the histogram
t_ms_r = [-extra_time_ms:t_bin_ms_r:sweep_duration_ms]; %The edges of the histogram for raster plot

num_triggers = length(start_ix_ms);


%Find the psths for the given channel

psth = zeros(num_triggers,numel(t_ms));
psth_r = zeros(num_triggers,numel(t_ms_r));

for trigger = 1:num_triggers
    psth(trigger,:) = histc(spike_times_ms,start_ix_ms(trigger) + t_ms);
    psth_r(trigger,:) = histc(spike_times_ms,start_ix_ms(trigger) + t_ms_r);
end

psth = psth(:,1:end-1); %Delete the last bin which is weird
[SP, NP, TP, SP_std_error] = sahani_quick(psth);
NPSP = NP/SP;
psth_r = psth_r(:,1:end-1); %Delete the last bin which is weird
avg_psth = mean(psth,1); %Find the average across repetitions
psth_final = avg_psth.*(1000/t_bin_ms); %Convert firing rate into spikes/s

for trial = 1:size(psth_r,1)
    raster(trial,:) = trial*psth_r(trial,:);
end