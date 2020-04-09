clear all; close all;
data_root = input('Data root: ');
channel_no = input('Choose a channel [1 - 384]: ');
start_time_s = input('Start time [s]: ');
end_time_s = input('End time [s]: ');
method = input('Spike extraction method (default is benware): ');
if isempty(method)
    method = 'b';
end

[spike_times_ms,synch_ch,fs] = extract_spikes_pixels_ind(data_root,channel_no,start_time_s,end_time_s,method);

t_bin_ms = 5; %Time bin in ms
t_bin_ms_r = 0.1; %Time bin for the raster plot

total_time_s = end_time_s - start_time_s;
extra_time_ms = 200;
sweep_duration_ms = 1200 + extra_time_ms; %The duration of one time interval of interest in ms. Actual duration can be slightly different
trig_min_length = (sweep_duration_ms - extra_time_ms)*(fs/1000);
t_ms = [-200:t_bin_ms:sweep_duration_ms]; %The edges of the histogram
t_ms_r = [-200:t_bin_ms_r:sweep_duration_ms]; %The edges of the histogram for raster plot

diff_sig = diff(synch_ch); %Find the difference between every n+1 sample - n sample. This tells us the the beginning/end of each sweep

%Find the sample no for the beginning of each sweep
start_ix = find(diff_sig==1); %Note that this diff will differ depending on which sync channel we use so always check it before analyzing
end_ix = find(diff_sig==-1); %Note that this diff will differ depending on which sync channel we use so always check it before analyzing
start_ix_ms = (start_ix/fs).*1000; %Convert the starting sample numbers to times in ms

if start_ix_ms(end) + sweep_duration_ms > total_time_s*1000
    start_ix(end) = [];
    start_ix_ms(end) = [];
end

diff_ix = end_ix - start_ix; %Find the length of each sweep in samples
start_ix_ms = start_ix_ms(diff_ix >= trig_min_length); %Keep only the triggers which have length >= minimum triger length
num_triggers = length(start_ix_ms);


%Find the psths for the given channel


for trigger = 1:num_triggers
    psth(trigger,:) = histc(spike_times_ms,start_ix_ms(trigger) + t_ms);
    psth_r(trigger,:) = histc(spike_times_ms,start_ix_ms(trigger) + t_ms_r);
end
psth = psth(:,1:end-1); %Delete the last bin which is weird
psth_r = psth_r(:,1:end-1); %Delete the last bin which is weird

avg_psth = mean(psth,1); %Find the average across repetitions

psth_final = avg_psth.*(1000/t_bin_ms); %Convert firing rate into spikes/s

%% Plotting Raster
figure;
for trial = 1:size(psth_r,1)
    raster(trial,:) = trial*psth_r(trial,:);
end
plot(t_ms_r(1:end-1),raster,'.k');
title(['Raster plot for channel #', num2str(channel_no),' bilateral noise L->R->Both']);
xlabel('Time [ms]');
ylabel('Trial number');
raster_max = max(raster(:));
left_line_r = line(zeros(2,1),[0,raster_max],'Color','r');
right_line_r = line(500*ones(2,1),[0,raster_max],'Color','r');
both_line_r = line(1000*ones(2,1),[0,raster_max],'Color','r');
ylim([1 raster_max]);
set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold');
%% Plotting PSTH
figure;
bar(t_ms(1:end-1), psth_final, 'hist');
xlim([t_ms(1) t_ms(end-1)]);
title(['PSTH for channel #', num2str(channel_no),' bilateral noise L->R->Both']);
xlabel('Time [ms]');
ylabel('Rate [spikes/s]');
psth_max = max(psth_final(:));
left_line = line(zeros(2,1),[0,psth_max],'Color','r');
right_line = line(500*ones(2,1),[0,psth_max],'Color','r');
both_line = line(1000*ones(2,1),[0,psth_max],'Color','r');
set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold');