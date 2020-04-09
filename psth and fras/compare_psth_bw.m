function [mean_spikes,t_ms] = compare_psth_bw(data_root,start_time,end_time)
plot_name = 'repeated-noise-6';
plot_name = [plot_name,'_',num2str(start_time),'s to ',num2str(end_time),'s'];
%A sweep is defined as one square wave in the 385th sync channel. This will
%contain n number of stimulus presentations. E.g. the probe noise sweep
%contains 39 stimulus presentations of 1s each. This 1s can be further divided into 500ms noise + 500ms silence

stimulus_period_ms = 1000; %The duration in ms of one stimulus presentation within a sweep
sweep_duration = 39; %The duration of one benware sweep in s


fs = 30000; %Sampling rate
total_chan = 384;

sweep_duration_ms = sweep_duration.*1000; %Convert to ms
n_cycles = sweep_duration_ms/stimulus_period_ms; %Find the number of stumulus presentations (cycles) in one sweep

[spikeTimes,synch_ch] = get_spikes(data_root,start_time,end_time);

diff_sig = diff(synch_ch); %Find the difference between every n+1 sample - n sample. This tells us the the beginning/end of each sweep

%Find the sample no for the beginning of each sweep
start_ix = find(diff_sig==1);
start_ix_ms = (start_ix/fs).*1000; %Convert the starting sample numbers to times in ms

stimulus_offset_ms = [0:n_cycles-1]*stimulus_period_ms; %Define the beginning of every stimulus presentation within a sweep

start_time_ms = [];

%Get all start times of every stimulus presentation in every sweep and
%concatenate them together in one big vector
for sweep_no = 1:length(start_ix_ms)
    start_time_ms = [start_time_ms start_ix_ms(sweep_no) + stimulus_offset_ms];
end

%Check whether the last stimulus presentation is interupted and discard it
%in thath case
if start_time_ms(end)+stimulus_period_ms > end_time*1000
    start_time_ms = start_time_ms(1:end-1);
end

binsize_ms = 5;
t_ms = 0:binsize_ms:stimulus_period_ms;
no_bins = stimulus_period_ms/5;

%Extract the spike counts for every stimulus presentation for every channel
% in a histogram with specified bins
spikes = [];
for stim_pres_no = 1:length(start_time_ms)
    for chan_no = 1:length(spikeTimes)
        spikes(:, stim_pres_no, chan_no) = histc(spikeTimes{chan_no}, start_time_ms(stim_pres_no)+t_ms);
    end
end
spikes = spikes(1:end-1,:, :); %Get rid of the last bin which is weird
mean_spikes = mean(spikes,2); %Take the mean across all stimulus presentations
mean_spikes = reshape(mean_spikes,no_bins,total_chan);
mean_spikes = mean_spikes.*(1000/binsize_ms); %Multiply to convert spike count into spikes/s
%% Discriminate between real neurons and noise
noise_duration_ms = 500; %This is the duration of the noise in one stimulus presentation
search_interval_ms = 50; %This is the duration of the interval in ms where we look for the maximum
binsize_ms = 5;

search_interal_bins = search_interval_ms/binsize_ms; %This is the search interval in bins where we look for the maximum
max_resp = max(mean_spikes(1:search_interal_bins,:)); %Find the maximum value of spike count in the search interval
noise_duration_bins = noise_duration_ms/binsize_ms; %This is the noise duration in bins
mean_val = mean(mean_spikes(1:noise_duration_bins,:)); %Find the mean value of the signal in the period of noise duration
reality_ix = (max_resp - mean_val)./mean_val; %Generate an index which will compare the first peak response to the overall signal
%% Plotting 1
save_dir = '/home/phant0msp1ke/Desktop/PSTHs_Ben/';
figure('units','normalized','outerposition',[0 0 1 1]);
channels = [1:384];
plot(channels,reality_ix,'r*');
title('Reality index for all channels');
xlabel('Channel #');
ylabel('(Max initial resp - Mean)/Mean');
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',20);
hold on;
plot(channels,reality_ix);
hold off;
save_name = [save_dir,plot_name,'_reality_ix_allchannels.jpg'];
export_fig(save_name,'-q100');
%% Plotting 2
channel_no = 304;
figure('units','normalized','outerposition',[0 0 1 1]);
bar(t_ms(1:end-1), mean_spikes(:,channel_no), 'hist');
title(['Histogram of channel no ',num2str(channel_no)]);
xlabel('Time [ms]');
ylabel('Spikes/s');
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',20);
save_name = [save_dir,plot_name,' channel ',num2str(channel_no),'.png'];
export_fig(save_name);
%% Plotting 3
per = 0.025;
edgel = per; edger = per; edgeh = per; edgeb =per; space_h = 0.025; space_v =0.05;
[pos]=subplot_pos(5,8,edgel,edger,edgeh,edgeb,space_h,space_v);
chan_list = [0:40:320];
for chunk = 1:9
    figure('units','normalized','outerposition',[0 0 1 1]);
    chan = chan_list(chunk);
    for ii = 1:40
        subplot('position',pos{ii});
        bar(t_ms(1:end-1), mean_spikes(:,chan + ii), 'hist');
        legend(['Channel ', num2str(chan + ii)],'Location','northeast');
    end
    save_name = [save_dir,plot_name,' Chunk #',num2str(chunk),'.png'];
    export_fig(save_name);
    close(gcf);
end
