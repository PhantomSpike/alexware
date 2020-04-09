function [spikeTimes,synch_chan] = get_spikes_old(data_root,start_time,end_time,plot_name)
plot_on = false;
fs = 30000; %In Hz
spikeThreshold = 3; %In standard deviations
%For average plotting
first_ch = 20;
last_ch = 120;
step_ch = 10;
chans = [first_ch:last_ch];
%For individual plotting
chan_list = [first_ch:step_ch:last_ch];
time_bin_size = 20; %The size of the time bin in ms for the PSTH

total_time = end_time - start_time;
time_bin_no = (1000.*total_time)./time_bin_size;

%Load the data
[cropped_data] = load_spikedata(data_root,start_time,end_time);

%Extract the sync channel so it is not filtered
synch_chan = cropped_data(385,:);
cropped_data(385,:) = [];
spikeTimes = find_spikes(cropped_data,fs,spikeThreshold);

%% Plotting
if plot_on
    save_dir = '/home/phant0msp1ke/Desktop/PSTHs_Ben/';
    if nargin == 4
        save_name = [save_dir,plot_name,'_avg.png'];
        figure('units','normalized','outerposition',[0 0 1 1]);
    else
        figure;
    end
    
    per = 0.04;
    edgel = per; edger = per; edgeh = per; edgeb = per; space_h = 0.01; space_v =0.06;
    [pos]=subplot_pos(2,1,edgel,edger,edgeh,edgeb,space_h,space_v);
    subplot('position',pos{1});
    hist(cat(1, spikeTimes{chans}),time_bin_no) %This plots PSTH
    title(['PSTH for Channels ' ,num2str(chans(1)), ':', num2str(chans(end))]);
    xlabel('Time [ms]');
    ylabel('Spike count');
    subplot('position',pos{2});
    time_points= [1/fs:1/fs:total_time].*1000;
    plot(time_points,synch_chan);
    title('Syncrhonization channel');
    xlabel('Time [ms]');
    ylabel('Stimulus On/Off');
    if nargin == 4
        export_fig(save_name);
    end
    
    if nargin == 4
        figure('units','normalized','outerposition',[0 0 1 1]);
    else
        figure;
    end
    sz = length(chan_list);
    per = 0.01;
    edgel = per; edger = per; edgeh = per; edgeb =0.05; space_h = 0.01; space_v =0.01;
    [pos]=subplot_pos(sz+1,1,edgel,edger,edgeh,edgeb,space_h,space_v);
    for ii = 1:sz
        subplot('position',pos{ii});
        hist(cat(1, spikeTimes{chan_list(ii)}),time_bin_no); %This plots PSTH
        legend(['Channel ', num2str(chan_list(ii))],'Location','northwest');
        axis off;
    end
    subplot('position',pos{sz+1});
    time_points= [1/fs:1/fs:total_time].*1000;
    plot(time_points,synch_chan);
    %         title('Syncrhonization channel');
    xlabel('Time [ms]');
    if nargin == 4
        save_name = [save_dir,plot_name,'_ind.png'];
        export_fig(save_name);
    end
    ylabel('Signal');
end
