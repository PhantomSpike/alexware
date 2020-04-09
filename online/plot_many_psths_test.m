clear all; close all;
pic_save_dir = '/home/alex/Desktop/Troubleshoot_pics/';
default_bin_path = '/media/alex/Alex SSD/Cork_24_07_2018/';
default_grid_path = '/media/alex/Alex SSD/Cork_24_07_2018/benware_gridinfo/';
[bin_filename,bin_path] = uigetfile('*.ap.bin','Select the bilateral noise .ap.bin file',default_bin_path);
[grid_filename,grid_path] = uigetfile('*Info.mat','Select the bilateral noise benware gridInfo file',default_grid_path);
probe = input('Select the probe option [3 or 4]: ');
file_name = bin_filename(1:end-12);
% no_channels = input('Select number of channels [9, 16 or 25]: ');
no_channels = 25;
switch probe
    case 3
        %         starting_channel = input('Starting channel [1 - 359]: ');
        %         end_channel = input('Ending channel [1 - 384]: ');
        starting_channel = 10;
        end_channel = 300;
    case 4
%         starting_channel = input('Starting channel [1 - 251]: ');
%         end_channel = input('Ending channel [1 - 276]: ');
        starting_channel = 5;
        end_channel = 270;
end


channels_per_side = sqrt(no_channels);
method = 'b'; %Method for detecting the spikes
fs = 30000; %The sampling rate of the Neuropixels

extra_time_ms = 200;
grid_fullname = fullfile(grid_path,grid_filename);
load(grid_fullname);
duration_ms = grid.stimGrid(1);
begin_left_ms = grid.stimGrid(2);
begin_right_ms = grid.stimGrid(3);
begin_both_ms = grid.stimGrid(4);
sweep_duration_ms = begin_both_ms + duration_ms + extra_time_ms; %The duration of one time interval of interest in ms. Actual duration can be slightly different

% start_time_s = input('Start time [s]: ');
start_time_s = 1;
% end_time_s = input('End time [s]: ');
end_time_s = (numel(dir([grid_path 'sweep.mat/'])) - 2)*(sweep_duration_ms/1000 + 0.14);
num_stim = numel(dir([grid_path 'sweep.mat/'])) - 12;

channels = round(linspace(starting_channel,end_channel,no_channels));
big_psth = cell(numel(channels),1);

bin_fullname = fullfile(bin_path, bin_filename);
fprintf('== Processing channels #%0.f:%0.f ==\n',starting_channel,end_channel);
tic;

total_time_s = end_time_s - start_time_s;
min_trig_length_s = (sweep_duration_ms - extra_time_ms)/1000;
min_inter_trig_length_s = grid.postStimSilence(1);

%Get the synch channel
[synch_ch] = get_synch_online(bin_fullname,start_time_s,end_time_s,probe);
%Get the triggers
% [start_ix_ms] = get_triggers_new(synch_ch,min_trig_length_s,min_inter_trig_length_s,fs);
sweep_dur_ms = 1553.5;
diff_sig = diff(synch_ch); 
start_ix = find(diff_sig>0);
first_start_ix = start_ix(1);
first_start_ms = (first_start_ix/fs)*1000;
start_ix_ms = first_start_ms*ones(1,num_stim) + [0:1:num_stim-1]*sweep_dur_ms;

parfor ch = 1:numel(channels)
    channel_no = channels(ch);
    fprintf('== Processing channel #%0.f ==\n',channel_no);
    [big_psth{ch},big_raster{ch},t_ms(ch,:),t_ms_r(ch,:),NPSP(ch)] =  plot_psth(bin_fullname,start_ix_ms,sweep_duration_ms,extra_time_ms,channel_no,start_time_s,end_time_s,method,probe);
end

fprintf('== Done! Processing took %0.fs ==\n',toc);

row = channels_per_side;
col = channels_per_side;
per = 0.005;
edgel = 0.035; edger = per; edgeh = per; edgeb = 0.05; space_h = per; space_v = 0.01;
[pos]=subplot_pos(row,col,edgel,edger,edgeh,edgeb,space_h,space_v);

%% Plot PSTH
figure('units','normalized','outerposition',[0 0 1 1]);
for ii = 1:numel(channels)
    subplot('position',pos{ii});
    bar(t_ms(1,1:end-1), big_psth{ii}, 'hist');
    xlim([t_ms(1,1) t_ms(1,end-1)]);
    psth_max = max(big_psth{ii}(:));
    left_line_on = line(zeros(2,1),[0,psth_max],'Color','b');
    left_line_off = line(duration_ms*ones(2,1),[0,psth_max],'Color','b');
    right_line_on = line(begin_right_ms*ones(2,1),[0,psth_max],'Color','r');
    right_line_off = line((begin_right_ms+duration_ms)*ones(2,1),[0,psth_max],'Color','r');
    both_line_on = line(begin_both_ms*ones(2,1),[0,psth_max],'Color','m');
    both_line_off = line((begin_both_ms+duration_ms)*ones(2,1),[0,psth_max],'Color','m');
    axis off;
    npsp_val = num2str(NPSP(ii));
    legend(['NPSP: ',npsp_val,' Ch ',num2str(channels(ii)),' ' ,num2str(channels(ii)*10),'um']);
    if ii == col*(row - 1) + 1
        axis on;
        xlabel('Time [ms]');
        ylabel('Rate [spikes/s]');
        set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold');
    end
end

save_name = [pic_save_dir,file_name,'_psth.png'];
export_fig(save_name);
close all;

%% Plotting Raster
row = channels_per_side;
col = channels_per_side;
per = 0.005;
edgel = 0.03; edger = per; edgeh = per; edgeb = 0.05; space_h = 0.0075; space_v = 0.01;
[pos]=subplot_pos(row,col,edgel,edger,edgeh,edgeb,space_h,space_v);
figure('units','normalized','outerposition',[0 0 1 1]);
for jj = 1:numel(channels)
    subplot('position',pos{jj});
    plot(t_ms_r(1,1:end-1),big_raster{jj},'.k');
    axis off;
    raster_max = size(big_raster{jj},1);
    left_line_r_on = line(zeros(2,1),[0,raster_max],'Color','b');
    left_line_r_off = line(duration_ms*ones(2,1),[0,raster_max],'Color','b');
    right_line_r_on = line(begin_right_ms*ones(2,1),[0,raster_max],'Color','r');
    right_line_r_off = line((begin_right_ms+duration_ms)*ones(2,1),[0,raster_max],'Color','r');
    both_line_r_on = line(begin_both_ms*ones(2,1),[0,raster_max],'Color','m');
    both_line_r_off = line((begin_both_ms+duration_ms)*ones(2,1),[0,raster_max],'Color','m');
    ylim([1 raster_max]);
    xlim([t_ms_r(1,1) t_ms_r(1,end-1)]);
    if jj == col*(row - 1) + 1
        axis on;
        xlabel('Time [ms]');
        ylabel('Trial number');
        set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold');
    end
    
end
save_name2 = [pic_save_dir,file_name,'_raster.png'];
export_fig(save_name2);
close all;