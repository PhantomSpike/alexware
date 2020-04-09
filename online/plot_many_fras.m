clear all; close all;
pic_save_dir = '/media/alex/5FC39EAD5A6AA312/MSc_students/280219/Plots/';
default_bin_path = '/media/alex/5FC39EAD5A6AA312/MSc_students/Mouse14B_26_02_19/';
default_grid_path =  '/media/alex/5FC39EAD5A6AA312/MSc_students/expt25/';
plot_psths = true;


start_window_ms = 5; %The start of the summation window for generating FRAs in ms
end_window_ms = 60; %The end of the summation window for generating FRAs in ms
t_bin_ms = 5; %The time bin for the histogram in ms

no_channels = 25;
start_time_s = 1;

starting_channel = input('Starting channel for analysis: ');
end_channel = input('End channel for analysis: ');
end_time_s = input('End time [s]: ');
[bin_filename,bin_path] = uigetfile('*.ap.bin','Select the bilateral noise .ap.bin file',default_bin_path);
[grid_filename,grid_path] = uigetfile('*Info.mat','Select the bilateral noise benware gridInfo file',default_grid_path);
meta_dir = dir(fullfile(bin_path, '*ap*.meta')); % meta file from spikeGLX specifically
meta_info = readSpikeGLXmeta(fullfile(meta_dir.folder, meta_dir.name));
probe = meta_info.imProbeOpt;
fs = meta_info.sRateHz;

file_name = bin_filename(1:end-12);


% end_time_s = 600;
method = 'b';

channels_per_side = sqrt(no_channels);
num_figs = channels_per_side;
bin_fullname = fullfile(bin_path, bin_filename);
channels = round(linspace(starting_channel,end_channel,no_channels));


grid_fullname = fullfile(grid_path,grid_filename);
load(grid_fullname);

trig_length_ms = grid.stimGrid(1,2); %The minimum length of one trigger in samples
min_trig_length_s = (trig_length_ms/1000);
dB_lvls = unique(grid.stimGrid(:,3));
dB_lvls = sort(dB_lvls,'descend');
num_dB_lvls = numel(dB_lvls); %Find the number of different dB levels used
num_stim = grid.nStimConditions; %Total number of stimuli
freqs = unique(grid.stimGrid(:,1)); %Find the specific frequencies used
num_freqs = numel(freqs); %Find the number of different frequencies used
total_time_s = end_time_s - start_time_s;
min_inter_trig_length_s = grid.postStimSilence(1);

fprintf('== Processing channels #%0.f:%0.f ==\n',starting_channel,end_channel);
tic;

%Get the synch channel
[synch_ch] = get_synch_online(bin_fullname,start_time_s,end_time_s,probe);
%Get the triggers
[start_ix_ms] = get_triggers_new(synch_ch,min_trig_length_s,min_inter_trig_length_s,fs);

%Initialize variables
big_fra = cell(no_channels,1);
NPSP = zeros(no_channels,1);
big_psth_f_t = cell(no_channels,1);
% t_ms = zeros(no_channels,2);

parfor ch = 1:numel(channels)
    channel_no = channels(ch);
    fprintf('== Processing channel #%0.f ==\n',channel_no);
    [big_fra{ch},big_psth_f_t{ch},t_ms(ch,:)] = plot_fra(bin_fullname,grid_fullname,start_ix_ms,channel_no,start_time_s,end_time_s,method,probe,start_window_ms,end_window_ms,t_bin_ms);
end
fprintf('== Done! Processing took %0.fs ==\n',toc);

big_psth_f_t = cell2mat(big_psth_f_t);
big_psth_f_t = big_psth_f_t*(1000/t_bin_ms); %Convert to Spikes/s [Hz]
%% Plot the FRAs
fprintf('== Plotting the results ==\n');
tic;
row = channels_per_side;
col = channels_per_side;
per = 0.005;
edgel = 0.03; edger = per; edgeh = per; edgeb = 0.05; space_h = per; space_v = per;
[pos]=subplot_pos(row,col,edgel,edger,edgeh,edgeb,space_h,space_v);

figure('units','normalized','outerposition',[0 0 1 1]);
for ii = 1:numel(channels)
    subplot('position',pos{ii});
    min_fra = min(big_fra{ii}(:));
    max_fra = max(big_fra{ii}(:));
    imagesc((big_fra{ii} - min_fra)./(max_fra - min_fra));
    colormap('jet');
    axis off;
    if ii == col*(row - 1) + 1
        axis on;
        xticks([1:num_freqs]);
        freqs = ceil(freqs)/1000; %Convert the frequencies into kHz
        
        for jj = 1:numel(freqs)
            x_labels{jj} = num2str(freqs(jj),'%.1f');
        end
        xticklabels(x_labels);
        yticks(1:num_dB_lvls);
        y_labels = string(dB_lvls);
        yticklabels(y_labels);
        set(gca,'FontName','Arial','FontSize',8,'FontWeight','Bold');
        xlabel('Frequency [kHz]','FontSize',10,'FontWeight','bold');
        ylabel('Sound Level [dB]','FontSize',16,'FontWeight','bold');
    end
end
save_name = [pic_save_dir,file_name,'.png'];
export_fig(save_name);
close all;

%% Plot the FRA psths

if plot_psths
    num_subplots = channels_per_side*num_freqs;
    ch_ix = [1:channels_per_side:no_channels];
    select_ix = [1:num_subplots:size(big_psth_f_t,1)];
    row = channels_per_side;
    col = num_freqs;
    per = 0.005;
    edgel = 0.03; edger = per; edgeh = per; edgeb = 0.07; space_h = per; space_v = per;
    [pos]=subplot_pos(row,col,edgel,edger,edgeh,edgeb,space_h,space_v);
    
    for k = 1:num_figs
        figure('units','normalized','outerposition',[0 0 1 1]);
        name_ix = channels(ch_ix(k):ch_ix(k)+channels_per_side-1);
        plot_ix = [select_ix(k):select_ix(k)+num_subplots-1];
        psth_max = max(max(big_psth_f_t(plot_ix,:)));
        count = 0;
        for jj = 1:num_subplots
            subplot('position',pos{jj});
            ix = plot_ix(jj);
            bar(t_ms(1,1:end-1),big_psth_f_t(ix,:));
            axis off;
            ylim([0 psth_max]);
            sound_on_line = line(zeros(2,1),[0,psth_max],'Color','r');
            sound_off_line = line(trig_length_ms*ones(2,1),[0,psth_max],'Color','r');
            
            if jj >= col*(row - 1) + 1 && jj <= col*row
                count = count + 1;
                axis on;
                if jj ==col*(row - 1) + 1
                    ylabel('Spike rate [Hz]');
                else
                    set(gca,'ytick',[]);
                end
                xlabel([num2str(freqs(count),'%.0f'),'kHz']);
                set(gca,'FontName','Arial','FontSize',10,'FontWeight','Bold');
            end
        end
        save_name = [pic_save_dir,file_name,'_frapsth_Channels:',num2str(name_ix),'_.png'];
        export_fig(save_name);
        close all;
    end
end
fprintf('== Done! Plotting took %0.fs ==\n',toc);