%This function gets the start times of stimuli, extracts the corresponding
%LFP and takes the average
%% Names of the files
lfp_folder = '/mnt/40086D4C086D41D0/Inhibitory_excitatory_analysis/For_analysis/PepeLP/P02/P02-quning';
grid_dir = '/mnt/40086D4C086D41D0/Inhibitory_excitatory_analysis/Meta/benware_grids/PLP/P02-quning';
dir_info = dir(fullfile(lfp_folder,'*lf.bin'));
dir_info2 = dir(fullfile(lfp_folder,'*lf.meta'));
lfpFilename = fullfile(lfp_folder,dir_info.name);
lfpMetaname = fullfile(lfp_folder,dir_info2.name);
meta_info = readSpikeGLXmeta(lfpMetaname);
%% Params
%Channel range for BF calculation
depth_start_um = 3500; %The topmost depth in um
depth_end_um = 300; %The bottommost depth in um
cortex_stop_um = 2000; %The extent of cortex in um from the beginning of L1 in um
l12_border_um = 3480; %The boundary between L1 and L2 in um
l1_start_um = l12_border_um + 150; %The beginning of L1 as defined from L1/L2 border
l2_end_um = l12_border_um - 400; %The end of L2/3 analysis window in um
%CSD analysis
spacing_ch_m = 40*10^-6; %The spacing between the channels in meters
uv_sample = meta_info.uV_per_bit_lfp; %The scaling factor of how many uv per unit of the raw signal 
electrode_radius_m = 12e-6; %The effective electrode radius in meters
%Ploting
bf_bw = 1; %How many freqeincies around BF to include in the mean LFP generation
notch = 1; %Whether to include a 50Hz notch filter of the data
bpass = 1; %Whether to bandpass the LFP
plot_fft_raw = 0; %Whether to plot the FFT spectrum of the raw data 
plot_fft_filt = 0; %Whether to plot the FFT spectrum of the filtered data
% PSTH
psth_end_ms = 150; %The end of the psth window in ms
psth_start_ms = -50; %The beginning of the psth window
sum_window_start_ms = 0; %The beginning of the summation window for generating FRAs in ms
sum_window_end_ms = 100; %The end of the summation window for generating FRAs in ms
stim_end_ms = 110;
%Sampling
lfp_fs = meta_info.imSampRate;
switch meta_info.imProbeOpt
    case {1,2,3}
        ref_chs = [37 76 113 152 189 228 265 304 341 380]; %The reference channels to be excluded from the recording
    case 4
        ref_chs = [37 76 113 152 189 228 265]; 
end
nChansInFile = meta_info.nChans;
%% Load the raw LFP, synch_ch

fprintf('== Loading the LFP ==\n');tic;
[lfp] = getLFPdata(lfpFilename,nChansInFile);
synch_ch = lfp(nChansInFile,:); %Save the synch_ch
lfp(nChansInFile,:) = []; %Delete the synch_ch from the data
lfp(ref_chs,:) = []; %Delete the ref channels 
n_func_chan = size(lfp,1);
fprintf('== Done! Loading took %.0fs ==\n',toc);
%Remove the mean from each channel to center at 0
lfp = lfp - mean(lfp,2);
N_time_points = length(lfp);
%% Filter the LFP with a notch filter to get rid of the 50Hz noise
% Design a filter with a Q-factor of Q=35 to remove a 50 Hz tone
if notch
    fprintf('== Notch filtering the LFP ==\n');tic;
    f0 = 50;
    Q = 35;
    Wo = f0/(lfp_fs/2);
    BW = Wo/Q;
    [b,a] = iirnotch(Wo,BW);
    lfp_filt = filtfilt(b,a,lfp');
    lfp_filt = lfp_filt';
    fprintf('== Done! Filtering took %.0fs ==\n',toc);
end
%% Bandpass filter the data to get rid of very slow oscilations and high freqeuncies which don't carry much information (Not sure?!)
if bpass
    fprintf('== Bandpass filtering the LFP ==\n');tic;
    fcutlow=2;   %low cut frequency in Hz
    fcuthigh=300;   %high cut frequency in Hz
    order = 2;
    [b,a]=butter(order,[fcutlow fcuthigh]/(lfp_fs/2));
    if ~exist('lfp_filt','var')
        lfp_filt = filtfilt(b,a,lfp');
    else
        lfp_filt = filtfilt(b,a,lfp_filt');
    end
    lfp_filt = lfp_filt';
    fprintf('== Done! Filtering took %.0fs ==\n',toc);
end
%% Generate the FFT of all raw channels
if plot_fft_raw
    fc = 150; %The center freqeuncy of the FFT
    fw = 150; %The bandwifth
    [P1,f] = fft_run(lfp',lfp_fs,[],fc,fw);
    f = f*1000; %Convert back to Hz
    P1 = P1';
    max_P1 = max(P1,[],2);
    P1 = 20*log10(P1./max_P1);
end
%% Generate the FFT of the filtered data
if plot_fft_filt
    fc = 150; %The center freqeuncy of the FFT
    fw = 150; %The bandwifth
    [P1_filt,f] = fft_run(lfp_filt',lfp_fs,[],fc,fw);
    f = f*1000; %Convert back to Hz
    P1_filt = P1_filt';
    max_P1 = max(P1_filt,[],2);
    P1_filt = 20*log10(P1_filt./max_P1);
end
%% Set things for the PSTH
t_bin_ms = (1/lfp_fs)*1000; %Time bin in ms
t_ms = [psth_start_ms:t_bin_ms:psth_end_ms]; %The edges of the histogram in ms
sum_window_start_bin = round((sum_window_start_ms - psth_start_ms)/t_bin_ms); %Specifies which bins to extract from the psth for generating the fra
sum_window_end_bin = round((sum_window_end_ms - psth_start_ms)/t_bin_ms);
t_samples = round((t_ms./1000)*lfp_fs); %The edges of the histogram in samples
num_time_points = length(t_samples); %Find the number of time points in the PSTH

 %% Load the benware info
dir_info = dir([grid_dir '/*Info.mat']); %Get the names of all files with this extension
grid_filename = fullfile(grid_dir, dir_info.name); %Form the name of the file to be opened
grid_load = load(grid_filename);
grid = grid_load.grid;
num_reps = grid.repeatsPerCondition; %Find the number of reps
min_trig_length_ms = grid.stimGrid(1,2); %The minimum length of one trigger in ms
min_trig_length_s = (min_trig_length_ms/1000); %The minimum length of one trigger in s
% min_inter_trig_length_s = grid.postStimSilence(1); %The time between two sweeps
min_inter_trig_length_s = 0.1; %Hack alex
dB_lvls = unique(grid.stimGrid(:,3));
dB_lvls = sort(dB_lvls,'descend');
num_dB_lvls = numel(dB_lvls); %Find the number of different dB levels used
freqs = unique(grid.stimGrid(:,1)); %Find the specific frequencies used
num_freqs = numel(freqs); %Find the number of different frequencies used
num_stim = grid.nStimConditions; %Total number of stimuli

%% Find the triggers
[start_time_ms] = get_triggers(synch_ch,min_trig_length_s,min_inter_trig_length_s,lfp_fs);
start_time_samples = round((start_time_ms./1000)*lfp_fs); %Convert to samples
num_triggers = numel(start_time_ms);
%% Make LFP PSTHs based on all stimuli for every channel
ch_start = round(depth_start_um/10); %Find the starting channel for the psth analysis. This is the one closer to the surface
ch_end = round(depth_end_um/10); %Find the starting channel for the psth analysis. This is the one closer to the tip
chans = [ch_end:ch_start];
num_chan_anal = (ch_start-ch_end)+1;
fprintf('== Extracting BF ==\n');tic;
for ch = 1:num_chan_anal
    current_chan = chans(ch);
%     X_dbft = zeros(num_stim,num_reps,num_time_points);
    X_dbft = [];
    for stim = 1:num_stim
        ix_rep = find(grid.randomisedGridSetIdx(1:num_triggers,1)==stim);
        ix_rep_samples = start_time_samples(ix_rep);
        num_actual_reps = length(ix_rep_samples); %Find out the number of repeats that were played
        for rep = 1:num_actual_reps
            X_dbft(stim,rep,:) = lfp_filt(current_chan,ix_rep_samples(rep) + t_samples);
        end
    end
    X_dbft(:,:,end) = []; %Delete the last bin which is weird
    X_dbft = mean(X_dbft,2);
    X_dbft = reshape(X_dbft,num_dB_lvls,num_freqs,numel(t_ms)-1);
    X_dbf = rms(X_dbft(:,:,sum_window_start_bin:sum_window_end_bin),3);
%     imagesc(flipud(X_dbf));
%     pause(0.3);
    X_f = mean(X_dbf); %Take the mean across levels
    [~,f_max(ch)] = max(X_f); %Find the ix of the max frequency
end
fprintf('== Done! BF extraction took %.0fs ==\n',toc);
ix_fmax = mode(f_max); %Find the most common BF
ix_fmax_l = ix_fmax - bf_bw;
ix_fmax_h = ix_fmax + bf_bw;
ix_keep_l = num_dB_lvls*(ix_fmax_l-1) + 1;
ix_keep_h = num_dB_lvls*(ix_fmax_h-1) + 1;
ix_keep = [ix_keep_l:(ix_keep_h+(num_dB_lvls-1))]; % Define the indices of the presentations to keep including only BF 
%% Different
[log_vec,~] = ismember(grid.randomisedGridSetIdx(1:num_triggers),ix_keep);
start_time_samples_bf = start_time_samples(log_vec);
num_bf_triggers = length(start_time_samples_bf);
LFP_chtr = zeros(n_func_chan,num_time_points,num_bf_triggers);
LFP_filt_chtr = zeros(n_func_chan,num_time_points,num_bf_triggers);
fprintf('== Generating mean LFP ==\n');
for tr = 1:num_bf_triggers
    LFP_chtr(:,:,tr) = lfp(:,start_time_samples_bf(tr) + t_samples); %
    LFP_filt_chtr(:,:,tr) = lfp_filt(:,start_time_samples_bf(tr) + t_samples); %
end
mean_lfp = mean(LFP_chtr,3)*uv_sample;
mean_lfp_filt= mean(LFP_filt_chtr,3)*uv_sample;
clear LFP_chtr;
clear LFP_filt_chtr;
%% Load probe geometry info
[fp1,fp2] = fileparts(lfpFilename);
[root_dir,file_name] = fileparts(fp1);
chan_map = load(fullfile(root_dir,['config_dir/chanMap.mat']));
ycoords = chan_map.ycoords;
depth_um = ycoords;
ycoords(end) = [];
depth_um(end) = []; %Remove the last channel
depth_um(ref_chs) = []; %Remove ref channels 

xcoords = chan_map.xcoords;
xcoords(end) = [];
xcoords(ref_chs) = []; 
% 22/07/19

% Define the indices of the columns
col1_ix = find(xcoords==xcoords(1));
depth_um_col1 = depth_um(col1_ix);
depth_um_col1 = flipud(depth_um_col1);

col2_ix = find(xcoords==xcoords(2));
depth_um_col2 = depth_um(col2_ix);
depth_um_col2 = flipud(depth_um_col2);

col3_ix = find(xcoords==xcoords(3));
depth_um_col3 = depth_um(col3_ix);
depth_um_col3 = flipud(depth_um_col3);

col4_ix = find(xcoords==xcoords(4));
depth_um_col4 = depth_um(col4_ix);
depth_um_col4 = flipud(depth_um_col4);

xcoords = flipud(xcoords);
ycoords = flipud(ycoords);
depth_um = flipud(depth_um);
%% Run the CSD analysis
%Convert to the data to volts
mean_lfpn = mean_lfp./1e+6; 
mean_lfp_filtn = mean_lfp_filt./1e+6;
%Run the inverse method
[CSDout]  = CSD(fliplr(mean_lfpn'),lfp_fs,spacing_ch_m,'inverse',electrode_radius_m);
CSDout = CSDout';

[CSD_filtout]  = CSD(fliplr(mean_lfp_filtn'),lfp_fs,spacing_ch_m,'inverse',electrode_radius_m);
CSD_filtout = CSD_filtout';

[CSD_filtout_col1]  = CSD(fliplr(mean_lfp_filtn(col1_ix,:)'),lfp_fs,spacing_ch_m,'inverse',electrode_radius_m);
CSD_filtout_col1 = CSD_filtout_col1';

[CSD_filtout_col2]  = CSD(fliplr(mean_lfp_filtn(col2_ix,:)'),lfp_fs,spacing_ch_m,'inverse',electrode_radius_m);
CSD_filtout_col2 = CSD_filtout_col2';

[CSD_filtout_col3]  = CSD(fliplr(mean_lfp_filtn(col3_ix,:)'),lfp_fs,spacing_ch_m,'inverse',electrode_radius_m);
CSD_filtout_col3 = CSD_filtout_col3';

[CSD_filtout_col4]  = CSD(fliplr(mean_lfp_filtn(col4_ix,:)'),lfp_fs,spacing_ch_m,'inverse',electrode_radius_m);
CSD_filtout_col4 = CSD_filtout_col4';

%% Plot the mean tone-evoked raw LFP
start_plot_ms = -40; %Where to start the plotting for the time axis
[~,start_plot_ix] = min(abs(start_plot_ms - t_ms)); %Where to start
line_width = 3;
font_size_axes = 50;
font_size = 35;
ticks_font_size = 40;
skip_y = 10;
y_step_um = 200; %The step in y in um
skip_x = 50;
depth_stop_ix = find(depth_um == (depth_um(1)-cortex_stop_um));
depth_stop_ix = depth_stop_ix(1);
depth_trans_um = abs(depth_um - depth_um(1));
[~,y_ix,~] = intersect(depth_trans_um,[0:y_step_um:cortex_stop_um]);

font = 'Liberation Sans';
figure('DefaultTextFontName', font, 'DefaultAxesFontName',font,'DefaultAxesFontSize',font_size);
plot_lfp = flipud(mean_lfp);
plot_lfp = plot_lfp(1:depth_stop_ix,:);
max_val = max(abs(plot_lfp(:)));
imagesc(plot_lfp);
hold on;
stim_on = find(t_ms==0); %Define beginning of the stimulus
stim_off = find(t_ms==stim_end_ms); %Define beginning of the stimulus
[~,l12_border] = min(abs(depth_um - l12_border_um)); %Define beginning of the stimulus
[~,l2_end] = min(abs(depth_um - l2_end_um)); %Define beginning of the stimulus
% xline(stim_off,'--k','LineWidth',line_width);
% xline(stim_on,'--k','LineWidth',line_width);
% yline(l12_border,'--k','LineWidth',line_width);
% yline(l2_end,'--k','LineWidth',line_width);
colormap(flipud(jet));
caxis([-max_val max_val]);
colorbar;
hold off;
%Make labels
for jj = 1:numel(y_ix)
    ix = y_ix(jj);
    y_labels{jj} = num2str(depth_trans_um(ix),'%.0f');
end
yticks(y_ix);
yticklabels(y_labels);

xticks([start_plot_ix:skip_x:num_time_points]);
x_labels = string(t_ms(start_plot_ix:skip_x:end));
xticklabels(x_labels);
xlabel('Time relative to stimulus onset [ms]','FontSize',font_size_axes,'FontWeight','Bold');
ylabel('Depth [\mum]','FontSize',font_size_axes,'FontWeight','Bold');
xl = get(gca,'XLabel');
xAX = get(gca,'XAxis');
set(xAX,'FontSize', ticks_font_size)
set(xl, 'FontSize', font_size_axes);
yl = get(gca,'YLabel');
yAX = get(gca,'YAxis');
set(yAX,'FontSize', ticks_font_size)
set(yl, 'FontSize', font_size_axes);
% title(['Mean tone evoked LFP F ',file_name],'FontSize',font_size,'FontWeight','Bold');
cb1 = colorbar;
colorTitleHandle = get(cb1,'Title');
titleString = '[\muV]';
set(colorTitleHandle ,'String',titleString);
set(gcf,'color','w');
%% Plot the mean tone-evoked filtered LFP
font = 'Liberation Sans';
figure('DefaultTextFontName', font, 'DefaultAxesFontName',font,'DefaultAxesFontSize',font_size);
plot_lfp_filt = flipud(mean_lfp_filt);
plot_lfp_filt = plot_lfp_filt(1:depth_stop_ix,:);
max_val = max(abs(plot_lfp_filt(:)));
imagesc(plot_lfp_filt);
hold on;
% xline(stim_on,'--k','LineWidth',line_width);
% xline(stim_off,'--k','LineWidth',line_width);
% yline(l12_border,'--k','LineWidth',line_width);
% yline(l2_end,'--k','LineWidth',line_width);
colormap(flipud(jet));
caxis([-max_val max_val]);
colorbar;
hold off;

yticks(y_ix);
yticklabels(y_labels);
xticks([start_plot_ix:skip_x:num_time_points]);
x_labels = string(t_ms(start_plot_ix:skip_x:end));
xticklabels(x_labels);
xlabel('Time relative to stimulus onset [ms]','FontSize',font_size_axes,'FontWeight','Bold');
ylabel('Depth [\mum]','FontSize',font_size_axes,'FontWeight','Bold');
xl = get(gca,'XLabel');
xlFontSize = get(xl,'FontSize');
xAX = get(gca,'XAxis');
set(xAX,'FontSize', ticks_font_size)
set(xl, 'FontSize', font_size_axes);
yl = get(gca,'YLabel');
ylFontSize = get(yl,'FontSize');
yAX = get(gca,'YAxis');
set(yAX,'FontSize', ticks_font_size)
set(yl, 'FontSize', font_size_axes);
% title('Average LFP response at BF','FontSize',font_size,'FontWeight','Bold');
cb2 = colorbar;
colorTitleHandle = get(cb2,'Title');
titleString = '[\muV]';
set(colorTitleHandle ,'String',titleString);
set(gcf,'color','w');
%% Plot the CSD for the raw tone-evoked LFP

max_val = max(abs(CSDout(:)));
imagesc(CSDout);
hold on;
xline(stim_on,'--k','LineWidth',2);
xline(stim_off,'--k','LineWidth',2);
colormap(flipud(jet));
caxis([-max_val max_val]);
colorbar;
hold off;

yticks([1:skip_y:n_func_chan]);
yticklabels(y_labels);
xticks([1:skip_x:num_time_points]);
x_labels = string(t_ms(1:skip_x:end));
xticklabels(x_labels);

xlabel('Time relative to stimulus onset [ms]','FontSize',font_size,'FontWeight','Bold');
ylabel('Depth [\mum]','FontSize',font_size,'FontWeight','Bold');
ax = gca;
ax.FontSize = ticks_font_size; 
title(['Mean tone-evoked CSD F ',file_name, ' [\color{blue}Source \color{red}Sink\color{black}]' ],'FontSize',font_size,'FontWeight','Bold');
cb2 = colorbar;
colorTitleHandle = get(cb2,'Title');
titleString = 'Amplitude [\muA/mm^3]';
set(colorTitleHandle ,'String',titleString);
set(gcf,'color','w');
%% Plot the CSD for the filtered tone-evoked LFP

max_val = max(abs(CSD_filtout(:)));
imagesc(CSD_filtout);
hold on;
xline(stim_on,'--k','LineWidth',2);
xline(stim_off,'--k','LineWidth',2);
colormap(flipud(jet));
caxis([-max_val max_val]);
colorbar;
hold off;

yticks([1:skip_y:n_func_chan]);
yticklabels(y_labels);
xticks([1:skip_x:num_time_points]);
x_labels = string(t_ms(1:skip_x:end));
xticklabels(x_labels);

xlabel('Time relative to stimulus onset [ms]','FontSize',font_size,'FontWeight','Bold');
ylabel('Depth [\mum]','FontSize',font_size,'FontWeight','Bold');
ax = gca;
ax.FontSize = ticks_font_size; 
title(['Mean tone-evoked filt CSD F ',file_name,' [\color{blue}Source \color{red}Sink\color{black}]'],'FontSize',font_size,'FontWeight','Bold');
cb2 = colorbar;
colorTitleHandle = get(cb2,'Title');
titleString = 'Amplitude [\muA/mm^3]';
set(colorTitleHandle ,'String',titleString);
set(gcf,'color','w');
%% Plot the CSD for the filtered tone-evoked LFP for probe columns 2,4,1,3
%Column 2
font = 'Liberation Sans';
figure('DefaultTextFontName', font, 'DefaultAxesFontName',font,'DefaultAxesFontSize',font_size);
%Find the um and ix of the channel closest to the beginning of L1
[~,l1_start_col2_ix] = min(abs(depth_um_col2 - l1_start_um));
[~,l12_border_col2] = min(abs(depth_um_col2(l1_start_col2_ix:end) - l12_border_um)); %Define the border of L1/L2
[~,l2_end_col2] = min(abs(depth_um_col2(l1_start_col2_ix:end) - l2_end_um)); %Define end of L2/3
depth_stop_col2_ix = find(depth_um_col2 == (depth_um_col2(l1_start_col2_ix)-cortex_stop_um));
depth_stop_col2_ix = depth_stop_col2_ix(1);
depth_trans_col2_um = abs(depth_um_col2(l1_start_col2_ix:end) - depth_um_col2(l1_start_col2_ix));
[~,y_ix,~] = intersect(depth_trans_col2_um,[0:y_step_um:cortex_stop_um]);
plot_lfp_col2 = CSD_filtout_col2(l1_start_col2_ix:depth_stop_col2_ix,:);
max_val = max(abs(plot_lfp_col2(:)));
imagesc(plot_lfp_col2);
hold on;
% xline(stim_on,'--k','LineWidth',2);
% xline(stim_off,'--k','LineWidth',2);
% yline(l12_border_col2,'--k','LineWidth',line_width);
% yline(l2_end_col2,'--k','LineWidth',line_width);
colormap(flipud(jet));
caxis([-max_val max_val]);
colorbar;
hold off;
%Make labels
for jj = 1:numel(y_ix)
    ix = y_ix(jj);
    y_labels{jj} = num2str(depth_trans_col2_um(ix),'%.0f');
end
yticks(y_ix);
yticklabels(y_labels);
xticks([start_plot_ix:skip_x:num_time_points]);
x_labels = string(t_ms(start_plot_ix:skip_x:end));
xticklabels(x_labels);

xlabel('Time relative to stimulus onset [ms]','FontSize',font_size_axes,'FontWeight','Bold');
ylabel('Depth [\mum]','FontSize',font_size_axes,'FontWeight','Bold');
xl = get(gca,'XLabel');
xAX = get(gca,'XAxis');
set(xAX,'FontSize', ticks_font_size)
set(xl, 'FontSize', font_size_axes);
yl = get(gca,'YLabel');
yAX = get(gca,'YAxis');
set(yAX,'FontSize', ticks_font_size)
set(yl, 'FontSize', font_size_axes);
% title(['Average CSD response at BF ',' [\color{blue}Source \color{red}Sink\color{black}]'],'FontSize',font_size,'FontWeight','Bold');
cb2 = colorbar;
colorTitleHandle = get(cb2,'Title');
titleString = '[\muA/mm^3]';
set(colorTitleHandle ,'String',titleString);
set(gcf,'color','w');
%% Column 4
font = 'Liberation Sans';
figure('DefaultTextFontName', font, 'DefaultAxesFontName',font,'DefaultAxesFontSize',font_size);
[~,l12_border_col4] = min(abs(depth_um_col4 - l12_border_um)); %Define beginning of the stimulus
[~,l2_end_col4] = min(abs(depth_um_col4 - l2_end_um)); %Define beginning of the stimulus
depth_stop_col4_ix = find(depth_um_col4 == (depth_um_col4(1)-cortex_stop_um));
depth_stop_col4_ix = depth_stop_col4_ix(1);
depth_trans_col4_um = abs(depth_um_col4 - depth_um_col4(1));
[~,y_ix,~] = intersect(depth_trans_col4_um,[0:y_step_um:cortex_stop_um]);
plot_lfp_col4 = CSD_filtout_col4(1:depth_stop_col4_ix,:);
max_val = max(abs(plot_lfp_col4(:)));
imagesc(plot_lfp_col4);
hold on;
stim_on = find(t_ms==0); %Define beginning of the stimulus
% xline(stim_on,'--k','LineWidth',2);
yline(l12_border_col4,'--k','LineWidth',line_width);
yline(l2_end_col4,'--k','LineWidth',line_width);
colormap(flipud(jet));
caxis([-max_val max_val]);
colorbar;
hold off;
%Make labels
for jj = 1:numel(y_ix)
    ix = y_ix(jj);
    y_labels{jj} = num2str(depth_trans_col4_um(ix),'%.0f');
end
yticks(y_ix);
yticklabels(y_labels);
xticks([start_plot_ix:skip_x:num_time_points]);
x_labels = string(t_ms(start_plot_ix:skip_x:end));
xticklabels(x_labels);

xlabel('Time relative to stimulus onset [ms]','FontSize',font_size_axes,'FontWeight','Bold');
ylabel('Depth [\mum]','FontSize',font_size_axes,'FontWeight','Bold');
xl = get(gca,'XLabel');
xAX = get(gca,'XAxis');
set(xAX,'FontSize', ticks_font_size)
set(xl, 'FontSize', font_size_axes);
yl = get(gca,'YLabel');
yAX = get(gca,'YAxis');
set(yAX,'FontSize', ticks_font_size)
set(yl, 'FontSize', font_size_axes);
% title(['Average CSD response at BF ',' [\color{blue}Source \color{red}Sink\color{black}]'],'FontSize',font_size,'FontWeight','Bold');
cb2 = colorbar;
colorTitleHandle = get(cb2,'Title');
titleString = '[\muA/mm^3]';
set(colorTitleHandle ,'String',titleString);
set(gcf,'color','w');
%Column 1
font = 'Liberation Sans';
figure('DefaultTextFontName', font, 'DefaultAxesFontName',font,'DefaultAxesFontSize',font_size);
[~,l12_border_col1] = min(abs(depth_um_col1 - l12_border_um)); %Define beginning of the stimulus
[~,l2_end_col1] = min(abs(depth_um_col1 - l2_end_um)); %Define beginning of the stimulus
depth_stop_col1_ix = find(depth_um_col1 == (depth_um_col1(1)-cortex_stop_um));
depth_stop_col1_ix = depth_stop_col1_ix(1);
depth_trans_col1_um = abs(depth_um_col1 - depth_um_col1(1));
[~,y_ix,~] = intersect(depth_trans_col1_um,[0:y_step_um:cortex_stop_um]);
plot_lfp_col1 = CSD_filtout_col1(1:depth_stop_col1_ix,:);
max_val = max(abs(plot_lfp_col1(:)));
imagesc(plot_lfp_col1);
hold on;
% xline(stim_on,'--k','LineWidth',2);
% xline(stim_off,'--k','LineWidth',2);
yline(l12_border_col1,'--k','LineWidth',line_width);
yline(l2_end_col1,'--k','LineWidth',line_width);
colormap(flipud(jet));
caxis([-max_val max_val]);
colorbar;
hold off;
%Make labels
for jj = 1:numel(y_ix)
    ix = y_ix(jj);
    y_labels{jj} = num2str(depth_trans_col1_um(ix),'%.0f');
end
yticks(y_ix);
yticklabels(y_labels);
xticks([start_plot_ix:skip_x:num_time_points]);
x_labels = string(t_ms(start_plot_ix:skip_x:end));
xticklabels(x_labels);

xlabel('Time relative to stimulus onset [ms]','FontSize',font_size_axes,'FontWeight','Bold');
ylabel('Depth [\mum]','FontSize',font_size_axes,'FontWeight','Bold');
xl = get(gca,'XLabel');
xAX = get(gca,'XAxis');
set(xAX,'FontSize', ticks_font_size)
set(xl, 'FontSize', font_size_axes);
yl = get(gca,'YLabel');
yAX = get(gca,'YAxis');
set(yAX,'FontSize', ticks_font_size)
set(yl, 'FontSize', font_size_axes);
% title(['Average CSD response at BF ',' [\color{blue}Source \color{red}Sink\color{black}]'],'FontSize',font_size,'FontWeight','Bold');
cb2 = colorbar;
colorTitleHandle = get(cb2,'Title');
titleString = '[\muA/mm^3]';
set(colorTitleHandle ,'String',titleString);
set(gcf,'color','w');
%Column 2
font = 'Liberation Sans';
figure('DefaultTextFontName', font, 'DefaultAxesFontName',font,'DefaultAxesFontSize',font_size);
[~,l12_border_col2] = min(abs(depth_um_col2 - l12_border_um)); %Define beginning of the stimulus
[~,l2_end_col2] = min(abs(depth_um_col2 - l2_end_um)); %Define beginning of the stimulus
depth_stop_col2_ix = find(depth_um_col2 == (depth_um_col2(1)-cortex_stop_um));
depth_stop_col2_ix = depth_stop_col2_ix(1);
depth_trans_col2_um = abs(depth_um_col2 - depth_um_col2(1));
[~,y_ix,~] = intersect(depth_trans_col2_um,[0:y_step_um:cortex_stop_um]);
plot_lfp_col2 = CSD_filtout_col1(1:depth_stop_col2_ix,:);
max_val = max(abs(plot_lfp_col2(:)));
imagesc(plot_lfp_col1);
hold on;
% xline(stim_on,'--k','LineWidth',2);
% xline(stim_off,'--k','LineWidth',2);
yline(l12_border_col2,'--k','LineWidth',line_width);
yline(l2_end_col2,'--k','LineWidth',line_width);
colormap(flipud(jet));
caxis([-max_val max_val]);
colorbar;
hold off;
%Make labels
for jj = 1:numel(y_ix)
    ix = y_ix(jj);
    y_labels{jj} = num2str(depth_trans_col2_um(ix),'%.0f');
end
yticks(y_ix);
yticklabels(y_labels);
xticks([start_plot_ix:skip_x:num_time_points]);
x_labels = string(t_ms(start_plot_ix:skip_x:end));
xticklabels(x_labels);

xlabel('Time relative to stimulus onset [ms]','FontSize',font_size_axes,'FontWeight','Bold');
ylabel('Depth [\mum]','FontSize',font_size_axes,'FontWeight','Bold');
xl = get(gca,'XLabel');
xAX = get(gca,'XAxis');
set(xAX,'FontSize', ticks_font_size)
set(xl, 'FontSize', font_size_axes);
yl = get(gca,'YLabel');
yAX = get(gca,'YAxis');
set(yAX,'FontSize', ticks_font_size)
set(yl, 'FontSize', font_size_axes); 
% title(['Average CSD response at BF ',' [\color{blue}Source \color{red}Sink\color{black}]'],'FontSize',font_size,'FontWeight','Bold');
cb2 = colorbar;
colorTitleHandle = get(cb2,'Title');
titleString = '[\muA/mm^3]';
set(colorTitleHandle ,'String',titleString);
set(gcf,'color','w');
%% Do a low-tech CSD version as a sanity check for column 3
% CSD_check_col3 = conv2(mean_lfp_filt(col3_ix,:),[1 -2 1],'same');
% %Column 3
% h(9) = figure;
% max_val = max(abs(CSD_check_col3(:)));
% imagesc(flipud(CSD_check_col3));
% hold on;
% stim_on = find(t_ms==0); %Define beginning of the stimulus
% xline(stim_on,'--k','LineWidth',2);
% colormap(flipud(jet));
% caxis([-0.3 0.3]);
% colorbar;
% hold off;
% %Make labels
% yticks([1:skip_y:length(depth_um_col3)]);
% yticklabels(y_labels);
% xticks([1:skip_x:num_time_points]);
% x_labels = string(t_ms(1:skip_x:end));
% xticklabels(x_labels);
% 
% xlabel('Time relative to stimulus onset [ms]','FontSize',font_size,'FontWeight','Bold');
% ylabel('Depth [\mum]','FontSize',font_size,'FontWeight','Bold');
% ax = gca;
% ax.FontSize = ticks_font_size; 
% title(['Average CSD response at BF ',file_name,' [\color{blue}Source \color{red}Sink\color{black}]'],'FontSize',font_size,'FontWeight','Bold');
% cb2 = colorbar;
% colorTitleHandle = get(cb2,'Title');
% titleString = 'Amplitude [\muA/mm^3]';
% set(colorTitleHandle ,'String',titleString);
%% Do a low-tech CSD version as a sanity check for all channels
% CSD_check_all = conv2(mean_lfp_filt,[1 -2 1],'same');
% %Column 3
% h(10) = figure;
% imagesc(flipud(CSD_check_all));
% hold on;
% stim_on = find(t_ms==0); %Define beginning of the stimulus
% xline(stim_on,'--k','LineWidth',2);
% colormap(flipud(jet));
% caxis([-0.3 0.3]);
% colorbar;
% hold off;
% %Make labels
% %Make labels
% count = 0;
% for jj = 1:skip_y:n_func_chan
%     count = count+1;
%     y_labels{count} = num2str(depth_um(jj),'%.0f');
% end
% yticks([1:skip_y:n_func_chan]);
% yticklabels(y_labels);
% 
% xticks([1:skip_x:num_time_points]);
% x_labels = string(t_ms(1:skip_x:end));
% xticklabels(x_labels);
% 
% xlabel('Time relative to stimulus onset [ms]','FontSize',font_size,'FontWeight','Bold');
% ylabel('Depth [\mum]','FontSize',font_size,'FontWeight','Bold');
% ax = gca;
% ax.FontSize = ticks_font_size; 
% title(['Mean tone-evoked filt All Chs  Manual CSD F ',file_name,' [\color{blue}Source \color{red}Sink\color{black}]'],'FontSize',font_size,'FontWeight','Bold');
% cb2 = colorbar;
% colorTitleHandle = get(cb2,'Title');
% titleString = 'Amplitude [\muA/mm^3]';
% set(colorTitleHandle ,'String',titleString);
%% Save params
params.depth_start = depth_start_um; %The topmost depth in um
params.depth_end = depth_end_um; %The bottommost depth in um
params.bf_bw = bf_bw; %How many freqeincies around BF to include in the mean LFP generation
params.notch = notch; %Whether to include a 50Hz notch filter of the data
params.notch_f0 = f0;
params.notch_Q = Q;
params.bpass = bpass; %Whether to bandpass the LFP
params.bpass_fcutlow=fcutlow;   %low cut frequency in Hz
params.bpass_fcuthigh=fcuthigh;   %high cut frequency in Hz
params.bpass_order = order;
%CSD analysis
params.spacing_ch_m = spacing_ch_m; %The spacing between the channels in meters
params.uv_sample = uv_sample; %The scaling factor of how many uv per unit of the raw signal 
params.electrode_radius_m = electrode_radius_m; %The effective electrode radius in meters
%PSTH
params.psth_end_ms = psth_end_ms; %The end of the psth window in ms
params.psth_start_ms = psth_start_ms; %The beginning of the psth window
params.sum_window_start_ms = sum_window_start_ms; %The beginning of the summation window for generating FRAs in ms
params.sum_window_end_ms = sum_window_end_ms; %The end of the summation window for generating FRAs in ms
%Sampling
params.lfp_fs = lfp_fs;
params.nChansInFile = nChansInFile;
params.ref_chs = ref_chs;

save_dir = fullfile(fp1,'/LFP_analysis');
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end

savefig(h,fullfile(save_dir,'LFP_analysis_figures.fig'));
save(fullfile(save_dir,'params.mat'),'params');
%% Plot the FFT of the raw data
if plot_fft_raw
    f_res = lfp_fs/N_time_points;
    f_res_desired = 10; %In Hz
    skip_x = round(f_res_desired/f_res);
    font_size = 25;
    ticks_font_size = 16;
    figure;
    num_freq_fft = length(f);
    imagesc(flipud(P1));caxis([-100 0]);
    xticks([1:skip_x:num_freq_fft]);
    f = round(f);
    x_labels = string(f(1:skip_x:end));
    xticklabels(x_labels);
    cb = colorbar;
    colorTitleHandle = get(cb,'Title');
    titleString = 'Power [dB]';
    set(colorTitleHandle ,'String',titleString);
    yticks([1:skip_y:n_func_chan]);
    yticklabels(y_labels);
    xlabel('Frequency [Hz]','FontSize',font_size,'FontWeight','Bold');
    ylabel('Depth [\mum]','FontSize',font_size,'FontWeight','Bold');
    ax = gca;
    ax.FontSize = ticks_font_size;
    title(['Unfiltered LFP spectrum for ferret ',file_name],'FontSize',font_size,'FontWeight','Bold');
end
%% Plot the FFT of the filtered data
if plot_fft_filt
    figure;
    imagesc(flipud(P1_filt));caxis([-100 0]);
    xticks([1:skip_x:num_freq_fft]);
    xticklabels(x_labels);
    cb = colorbar;
    colorTitleHandle = get(cb,'Title');
    titleString = 'Power [dB]';
    set(colorTitleHandle ,'String',titleString);
    yticks([1:skip_y:n_func_chan]);
    yticklabels(y_labels);
    xlabel('Frequency [Hz]','FontSize',font_size,'FontWeight','Bold');
    ylabel('Depth [\mum]','FontSize',font_size,'FontWeight','Bold');
    ax = gca;
    ax.FontSize = ticks_font_size;
    title(['Filtered LFP spectrum for ferret ',file_name],'FontSize',font_size,'FontWeight','Bold');
end