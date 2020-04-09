function [ch_data] = plot_channel_tmp(ch_no,bin_fullname,start_time_s,end_time_s)
fs = 30000; %Sampling frequncy
%Select the total number of channels based on the probe option 
no_channels = 271;

bytes_sample = 2; %Number of bytes per one sample
total_time_s = end_time_s - start_time_s;
%% Load the synch channel
total_samples_synch = total_time_s*fs;

offset_synch = (start_time_s*fs*no_channels + ch_no)*bytes_sample; %This should move the pointer to the 385 channel i.e. the synch channel

fid = fopen(bin_fullname, 'r'); %Open an fid for the file of interest. This is something like a pointer for this file

fseek(fid, offset_synch, 'bof'); %Move the place to which fseek is pointing

precision = '*int16';
skip = (no_channels-1)*bytes_sample;
ch_data = fread(fid,[total_samples_synch,1],precision,skip);
ch_data = ch_data - mean(ch_data);
rms_ch = rms(ch_data);

fclose(fid); %Close the fid so it doesn't interfere with future calls of fopen

%% Plot the channel
font_axis = 12;
t = [start_time_s:1/fs:end_time_s-1/fs];
figure;
plot(t,ch_data);
dim = [.8 .6 .3 .3];
anot = annotation('textbox',dim,'String',sprintf('RMS = %.1fuV',rms_ch),'FitBoxToText','on');
anot.FontSize = 12;
anot.FontWeight = 'bold';
set(gca,'FontName','Arial','FontSize',font_axis,'FontWeight','Bold');
xlabel('Time [s]');
ylabel('Amplitude [uV]');
title(['Channel # ',num2str(ch_no)]);