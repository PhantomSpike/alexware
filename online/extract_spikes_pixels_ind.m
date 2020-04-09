function [spike_times_ms,synch_ch,fs] = extract_spikes_pixels_ind(data_root,channel_no,start_time_s,end_time_s,method)

%% Parameters
%choi = chanel of interest
spike_th = 3; %Spike threshold expressed as number of standard deviations
ch_radius = 15; %Number of channels above and below the choi used for CRA
fs = 30000; %Sampling frequncy
no_channels = 385; %Total number of channels
bytes_sample = 2; %Number of bytes per one sample

first_ch = channel_no - ch_radius;
last_ch = channel_no + ch_radius;
total_channels = 2*ch_radius + 1;
total_time_s = end_time_s - start_time_s;

%% Load the synch channel
total_samples_synch = total_time_s*fs;

dir_info = dir([data_root '/*ap.bin']); %Get the names of all files with this extension
bin_filename = fullfile(dir_info.folder, dir_info.name); %Form the name of the file to be opened

offset_synch = (start_time_s*fs*no_channels + 384)*bytes_sample; %This should move the pointer to the 385 channel i.e. the synch channel

tic;
fprintf('== Loading Channel #%.0f ==\n',channel_no);
fid = fopen(bin_filename, 'r'); %Open an fid for the file of interest. This is something like a pointer for this file

fseek(fid, offset_synch, 'bof'); %Move the place to which fseek is pointing

precision = '*int16';
skip = 384*bytes_sample;
synch_ch = fread(fid,[total_samples_synch,1],precision,skip);

fclose(fid); %Close the fid so it doesn't interfere with future calls of fopen

%% Load the choi + ref channels
total_samples_chs = total_time_s*fs;

offset_choi = (start_time_s*fs*no_channels + first_ch - 1)*bytes_sample; %This should move the pointer to the first channel included in the averaging radius around the choi

fid = fopen(bin_filename, 'r'); %Open an fid for the file of interest. This is something like a pointer for this file

fseek(fid, offset_choi, 'bof'); %Move the place to which fseek is pointing

precision = [num2str(total_channels),'*int16=>int16'];
skip = (no_channels - total_channels)*bytes_sample;

channels_matrix = fread(fid,[total_channels,total_samples_chs],precision,skip);

fclose(fid); %Close the fid so it doesn't interfere with future calls of fopen
fprintf('== Done! Loading took %0.f sec ==\n',toc);
%% Common Reference Averaging
tic;
fprintf('== Performing CRA ==\n');
%Subtract the median values across all times points for each channel to
%bring mean to 0
med_vals_time = median(channels_matrix,2);
channels_matrix = channels_matrix - med_vals_time;

%Do the CRA
med_val_chs = median(channels_matrix,1);
choi = channels_matrix(ch_radius+1,:); 
choi = choi - med_val_chs;

fprintf('== Done! CRA took %0.f sec ==\n',toc);
%% Filtering 
fprintf('== Applying eliptic filter to the data ==\n');tic;
Wp = [300 3000]; %High and low-pass frequencies
n = 7; %Filter order. Choose a low one so it's not comp expensive
[spikeFilter.B,spikeFilter.A] = ellip(n, 0.01, 40, Wp/(fs/2));

choi = double(choi); %Convert choi to type double because filtering doesn't work otherwise
choi = filtfilt(spikeFilter.B, spikeFilter.A, choi); %Filter the signal
fprintf('== Done! Filtering took %.1f sec ==\n',toc);

%% Extract the spike times usning Quentin's Peak Detector Code or Benware method

fprintf('== Extracting spikes from the data ==\n');tic;

if ~exist('method','var')
    method = 'b';
end

switch method
    case 'q'
        [bin_peak_max,val_peak_max,peak_val,peak_abs] = peak_detector_general(choi,'mean',spike_th);
        spike_times_ms = (bin_peak_max/fs)*1000;
    case 'b'
        mn = mean(choi); %mean across all time points for each channel
        sd = std(choi);
          
        choi_processed = ((choi-mn)/sd) - spike_th;
        
        % find threshold crossings
        choi_sign = sign(choi_processed);
        spike_ix = find(diff(choi_sign)<0);
        
        % append new spike times to spikeTimes
        spike_times_ms = (spike_ix/fs)*1000;
end
        
fprintf('== Done! Extraction took %.1f sec ==\n',toc);