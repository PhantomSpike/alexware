function [spike_times_ms] = extract_spikes_pixels(bin_fullname,channel_no,start_time_s,end_time_s,method,probe)

%% Parameters
%choi = chanel of interest
spike_th = 3; %Spike threshold expressed as number of standard deviations
ch_radius = 15; %Number of channels above and below the choi used for CRA
fs = 30000; %Sampling frequncy
%Select the total number of channels based on the probe option
switch probe
    case 1
        no_channels = 385; %Total number of channels
    case 2
        no_channels = 385; %Total number of channels
    case 3
        no_channels = 385; %Total number of channels
    case 4
        no_channels = 277; %Total number of channels
end
bytes_sample = 2; %Number of bytes per one sample

first_ch = channel_no - ch_radius;
last_ch = channel_no + ch_radius;
total_channels = 2*ch_radius + 1;
total_time_s = end_time_s - start_time_s;

%% Load the choi + ref channels
total_samples_chs = total_time_s*fs;

offset_choi = (start_time_s*fs*no_channels + first_ch - 1)*bytes_sample; %This should move the pointer to the first channel included in the averaging radius around the choi

fid = fopen(bin_fullname, 'r'); %Open an fid for the file of interest. This is something like a pointer for this file

fseek(fid, offset_choi, 'bof'); %Move the place to which fseek is pointing

precision = [num2str(total_channels),'*int16=>int16'];
skip = (no_channels - total_channels)*bytes_sample;

channels_matrix = fread(fid,[total_channels,total_samples_chs],precision,skip);

fclose(fid); %Close the fid so it doesn't interfere with future calls of fopen
%% Common Reference Averaging
%Subtract the median values across all times points for each channel to
%bring mean to 0
med_vals_time = median(channels_matrix,2);
channels_matrix = channels_matrix - med_vals_time;

%Do the CRA
med_val_chs = median(channels_matrix,1);
choi = channels_matrix(ch_radius+1,:); 
choi = choi - med_val_chs;
%% Filtering 
Wp = [300 3000]; %High and low-pass frequencies
n = 7; %Filter order. Choose a low one so it's not comp expensive
[spikeFilter.B,spikeFilter.A] = ellip(n, 0.01, 40, Wp/(fs/2));

choi = double(choi); %Convert choi to type double because filtering doesn't work otherwise
choi = filtfilt(spikeFilter.B, spikeFilter.A, choi); %Filter the signal

%% Extract the spike times usning Quentin's Peak Detector Code or Benware method

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
