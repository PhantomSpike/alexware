function [synch_ch] = get_synch_online(bin_fullname,start_time_s,end_time_s,probe)
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
total_time_s = end_time_s - start_time_s;
%% Load the synch channel
total_samples_synch = total_time_s*fs;

offset_synch = (start_time_s*fs*no_channels + (no_channels-1))*bytes_sample; %This should move the pointer to the 385 channel i.e. the synch channel

fid = fopen(bin_fullname, 'r'); %Open an fid for the file of interest. This is something like a pointer for this file

fseek(fid, offset_synch, 'bof'); %Move the place to which fseek is pointing

precision = '*int16';
skip = (no_channels-1)*bytes_sample;
synch_ch = fread(fid,[total_samples_synch,1],precision,skip);

fclose(fid); %Close the fid so it doesn't interfere with future calls of fopen