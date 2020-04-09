function [cropped_data,total_time_samples] = load_spikedata(data_root,start_time,end_time,probe)
%cropped_data = load_spikedata(data_root)
%This function loads a selected chunk of a SpikeGLX-recorded spike data
%which is in binary format.
%>>>Input>> 
%data_root = The location of the file to be loaded
%start_time = The start time of the chunk in sec relative to the
%beginning of the recording
%end_time = The end time of the chunk in sec relative to the
%beginning of the recording
%<<<Output<<<
%cropped_data = The cropped data for the specified time period

switch probe
    case 1
        channels = 385; %Total number of channels
    case 2
        channels = 385; %Total number of channels
    case 3
        channels = 385; %Total number of channels
    case 4
        channels = 277; %Total number of channels
end
fs = 30000; %Sampling frequncy
bytes_samp = 2; %Number of bytes per one sample

total_time = end_time - start_time; %Find the total time in sec that will be extracted

dir_info = dir([data_root '/*ap.bin']); %Get the names of all files with this extension
bin_filename = fullfile(dir_info.folder, dir_info.name); %Form the name of the file to be opened

total_time_samples = total_time.*fs; %Find the number of samples corresponding to that time
start_samples = start_time.*fs.*channels; %Find the number of samples that start is offset by
start_offset = start_samples.*bytes_samp; %Convert to bytes as fseek needs the number of bytes for the offset

fprintf('== Opening file %s ==\n',bin_filename);tic;
fid = fopen(bin_filename, 'r'); %Open an fid for the file of interest. This is something like a pointer for this file
offset = start_offset;
fseek(fid, offset, 'bof'); %Move the place to which fseek is pointing
cropped_data = fread(fid, [channels,total_time_samples],'*int16'); %Get the data for the specified time period from the bin file for all channels
fclose(fid); %Close the fid so it doesn't interfere with future calls of fopen
fprintf('== Done! Loading took %.1f sec ==\n',toc);
end