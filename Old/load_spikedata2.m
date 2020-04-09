function [cropped_data,total_time_samples] = load_spikedata2(data_root)
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

dir_info = dir([data_root '/*ap.bin']); %Get the names of all files with this extension
bin_filename = fullfile(dir_info.folder, dir_info.name); %Form the name of the file to be opened
nChansTotal = 385; %Total number of channels
bytes_samp = 2; %Number of bytes per one sample. This is for an int16 data type
fs = 30000; %Sampling rate in Hz
%Define the size of one chunk in samples
chunkSize = 1000000; %This is the size recommended by Nick Steinmetz
%Find the total number of samples in the file to be analyzed
nSampsTotal = dir_info.bytes/nChansTotal/bytes_samp;
%Find the number of chunks necessary to process the file
nChunksTotal = ceil(nSampsTotal/chunkSize);
%This will be the number of chunks which have the standard size
normal_chunks = nChunksTotal - 1;
%Find out the size of the last chunk which will be normalChunksize + extra samples left
final_chunkSize = floor((nSampsTotal/chunkSize - normal_chunks)*chunkSize);

fprintf('== Opening file %s ==\n',bin_filename);tic;
fid = fopen(bin_filename, 'r'); %Open an fid for the file of interest. This is something like a pointer for this file
offset = start_offset;
fseek(fid, offset, 'bof'); %Move the place to which fseek is pointing
cropped_data = fread(fid, [channels,total_time_samples],'*int16'); %Get the data for the specified time period from the bin file for all channels
fclose(fid); %Close the fid so it doesn't interfere with future calls of fopen
fprintf('== Done! Loading took %.1f sec ==\n',toc);
end