 function synch_ch = get_synch(data_root) 

%Load the meta info
%Load the Spike GLX Meta file
meta_dir = dir(fullfile(data_root, '*ap*.meta')); % meta file from spikeGLX specifically
meta_info = readSpikeGLXmeta(fullfile(meta_dir.folder, meta_dir.name));

%Select the total number of channels based on the probe option 
probe = meta_info.imProbeOpt;
switch probe
    case 1
        nChansTotal = 385; %Total number of channels
    case 2
        nChansTotal = 385; %Total number of channels
    case 3
        nChansTotal = 385; %Total number of channels
    case 4
        nChansTotal = 277; %Total number of channels
end

%% Load the synch channel
fprintf('== Loading the synch channel ==\n');tic;
dir_info = dir([data_root '/*ap.bin']); %Get the names of all files with this extension
bin_filename = fullfile(dir_info.folder, dir_info.name); %Form the name of the file to be opened

bytes_sample = 2;
nSampsTotal = dir_info.bytes/nChansTotal/bytes_sample; %Find the total number of samples in the file to be analyzed

offset_synch = (nChansTotal-1)*bytes_sample; %This should move the pointer to the 385 channel i.e. the synch channel

fid = fopen(bin_filename, 'r'); %Open an fid for the file of interest. This is something like a pointer for this file

fseek(fid, offset_synch, 'bof'); %Move the place to which fseek is pointing

precision = '*int16';
skip = (nChansTotal-1)*bytes_sample;
synch_ch = fread(fid,[nSampsTotal,1],precision,skip);

fclose(fid); %Close the fid so it doesn't interfere with future calls of fopen
save(fullfile(data_root,'synch_ch.mat'),'synch_ch');
fprintf('== Done! Loading took %.0fs ==\n',toc);
