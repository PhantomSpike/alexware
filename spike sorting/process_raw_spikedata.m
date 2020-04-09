function process_raw_spikedata(fpath,ChanTotal) 

processed_folder_name = '/CRA';
%This function loads road raw spike data from a SpikeGLX which is in .bin
%form. It then does common reference averaging (CRA) and writes the
%processed data to a new .bin file.
%>>> Input >>>
%data_root = The path to your folder of interest
%output_root = The path to the output folder
%processed_name = The name of the porcessed file

[PathStr,FolderName]=fileparts(fpath);
processed_name = [FolderName,'_cra'];

dir_info = dir([fpath '/*ap.bin']); %Get the names of all files with this extension
bin_filename = fullfile(dir_info.folder, dir_info.name); %Form the name of the file to be opened

processed_folder_fullname = [fpath,processed_folder_name];

if exist(processed_folder_fullname, 'dir')
    rmdir(processed_folder_fullname);
end

mkdir(processed_folder_fullname);

nChansTotal = ChanTotal; %Total number of channels

bytes_samp = 2; %Number of bytes per one sample. This is for an int16 data type

fs = 30000; %Sampling rate in Hz

chunkSize = 1000000; %Define the size of one chunk in samples. This is the size recommended by Nick Steinmetz and is equa to 33 sec of data

nSampsTotal = dir_info.bytes/nChansTotal/bytes_samp; %Find the total number of samples in the file to be analyzed

nChunksTotal = ceil(nSampsTotal/chunkSize); %Find the number of chunks necessary to process the file

normal_chunks = nChunksTotal - 1; %This will be the number of chunks which have the standard size

final_chunkSize = floor((nSampsTotal/chunkSize - normal_chunks)*chunkSize); %Find out the size of the last chunk which will be normalChunksize + extra samples left

fprintf('== Opening file %s ==\n',bin_filename);
fid = fopen(bin_filename, 'r'); %Open an fid for the file of interest. This is something like a pointer for this file
offset = 0;
fseek(fid, offset ,'bof'); %Move the pointer to the beginning of the file

new_name = fullfile(processed_folder_fullname,processed_name);
new_name = [new_name,'imec.ap.bin'];
if exist(new_name, 'file')
    delete(new_name);
end
fid2 = fopen(new_name,'a');

for chunk = 1:nChunksTotal
    fprintf('== Chunk no %.0f/%.0f ==\n',chunk,nChunksTotal);
    %% Load spike data
    if chunk == nChunksTotal
        chunkSize = final_chunkSize; %When the final chunk has to be loaded adjust the size 
    end
    cropped_data = fread(fid, [nChansTotal,chunkSize],'*int16'); %Get the data for the specified time period from the bin file for all channels
    %% CRA
    cropped_data_cra = crefavg(cropped_data);
    %% Write processed data
    fprintf('== Writing file %s, chunk no %.0f  ==\n',new_name,chunk);tic;
    fwrite(fid2,cropped_data_cra,'int16');
    fprintf('== Done! Writing took %.1f sec ==\n',toc); 
    
    fseek(fid, offset ,'cof'); %Move the pointer to where the last one finished 
end

fclose(fid);
fclose(fid2);

end