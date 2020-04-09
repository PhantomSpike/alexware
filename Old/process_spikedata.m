function process_spikedata(data_root,output_root,processed_name,start_time_s,end_time_s,chunk_size_s)
%process_spikedata(data_root,processed_name,start_time,end_time)
%This function loads, processes and saves binary spike data
fs = 30000;
total_time = end_time_s - start_time_s; %Find the total time in sec that will be extracted
num_chunks = total_time./chunk_size_s; %Find the number of chunks necessary to load the data

if ~(floor(num_chunks)==num_chunks)
    error('The total time in seconds should be a multiple of the chunk size in seconds');
end

%Define the start and end times of all chunks
chunk_start_times_s = [start_time_s:chunk_size_s:end_time_s - chunk_size_s];
chunk_end_times_s = chunk_start_times_s + chunk_size_s;

new_name = fullfile(output_root,processed_name);
new_name = [new_name,'.ap.bin'];
if exist(new_name, 'file')
    delete(new_name);
end
fid2 = fopen(new_name,'a');
zf = [];
for chunk = 1:num_chunks
    fprintf('== Chunk no %.0f/%.0f ==\n',chunk,num_chunks);
    %% Load spike data
    [cropped_data] = load_spikedata(data_root,chunk_start_times_s(chunk),chunk_end_times_s(chunk));
    %% Filter spike data
%     ch_385 = cropped_data(385,:);
%     cropped_data(385,:) = [];
%     cropped_data = cropped_data';
%     low_cutoff = 300;
%     high_cutoff = 6000;
%     stopband_wide_low = 100;
%     stopband_wide_high = 2500;
%     [cropped_data_filt, zf] = butter_filt(cropped_data,fs,low_cutoff,high_cutoff,stopband_wide_low,stopband_wide_high, zf);
%     cropped_data_filt = int16(cropped_data_filt);
%     cropped_data_filt = cropped_data_filt';
%     cropped_data_filt(385,:) = ch_385;
    %% CRA
    cropped_data_filt_cra = cras(cropped_data_filt);
    %% Write filtered data
    fprintf('== Writing file %s, chunk no %.0f  ==\n',new_name,chunk);tic;
    fwrite(fid2,cropped_data_filt_cra,'int16');
    fprintf('== Done! Writing took %.1f sec ==\n',toc); 
end
fid2 = fopen(new_name,'a');