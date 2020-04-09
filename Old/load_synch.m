function synch_ch = load_synch(data_root,start_time_s,end_time_s,chunk_s)

total_time = end_time_s - start_time_s; %Find the total time in sec that will be extracted
num_chunks = total_time./chunk_s; %Find the number of chunks necessary to load the data

if ~(floor(num_chunks)==num_chunks)
    error('The total time in seconds should be a multiple of the chunk size in seconds');
end

%Define the start and end times of all chunks
chunk_start_times_s = [start_time_s:chunk_s:end_time_s - chunk_s];
chunk_end_times_s = chunk_start_times_s + chunk_s;

parfor chunk = 1:num_chunks
    fprintf('== Processing chunk %.0f/%.0f ==\n',chunk,num_chunks);
    [cropped_data] = load_spikedata(data_root,chunk_start_times_s(chunk),chunk_end_times_s(chunk));
    synch_ch(:,chunk) = cropped_data(385,:);
end

synch_ch = synch_ch(:);
    