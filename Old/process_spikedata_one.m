function process_spikedata_one(data_root,output_root,processed_name,start_time_s,end_time_s)
%process_spikedata(data_root,processed_name,start_time,end_time)
%This function loads, processes and saves binary spike data
fs = 30000;
new_name = fullfile(output_root,processed_name);
new_name = [new_name,'.ap.bin'];
if exist(new_name, 'file')
    delete(new_name);
end
fid2 = fopen(new_name,'w');
%% Load spike data
[cropped_data] = load_spikedata(data_root,start_time_s,end_time_s);
%% Filter spike data
ch_385 = cropped_data(385,:);
cropped_data(385,:) = [];
cropped_data = cropped_data';
low_cutoff = 300;
high_cutoff = 6000;
stopband_wide_low = 100;
stopband_wide_high = 2500;
[cropped_data_filt] = butter_filt_one(cropped_data,fs,low_cutoff,high_cutoff,stopband_wide_low,stopband_wide_high);
cropped_data_filt = int16(cropped_data_filt);
cropped_data_filt = cropped_data_filt';
cropped_data_filt(385,:) = ch_385;
%% CRA
cropped_data_filt_cra = cras(cropped_data_filt);
%% Write filtered data
fprintf('== Writing file %s ==\n',new_name);tic;
fwrite(fid2,cropped_data_filt_cra,'int16');
fprintf('== Done! Writing took %.1f sec ==\n',toc);