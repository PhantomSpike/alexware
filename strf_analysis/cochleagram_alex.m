function [coch,t_s] = cochleagram_alex(sound_files_root,time_bin_ms,n_h)

dir_info = dir([sound_files_root '/*.wav']); %Get the names of all files with this extension
ix_files = [1,11,21,31,41,51];
ear = 2; %Choose 1 for left and 2 for right
type_coch = "log"; %The type of cochleagram we want
f_min = 500; %The minmum freqeuncy in Hz
f_max = 22627; %The maximum freqeuncy in Hz
f_steps = 20; %Number of freqeuncy steps that we want
end_time_s = 36;

for stim = 1:numel(ix_files)
    ix = ix_files(stim);
    filename = fullfile(dir_info(1).folder, dir_info(ix).name); %Form the name of the file to be opened
    coch(stim).filename = filename;
    [data,fs] = audioread(filename);
    end_time_bin = floor(end_time_s*fs);
    input_data = data(1:end_time_bin,ear);
    [X_ft, t_s, params] = cochleagram(input_data, fs, time_bin_ms, type_coch, f_min, f_max, f_steps);
    params.n_h = n_h;
    params.time_bin_ms = time_bin_ms;
    params.end_time_s = end_time_s;
    coch(stim).params = params;
    X_fht = tensorize(X_ft, n_h); % n_h = number of history steps desired
    coch(stim).X_fht = X_fht(:,:,n_h:end);
end

coch(1).anech = cat(3,coch(1).X_fht,coch(2).X_fht);
coch(1).reverb1 = cat(3,coch(3).X_fht,coch(4).X_fht);
coch(1).reverb2 = cat(3,coch(5).X_fht,coch(6).X_fht);
% save('/media/phant0msp1ke/4TB SSD/DATA/Ferret_exp_17_10_2017/Penetration2/P02-reverb_with_noise/CRA/P02-reverb_with_noise_cra_allch_3000Hz/cochleagram/coch','coch');
end