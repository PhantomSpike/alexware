function results = create_sounds(ir_folder,save_dir)

lowfreq = 200; %The low frequency cut-off
highfreq = 20000; %The high freuqncy cut-off
filt_order = 8; %The order of the filter
tdt_50k = 48828; %The sampling rate of the TDT
db_target = 73; %The target dB 
aten_factor = 10;
ramp_length_s = 0.2;
sound_path = '/home/alex/Desktop/Figures/IR_investigate/Sounds/DRC/original/drc1_dur25_contrast20.wav';

[sound_orig,fs] = audioread(sound_path);

ir_list = dir(ir_folder);
ir_list = ir_list(~ismember({ir_list.name},{'.','..'}));

num_ir = length(ir_list);

sz = length(sound_orig);
total_length_s = sz/fs;
%Make the cosine ramp envelope
env = cosrampenv(total_length_s,ramp_length_s,tdt_50k);
env = env';
tic;
for jj = 1:num_ir
    fprintf('== Processing sound %0.f/%0.f ==\n',jj,num_ir);
    file_path = fullfile(ir_list(jj).folder,ir_list(jj).name);
    slash_ix = find(file_path == '/', 1, 'last');
    var_name = file_path(slash_ix+1:end-4);
    ir_bank{jj} = load(file_path);
    raw_snd = [];
    resamp_snd = [];
    bp_resamp_snd = [];
    ramp_bp_resamp_snd = [];
     adj_coeff = [];
    for ear = 1:2
        %Convolve with the correct IR
        raw_snd(:,ear) = conv(sound_orig,ir_bank{jj}.data(:,ear),'same');
        %Resample to tdt50k
        resamp_snd(:,ear) = resample(raw_snd(:,ear),tdt_50k,fs);
        %Bandpass between lowfreq and highfreq
        bp_resamp_snd(:,ear) = bandpass(resamp_snd(:,ear),tdt_50k,lowfreq,highfreq,filt_order);
        %Bandpass between lowfreq and highfreq
        ramp_bp_resamp_snd(:,ear) = env.*bp_resamp_snd(:,ear);
    end
    adj_coeff = db_adjust(ramp_bp_resamp_snd(:),db_target);
    results(jj).sound = ramp_bp_resamp_snd.*adj_coeff;
    results(jj).sound = results(jj).sound./aten_factor;
    results(jj).dB_level = db_calc(results(jj).sound(:));
    results(jj).name = var_name;
    results(jj).ir = ir_bank{jj};
    stimname = [save_dir,var_name,'.wav'];
    audiowrite(stimname,results(jj).sound,tdt_50k,'BitsPerSample',24);
end
save([save_dir,'results.mat'],'results');
fprintf('== Done! Processing took %1.fs ==\n',toc);

