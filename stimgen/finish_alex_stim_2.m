stimlist = dir('ferret_stim*.mat');
type_sound = 'same';
lowfreq = 200; %The low frequency cut-off
highfreq = 20000; %The high freuqncy cut-off
tdt_50k = 48828; %The sampling rate of the TDT
filt_order = 8; %The filter order
numrepeats = 10; %The number of repeats 
db_target = 73; %How loud is the middle condition
atten_factor = 10; %How much to attenuate the sound so that there is no clipping 
noiselevel = 5; %How many dB louder than the actual sound the noise to be
total_length_s = 40; %The length of all stimuli in one presentation in s
ramp_length_s = 0.2; %The length of the ramp in s
num_rep = 40; %How many times to repeat the ramp
noiselen = 0.5;
noise_start_s = 36;
noise_timewindow_s = 3.5;
end_s = 40;

ahs = [0 0];
bhs = [0 0];
chs = [0 0];
while (ahs(1)~=7) || (bhs(1)~=7) || (chs(1)~=7)
    %make matrix which determines which stim has noise for each repeat
    hasnoise = zeros(numrepeats,length(stimlist));
    hs = [ 1 1 1 0 0 0];
    for r = 1:numrepeats
        for n=1:(length(stimlist)./2)
            nn = n+3.*(rand>0.5);
            hasnoise(r,nn) = 1;
        end
    end

    
    ens = randperm(10);
    hasnoise(ens(1:4),[1 4]) = 1;
    ens = randperm(10);
    hasnoise(ens(1:4),[2 5]) = 1;
    ens = randperm(10);
    hasnoise(ens(1:4),[3 6]) = 1;
    
    ahs = sum(hasnoise(:,[1 4]));
    bhs = sum(hasnoise(:,[2 5]));
    chs = sum(hasnoise(:,[3 6]));
end
ahs = sum(hasnoise(:,[1 4]));
bhs = sum(hasnoise(:,[2 5]));
chs = sum(hasnoise(:,[3 6]));

%Get the adjustment coefficients to make all sounds as loud as db_target
[coeffs,db] = calc_adj_coeff(db_target,total_length_s,ramp_length_s);

%main loop
for n=1:length(stimlist)
    yy = [];
    yy_50k = [];
    y_50k_bp = [];
    y_50k_bp_norm = [];
    y_50k_bp_norm_bp = [];
    y_50k_bp_norm_bp_ramp = [];
    y_50k_bp_norm_bp_ramp_norm = [];
    y_50k_bp_norm_bp_ramp_norm_atten = [];
    %load
    load(stimlist(n).name)
    %figure(1)
    %subplot(3,1,1)
    %spectrogram(data(:,channel),fs);
    %imagesc(cg)
    switch type_sound
        case 'same'
            if n == 1 || n == 4
                coeff = coeffs.anech_init;
                coeff_final = coeffs.anech_final;
            elseif n == 2 || n == 5
                coeff = coeffs.reverb1_init;
                coeff_final = coeffs.reverb1_final;
            elseif n == 3 || n == 6
                coeff = coeffs.reverb2_init;
                coeff_final = coeffs.reverb2_final;
            end
        case 'relative'
            coeff = coeffs.reverb2_init;
    end
    
    %add noise bursts
    for r = 1:numrepeats
        
        yy = data;
        
        %truncate to 40s
        yy((end_s*fs + 1):end,:)=[];

   
        
        noisepos = 0;
        snoisepos = 0;
        snoiselen = 0;
        rmsL = 0;
        rmsR = 0;
        rmsB = 0;
        dbL = 0;
        dbR = 0;
        dbB = 0;
        


        
        %Resample to tdt50k
        
        for channel = 1:2
            yy_50k(:,channel) = resample(yy(:,channel),tdt_50k,fs);
        end

        %Bandpass between lowfreq and highfreq
        for channel = 1:2
            y_50k_bp(:,channel) = bandpass(yy_50k(:,channel),tdt_50k,lowfreq,highfreq,filt_order);
        end
        
        
        y_50k_bp_norm = y_50k_bp*coeff;
        
        %get rms and db
        dbL = db_calc(y_50k_bp_norm(:,1));
        dbR = db_calc(y_50k_bp_norm(:,2));
        dbB = db_calc(y_50k_bp_norm(:));
        rmsB = rms(y_50k_bp_norm(:));
        
        hasnoise(r,n);
        if hasnoise(r,n)          
            %make +5dB noise
            snoiselen = round(tdt_50k.*noiselen);
            noise = sqrt(noiselevel).*rmsB.*randn(snoiselen,1); %think about frozen vs. not frozen
            
            %insert noise
            %withnoise = data;
            noisepos = noise_start_s + noise_timewindow_s.*rand;
            snoisepos = round(tdt_50k.*noisepos);
            
            for channel = 1:2
                y_50k_bp_norm(snoisepos:(snoisepos+snoiselen-1),channel) = noise;
            end
        end
        
        
        %Bandpass again vecause of the noise
        for channel = 1:2
            y_50k_bp_norm_bp(:,channel) = bandpass(y_50k_bp_norm(:,channel),tdt_50k,lowfreq,highfreq,filt_order);
        end
        
        %Make the cosine ramp envelope
        env = cosrampenv(total_length_s,ramp_length_s,tdt_50k);
        env = env';
        env = env(1:length(y_50k_bp_norm_bp),:); %Make the ramp exactly 40s
        
        
        %Multiply the sound to introduce the ramp
        y_50k_bp_norm_bp_ramp = y_50k_bp_norm_bp.*env;
        
        
        %Normalize the sound to be the right db
        y_50k_bp_norm_bp_ramp_norm = y_50k_bp_norm_bp_ramp*coeff_final;

        %Divide by 10 to decrease sound by 20dB so that there is no
        %clipping
        y_50k_bp_norm_bp_ramp_norm_atten = y_50k_bp_norm_bp_ramp_norm./atten_factor;

        ix_clip_L = find(abs(y_50k_bp_norm_bp_ramp_norm_atten(:,1))>=1); 
        ix_clip_R = find(abs(y_50k_bp_norm_bp_ramp_norm_atten(:,2))>=1);
        per_clip_L = (numel(ix_clip_L)./length(y_50k_bp_norm_bp_ramp_norm_atten)).*100;
        per_clip_R = (numel(ix_clip_R)./length(y_50k_bp_norm_bp_ramp_norm_atten)).*100;
        
        switch stimlist(n).name(14:20)
            case 'anech.m'
                roomnum = 0;
            case 'burrow1'
                roomnum = 1;
            case 'burrow5'
                roomnum = 2;
        end
        stimnum = num2str(stimlist(n).name(12));
        stimname = ['wav_files/' 'alexreverb.' num2str(roomnum) '.' num2str(stimnum) '.' num2str(r) '.wav'];
        audiowrite(stimname,y_50k_bp_norm_bp_ramp_norm_atten,tdt_50k,'BitsPerSample',24);
        alexreverbmetadata(r,n).stimname = stimname;
        alexreverbmetadata(r,n).originalstim = stimlist(n).name;
        alexreverbmetadata(r,n).hasnoise = hasnoise(r,n);
        alexreverbmetadata(r,n).noisepos = noisepos;
        alexreverbmetadata(r,n).snoisepos = snoisepos;
        alexreverbmetadata(r,n).noiselen = noiselen;
        alexreverbmetadata(r,n).snoiselen = snoiselen;
        alexreverbmetadata(r,n).fs = tdt50k;
        alexreverbmetadata(r,n).bitspersample = 24;
        alexreverbmetadata(r,n).bandpasslowfreq = lowfreq;
        alexreverbmetadata(r,n).bandpasshighfreq = highfreq;
        alexreverbmetadata(r,n).bandpassorder = filt_order;
        alexreverbmetadata(r,n).dbL = dbL;
        alexreverbmetadata(r,n).dbR = dbR;
        alexreverbmetadata(r,n).dbB = dbB;
        alexreverbmetadata(r,n).noiselevel = 5;
        alexreverbmetadata(r,n).per_clip_L = per_clip_L;
        alexreverbmetadata(r,n).per_clip_R = per_clip_R;
        alexreverbmetadata(r,n).db_target = db_target;
        alexreverbmetadata(r,n).type = type_sound;
        
    end

end
save('wav_files/alexreverbmetadata',alexreverbmetadata);
