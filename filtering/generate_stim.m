fprintf('== Processing the file ==\n');
%Load the file
sound_file = '/home/alex/Desktop/Cricket/crickets.wav'; %Path to the file
[data,fs] = audioread(sound_file);
%Choose one of the channels
ear = 'left';

switch ear
    case 'left'
        data = data(:,1);
    case 'right'
        data = data(:,2);
end

%Plot the original sound file
figure;
example_coch(data,fs);
figure;
[P1,f] = fft_run(data,fs);
%Make a high-pass filter infinite response filter
Fstop = 1000; %Where the filter starts in Hz
Fpass = 1100; %Where the filter ends in Hz
Astop = 60; %How much attenuation there is from Fpass to Fstop in dB 
pb_ripple = 0.5; %How much ringing is acceptable in the passband in dB

%Generate the filter with these characteristics 
d_highpass = highpass_iir(Fpass,Fstop,Astop,pb_ripple,fs);

%Filter the data with zero-order filtering compensating for phase shifts
fprintf('== Performing filtering ==\n');tic;
filt_data = filtfilt(d_highpass,data);
fprintf('== Done! Filtering took %1.fs ==\n',toc);

%Apply a cosine ramp to avoid clicks at the beginning and end of the
%recording 
envelope_length = numel(data)/fs; %Find the length of the sound file in seconds 
ramp_length = 0.01; %Define the ramping time in seconds
ramp_envelope = cosrampenv(numel(data)/fs,0.5,fs); %Generate the ramp
ramp_filt_data = filt_data.*ramp_envelope; %Apply the ramp

%Write as a wav file with twice the sampling rate so that all freqeuncies
%are made twice higher

double_fs = 2*fs;
file_name = 'final_crickets.wav';
audiowrite('final_crickets.wav',ramp_filt_data,double_fs);

%Open the saved file and plot Freqeuncy spectrum and spectrogram to check
%that it is correct
[test_data,test_fs] = audioread('final_crickets.wav');

figure;
example_coch(test_data,test_fs);
figure;
[P1,f] = fft_run(test_data,test_fs);