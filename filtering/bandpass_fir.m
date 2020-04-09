function filt_data = bandpass_fir(data,Fstop1,Fpass1,Fpass2,Fstop2,Astop1,Apass,Astop2,fs)

% Fstop1 = 150;
% Fpass1 = 200;
% Fpass2 = 300;
% Fstop2 = 350;
% Astop1 = 65;
% Apass  = 0.5;
% Astop2 = 65;
% Fs = 1e3;


d_bandpass = designfilt('bandpassiir', ...
  'StopbandFrequency1',Fstop1,'PassbandFrequency1', Fpass1, ...
  'PassbandFrequency2',Fpass2,'StopbandFrequency2', Fstop2, ...
  'StopbandAttenuation1',Astop1,'PassbandRipple', Apass, ...
  'StopbandAttenuation2',Astop2, ...
  'DesignMethod','butter','SampleRate', fs);


fvtool(d_bandpass);

filt_data = filtfilt(d_bandpass,data);