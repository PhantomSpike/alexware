function filt_data = lowpass_fir(data,Fpass,Fstop,Apass,Astop,fs)

d_low = designfilt('lowpassiir', ...
  'PassbandFrequency',Fpass,'StopbandFrequency',Fstop, ...
  'PassbandRipple',Apass,'StopbandAttenuation',Astop, ...
  'DesignMethod','butter','SampleRate',fs);

fvtool(d_low);

filt_data = filtfilt(d_low,data);
