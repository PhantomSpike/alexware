function d_bandpass = bandpass_iir(Fpass1,Fstop1,Fpass2,Fstop2,Astop1,Astop2,pb_ripple,fs)

d_bandpass = designfilt('bandpassiir', ...
  'StopbandFrequency1',Fstop1,'PassbandFrequency1', Fpass1, ...
  'PassbandFrequency2',Fpass2,'StopbandFrequency2', Fstop2, ...
  'StopbandAttenuation1',Astop1,'PassbandRipple', pb_ripple, ...
  'StopbandAttenuation2',Astop2, ...
  'DesignMethod','butter','SampleRate', fs);

fvtool(d_bandpass);
