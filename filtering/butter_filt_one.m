function [dataOut] = butter_filt_one(dataIn,fs,low_cutoff,high_cutoff,stopband_wide_low,stopband_wide_high)
%dataOut = butter_filt(dataIn,fs)
%This function performs a Fast-Fourier Tranform on a given input data and
%returns the single sided spectrum P1 along with its frequncy domain values
%>>>Input>> 
% stopband_wide_low = How wide the stopband is on the lower end of the spectrum in Hz
% stopband_wide_high = How wide the stopbands are on the higher end of the spectrum in Hz
%fs = Sampling frequncy in Hz
%<<<Output<<<
%P1 = The single-sided spectrum P1 which is derived from P2. This is the y-axis of the fft graph. Because P2 is
%symmetric we can take 1/2 without losing information.
%f = The frequnecies which feature in the fft analysis. This is the x-axis
%of the fft graph.
fprintf('== Performing Butterworth filtering ==\n');tic;
Nyquist = fs/2; %Find the Nyquist frequency
Rp = 3; %Maximum passband ripple in dB
Rs = 30; %Minimum attenuation in stopbands in dB
Wp = [low_cutoff high_cutoff]/Nyquist;
Ws_low = low_cutoff-stopband_wide_low;
Ws_high = high_cutoff+stopband_wide_high;
Ws = [Ws_low  Ws_high]/Nyquist;
[n,Wn] = buttord(Wp,Ws,Rp,Rs);
[b,a] = butter(n,Wn);
[dataOut] = filter(b,a,dataIn);
fprintf('== Done! Filtering took %f sec ==\n',toc);
end