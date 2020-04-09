function y = bandpass(x,Fs,clow,chigh,steepness)
% y = lowpass(x,Fs,clow,chigh);
%bandpass butterworth filtering of data x, with sampling rate Fs, and
%frequency cutoff clow and chigh in Hz, and the filter order (default 8)
if nargin<5
steepness = 8;
end
nyquist = Fs./2;
%c/nyquist
%[clow chigh]/nyquist
[b,a] = butter(steepness,[clow chigh]/nyquist);
y = filter(b,a,x);
%freqz(b,a,1024,Fs)
