function [dB,time] = db_conv_calc(data,t_window_ms,fs,plot_on)

if nargin < 4
    plot_on = false;
end

t_window_s = t_window_ms/1000;
Q = 2*(10.^-5);
dB = 20.*log10(sqrt(conv(data.^2,ones(round(t_window_s.*fs),1)./round(t_window_s.*fs),'same'))/(Q));
time = [1/fs:1/fs:numel(data)/fs];

if plot_on
    figure;
    plot(time,dB);
    xlabel('Time [s]');
    ylabel('dB');
end