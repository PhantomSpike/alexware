function [P1,f] = fft_run(data,fs,plot_on,fc,fw)
%[P1,f] = fft_run(input_data,time_samples)
%This function performs a Fast-Fourier Tranform on a given input data and
%returns the single sided spectrum P1 along with its frequncy domain values
%>>>Input>> 
%input_data = This is the data that you want to perform the fft on.
%input_data must be a vector or an array. If it is an array, a separate fft
%is computed along each column of the array
%time_samples = These are the time samples of the input data
%fs = Sampling frequncy in Hz
%<<<Output<<<
%P1 = The single-sided spectrum P1 which is derived from P2. This is the y-axis of the fft graph. Because P2 is
%symmetric we can take 1/2 without losing information.
%f = The frequnecies which feature in the fft analysis. This is the x-axis
%of the fft graph.

lwidth = 1;

if ~exist('plot_on','var')
    plot_on = false;
end

if ~exist('fc', 'var')
    fc = 0;
end

if ~exist('fw', 'var')
    fw = 100;
end

time_samples = length(data);
fprintf('== Performing FFT ==\n');tic;
fft_data = fft(data); %Peform FFT on each channel for the whole duration of the chunk
%Compute the two-sided spectrum P2. Then compute the single-sided spectrum P1 based on P2 and the even-valued total number of samples.
P2 = abs(fft_data./time_samples);
clear fft_data; %As this can be quite large there is no point to keep it in memory
P1 = P2(1:time_samples/2+1,:);
clear P2;
P1(2:end-1,:) = 2.*P1(2:end-1,:);

f = fs*(0:(time_samples/2))/time_samples; %Define the frequency domain f
f = f/1000;

if fc ~= 0
    [~, ix_match] = min(abs(f-fc/1000));
    window_ix = fw./(mean(diff(f))*1000); %Find the ix for the frqeuncy plotting window
    f = f(ix_match-window_ix:ix_match+window_ix);
    P1 = P1(ix_match-window_ix:ix_match+window_ix,:);
end

if plot_on
    figure('units','normalized','outerposition',[0 0 1 1]);
    plot(f,P1,'LineWidth',lwidth);
    set(gca,'FontName','Arial','FontSize',10,'FontWeight','Bold');
    xlabel('Freqeuncy [kHz]','FontSize',16,'FontWeight','bold');
    ylabel('|P1(f)|','FontSize',16,'FontWeight','bold');
    title(['Single-Sided Amplitude Spectrum of X(t)'],'FontSize',16,'FontWeight','bold');
end
fprintf('== Done! FFT took %.1f sec ==\n',toc);
end