function fft_spikedata(start_time,end_time)
%% Load selected chunk of the data
fs = 30000;
data_root = '/data/alex/exp_10_04_2017_sylvester/P03-QCPitch';
plot_avg = true; %Whether to plot average results
plot_all = true; %Whether to plot all results
export = false; %Whether to save the figures
filter_step = true;%Whether to filter the data

save_path = '/home/phant0msp1ke/Desktop/Figures/FFT_Spike_Data/'; %Where to save the figures
processed_name = 'P01-QGPitch-az_fft_300_600.png';

[cropped_data,time_samples] = load_spikedata(data_root,start_time,end_time);
cropped_data = cropped_data'; %Flip the data so each column is one channel and each row one sample
%% Extra filtering step
if filter_step
    low_cutoff = 300;
    high_cutoff = 4000;
    stopband_wide_low = 100;
    stopband_wide_high = 1000;
    cropped_data_filt = butter_filt(cropped_data,fs,low_cutoff,high_cutoff,stopband_wide_low,stopband_wide_high);
end
%% Perform FFT on the data
[P1,f] = fft_run(cropped_data_filt,time_samples,fs);
[P_unfilt,f] = fft_run(cropped_data,time_samples,fs);
%% Plot the average values for channels
if plot_avg
    chan_no = (1:200); %Which channels to use for the averaging
    top_freq = (1:20); %How many frequencies to use, in this case the val of the top 20 frequencies
    
    med_val = median(P1(:,chan_no),2); %Median value for the frequencies' power for the selected channels
    med_val_unfilt = median(P_unfilt(:,chan_no),2); %Median value for the frequencies' power for the selected channels
    [sorted_val, sorted_ix] = sort(med_val,'descend'); %Sort these median values in descending order and get their indices
    [sorted_val_unfilt, sorted_ix_unfilt] = sort(med_val_unfilt,'descend'); %Sort these median values in descending order and get their indices
    freq = f(sorted_ix(top_freq));
    freq_unfilt = f(sorted_ix_unfilt(top_freq));
    val = sorted_val(top_freq);
    val_unfilt = sorted_val_unfilt(top_freq);
    figure;
    bar(freq,val);
    figure;
    bar(freq_unfilt,val_unfilt);
    figure;
    plot(f,med_val);
    xticks((0:1000:15000));
    xticklabels({'0','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15'});
    ylim([0 150]);
    title(['FIltered data','Median value of channels ',num2str(chan_no(1)),' to ',num2str(chan_no(end))]);
    xlabel('f(KHz)');
    ylabel('|P1(f)|');
    figure;
    plot(f,med_val_unfilt);
    xticks((0:1000:15000));
    xticklabels({'0','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15'});
    ylim([0 150]);
    title(['Unfiltered data','Median value of channels ',num2str(chan_no(1)),' to ',num2str(chan_no(end))]);
    xlabel('f(KHz)');
    ylabel('|P1(f)|');
end
%% Plot the single-sided amplitude spectrum P1
if plot_all
    ch_no = (1:20:300); %Which channels to plot
    per = 0.025;
    edgel = per; edger = per; edgeh = per; edgeb =per; space_h = 0.025; space_v =0.05;
    [pos]=subplot_pos(3,5,edgel,edger,edgeh,edgeb,space_h,space_v);
    figure;
    %     figure('units','normalized','outerposition',[0 0 1 1]);
    for ii = 1:15
        subplot('position',pos{ii});
        plot(f,P1(:,ch_no(ii)));
        xticks((0:1000:15000));
        xticklabels({'0','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15'});
        ylim([0 150]);
        title(['Channel no', num2str(ch_no(ii))]);
        xlabel('f(KHz)');
        ylabel('|P1(f)|');
    end
    if export
        save_name = [save_path,processed_name];
        export_fig(save_name,'-r200');
    end
end
end