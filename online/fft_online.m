function fft_online(start_time,end_time,data_root)
%% Load selected chunk of the data
fs = 30000;
plot_avg = true; %Whether to plot average results
plot_all = true; %Whether to plot all results

[cropped_data,time_samples] = load_spikedata(data_root,start_time,end_time);
cropped_data = cropped_data'; %Flip the data so each column is one channel and each row one sample
%% Perform FFT on the data
[P1,f] = fft_run(cropped_data,time_samples,fs);
%% Plot the average values for channels
if plot_avg
    chan_no = (1:200); %Which channels to use for the averaging
    top_freq = (1:30); %How many frequencies to use, in this case the val of the top 20 frequencies
    num_freq = num2str(length(top_freq));
    
    med_val = median(P1(:,chan_no),2); %Median value for the frequencies' power for the selected channels
    [sorted_val, sorted_ix] = sort(med_val,'descend'); %Sort these median values in descending order and get their indices
    freq = f(sorted_ix(top_freq));
    val = sorted_val(top_freq);
    figure;
    bar(freq,val);
    title(['Unfiltered data',' - Top ',num_freq,' frequencies of channels ',num2str(chan_no(1)),' to ',num2str(chan_no(end))]);
    xlabel('f(KHz)');
    ylabel('|P1(f)|');
    xlim([0 15000]);
    xticks((0:1000:15000));
    xticklabels({'0','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15'});
    figure;
    plot(f,med_val);
    xticks((0:1000:15000));
    xticklabels({'0','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15'});
%     ylim([0 150]);
    title(['Unfiltered data',' - Median value of channels ',num2str(chan_no(1)),' to ',num2str(chan_no(end))]);
    xlabel('f(KHz)');
    ylabel('|P1(f)|');
    fprintf('The top %s frequencies in Hz are: \n',num_freq);
    fprintf('%.2f \n',freq);
end
%% Plot the single-sided amplitude spectrum P1
if plot_all
    ch_no = (1:20:300); %Which channels to plot
    per = 0.027;
    edgel = per; edger = per; edgeh = per; edgeb =per; space_h = 0.03; space_v =0.06;
    [pos]=subplot_pos(3,5,edgel,edger,edgeh,edgeb,space_h,space_v);
    figure;
    num_ch = length(ch_no);
    for ii = 1:num_ch
        subplot('position',pos{ii});
        plot(f,P1(:,ch_no(ii)));
        xticks((0:1000:15000));
        xticklabels({'0','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15'});
        %         ylim([0 150]);
        title(['Channel no ', num2str(ch_no(ii))]);
        xlabel('f(KHz)');
        ylabel('|P1(f)|');
    end
    figure;
    channel = 77;
    plot(f,P1(:,channel));
    xticks((0:1000:15000));
    xticklabels({'0','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15'});
    %         ylim([0 150]);
    title(['Channel no ', num2str(channel)]);
    xlabel('f(KHz)');
    ylabel('|P1(f)|');
end
end