ferret = 1;
root_path = '/home/alex/Desktop/Code/roomsim/Impulse_response';
% ir_file = '/mnt/40086D4C086D41D0/Reverb_analysis/IRs/Open_IRs/exp_pos/ferret_small.mat';
file_path = '/mnt/40086D4C086D41D0/Reverb_analysis/IRs/Open_IRs/exp_pos/ferret_big.mat';
% file_path = fullfile(root_path,ir_file);
% file_path = '/media/alex/900f118b-953c-405d-8127-f5f5914ad078/alex/Reverb_MSc_Matlab/IRs/All/Full_new/ferret_reverb2_full.mat';
save_dir = '/home/alex/Desktop/Figures/IR_investigate/Ferret_fest/';

if ferret
    load(file_path);
    data(:,2) = [];
    file_type = 'mat';
else
    [data,Fs] = audioread(file_path);
    file_type = 'wav';
end
%Plot the IR
plot_IR_sing(file_path,save_dir,file_type);

%Plot the FFT
fc = 1250;
fw = 1250;
plot_on = 1;
[P1,f] = fft_run(data,Fs,plot_on,fc,fw);
var_name = get_names(file_path);
save_name = [save_dir,var_name,'_FFT_plot.jpg'];
export_fig(save_name);
close all;

%Plot the ACC
t_ms = 50;
font_axis = 20;
no_time_lags = round((t_ms/1000)*Fs);
[acf,time_lags] = autocorr(data,no_time_lags);
time_lag_ms = (time_lags/Fs)*1000;
figure('units','normalized','outerposition',[0 0 1 1]);
stem(time_lag_ms,acf,'filled');
xlabel('Time [ms]');
ylabel('Autocorrelation (r)');
title([var_name, ' Autocorrelogram of the IR']);
set(gca,'FontName','Arial','FontSize',font_axis,'FontWeight','Bold');
set(gcf,'color','w');
save_name = [save_dir,var_name,'ACC_of_IR.jpg'];
export_fig(save_name);
close all;

%Plot the kurtosis
t_win_ms = 10;
kurt = ir_kurtosis(data,Fs,t_win_ms);
t_s = [1/Fs:1/Fs:length(kurt)/Fs];
figure('units','normalized','outerposition',[0 0 1 1]);
plot(t_s,kurt,'LineWidth',3);
line([0 length(kurt)/Fs],[3 3],'Color','red','LineStyle','--','Linewidth',3);
xlabel('Time [s]');
ylabel('Kurtosis');
title([var_name, ' Kurtosis of the IR']);
set(gca,'FontName','Arial','FontSize',font_axis,'FontWeight','Bold');
set(gcf,'color','w');
save_name = [save_dir,var_name,'Kurtosis_of_IR.jpg'];
export_fig(save_name);
close all;


