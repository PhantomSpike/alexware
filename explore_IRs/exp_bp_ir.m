file_path = '/media/alex/900f118b-953c-405d-8127-f5f5914ad078/alex/Reverb_MSc_Matlab/IRs/All/Full_new/ferret_reverb1_full';
save_dir = '/home/alex/Desktop/Figures/IR_investigate/McDermoth_IRs/';
lowfreq = 200; %Low frequency cut-off
highfreq = 20000; %High frqeuncy cut-off
short_fft = lowfreq; %Whether to cut the FFT and which value to center on
f_window = 1000; %How wide is the freqeuncy window for the FFT
filter_order = 10;

load(file_path);
%Fix names
var_name = get_filename(file_path);
%Load meta info
meta.fs = Fs;
meta.lowf = num2str(lowfreq);
meta.highf = num2str(highfreq);
meta.var_name = var_name;
meta.save_dir = save_dir;
meta.short_fft = short_fft;
meta.f_window = f_window;
%Filter the IR
data_bp = filter_ir(data,Fs,lowfreq,highfreq,filter_order);
%Plot the comparison for left ear 
meta.ear = 'left';
plot_compare_irs2(data(:,1),data_bp(:,1),meta);
plot_compare_fft2(data(:,1),data_bp(:,1),meta);
%Plot the comparison for right ear 
meta.ear = 'right';
plot_compare_irs2(data(:,2),data_bp(:,2),meta);
plot_compare_fft2(data(:,2),data_bp(:,2),meta);
