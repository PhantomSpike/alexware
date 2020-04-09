function plot_compare_fft(file_path1,file_path2,save_dir,ear)

font_axis = 12;
line_width = 2;
if strcmp(ear,'left')
    chan = 1;
elseif strcmp(ear,'right')
    chan = 2;
end

var_name1 = get_filename(file_path1);
var_name2 = get_filename(file_path2);


stim1 = load(file_path1);
stim2 = load(file_path2);
fs = stim1.Fs;

data1 = stim1.data(:,chan);
data2 = stim2.data(:,chan);
sz1 = length(data1);
sz2 = length(data2);

if sz1 > sz2
    data2(sz2:sz1) = 0;
else
    data1(sz1:sz2) = 0;
end

[P1,f] = fft_run(data1,fs,false);
[P2,f] = fft_run(data2,fs,false);
max_P = max([P1(:);P2(:)]);
figure('units','normalized','outerposition',[0 0 1 1]);
row = 1;
col = 2;
per = 0.005;
edgel = 0.03; edger = per; edgeh = 0.03; edgeb = 0.05; space_h = per; space_v = per;
[pos]=subplot_pos(row,col,edgel,edger,edgeh,edgeb,space_h,space_v);
subplot('position',pos{1});
plot(f,P1,'LineWidth',line_width);
ylim([0 max_P]);
legend([var_name1,' ',ear,' ear']);
xlabel('Freqeuncy [kHz]');
ylabel('Fourier magnitude');
set(gca,'FontName','Arial','FontSize',font_axis,'FontWeight','Bold');
subplot('position',pos{2});
plot(f,P2,'Color',[0.8500 0.3250 0.0980],'LineWidth',line_width);
ylim([0 max_P]);
set(gca, 'YTick', []);
legend([var_name2,' ',ear,' ear']);
xlabel('Freqeuncy [kHz]');
set(gca,'FontName','Arial','FontSize',font_axis,'FontWeight','Bold');
set(gcf,'color','w');
save_name = [save_dir,var_name2,'_vs_',var_name1,'_FFT_spectrum_',ear,'_ear.jpg'];
export_fig(save_name);
close all;