function plot_compare_db(file_path1,file_path2,save_dir,ear)

window_ms = 5;
font_axis = 20;
line_width = 3;
if strcmp(ear,'left')
    chan = 1;
elseif strcmp(ear,'right')
    chan = 2;
end

slash_ix1 = find(file_path1 == '/', 1, 'last');
slash_ix2 = find(file_path2 == '/', 1, 'last');
var_name1 = file_path1(slash_ix1+1:end);
var_name2 = file_path2(slash_ix2+1:end);
var_name1 = strrep(var_name1,'_',' ');
var_name2 = strrep(var_name2,'_',' ');


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

[dB1,time] = db_conv_calc(data1,window_ms,fs,false);
[dB2,time] = db_conv_calc(data2,window_ms,fs,false);

figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
plot(time,dB1,'LineWidth',line_width);
plot(time,dB2,'--','LineWidth',line_width);
ylim([0 100]);
xlabel('Time [s]');
ylabel('Sound level [dB SPL]');
title([var_name2,' vs ',var_name1,' ',ear,' Ear dB Level across time (',num2str(window_ms),' ms window)']);
set(gca,'FontName','Arial','FontSize',font_axis,'FontWeight','Bold');
legend(var_name1,var_name2);
set(gcf,'color','w');
save_name = [save_dir,var_name2,'_vs_',var_name1,'_dB_level_across_time_(',num2str(window_ms),'ms)_',ear,'_ear.jpg'];
export_fig(save_name);
close all;