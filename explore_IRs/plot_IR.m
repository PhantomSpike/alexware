function plot_IR(file_path,save_dir)

font_axis = 20;
load(file_path);
sz = length(data);
t = [1/Fs:1/Fs:sz/Fs];
figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
left_ear = data(:,1);
right_ear = data(:,2);
% plot(t,left_ear);
% plot(t,right_ear + 1.5);
plot(t,right_ear);
hold off;
var_name = get_names(file_path);
xlabel('Time [s]');
ylabel('Amplitude');
title([var_name, ' Impulse response plot']);
set(gca,'FontName','Arial','FontSize',font_axis,'FontWeight','Bold');
% legend('Left ear','Right ear');
set(gcf,'color','w');
save_name = [save_dir,var_name,'.jpg'];
export_fig(save_name);
close all;