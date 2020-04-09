function plot_compare_irs2(data1,data2,meta)

font_axis = 20;
fs = meta.fs;
sz = length(data1);
t = [1/fs:1/fs:sz/fs];
figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
plot(t,data1);
plot(t,data2 + 1.5);
hold off;
var_name = meta.var_name;
xlabel('Time [s]');
ylabel('Amplitude');
title([var_name,' normal vs bp filtered ',meta.lowf,'-',meta.highf,'Hz ',meta.ear,' Ear Impulse response plot']);
set(gca,'FontName','Arial','FontSize',font_axis,'FontWeight','Bold');
legend('Normal','bp filtered ',meta.lowf,'-',meta.highf,'Hz');
set(gcf,'color','w');
save_name = [meta.save_dir,var_name,' normal vs bp filtered ',meta.lowf,'-',meta.highf,'Hz ',meta.ear,'_ear.jpg'];
export_fig(save_name);
close all;