function plot_IR_sing(file_path,save_dir,file_type)

font_axis = 20;
lwidth = 1;
switch file_type
    case 'wav'
        [data,Fs] = audioread(file_path);
    case 'mat'
        load(file_path);
        data(:,2) = [];
end
sz = length(data);
t = [1/Fs:1/Fs:sz/Fs];
figure('units','normalized','outerposition',[0 0 1 1]);
plot(t,data,'LineWidth',lwidth);
var_name = get_names(file_path);
xlabel('Time [s]');
ylabel('Amplitude');
title([var_name, ' Impulse response plot']);
set(gca,'FontName','Arial','FontSize',font_axis,'FontWeight','Bold');
set(gcf,'color','w');
save_name = [save_dir,var_name,'.jpg'];
export_fig(save_name);
close all