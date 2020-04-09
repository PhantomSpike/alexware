function plot_compare_irs(file_path1,file_path2,save_dir,ear,pow_leg)

font_axis = 20;

if nargin < 5
    pow_leg = false;
end

if strcmp(ear,'left')
    chan = 1;
elseif strcmp(ear,'right')
    chan = 2;
end


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

t = [1/fs:1/fs:length(data1)/fs];

figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
plot(t,data1);
plot(t,data2 + 1.5);
hold off;
var_name1 = get_filename(file_path1);
var_name2 = get_filename(file_path2);
xlabel('Time [s]');
ylabel('Amplitude');
title([var_name2,' vs ',var_name1,' ',ear,' Ear Impulse response plot']);
set(gca,'FontName','Arial','FontSize',font_axis,'FontWeight','Bold');
legend(var_name1,var_name2);

if pow_leg
    %Compute the power ratios
    [Lp,power_ratio] = comp_power(data1,data2);
    %Create the textbox
    dim = [.725 .5 .3 .3];
    anot = annotation('textbox',dim,'String',sprintf('L_p =10*log_{10}(^{P_1}/_{P_2}) = %.3f\n(^{P_1}/_{P_2})*100 = %.3f%%',Lp,power_ratio),'FitBoxToText','on');
    anot.FontSize = 12;
    anot.FontWeight = 'bold';
end

set(gcf,'color','w');
save_name = [save_dir,var_name2,'_vs_',var_name1,'_',ear,'_ear.jpg'];
export_fig(save_name);
close all;