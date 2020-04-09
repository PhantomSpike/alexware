function plot_IR2(data,Fs)

font_axis = 20;
sz = length(data);
t = [1/Fs:1/Fs:sz/Fs];
figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
left_ear = data(:,1);
right_ear = data(:,2);
plot(t,left_ear);
plot(t,right_ear + 1.5);
hold off;
xlabel('Time [s]');
ylabel('Amplitude');
% title([var_name, ' Impulse response plot']);
set(gca,'FontName','Arial','FontSize',font_axis,'FontWeight','Bold');
legend('Left ear','Right ear');
set(gcf,'color','w');
