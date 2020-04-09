function example_coch(data,fs,time_step)

fmin = 200;
fmax = fs/2;
total_no_freq = 70;
freq_step = 3;
dt = 5;
[X_ft, t, params] = cochleagram(data, fs, dt, 'log', fmin, fmax, total_no_freq);
figure('units','normalized','outerposition',[0 0 1 1]);
imagesc(flipud(X_ft));

time_s = t(1:time_step:numel(t));
xticks(1:time_step:numel(t));
for jj = 1:numel(time_s)
    x_labels{jj} = num2str(time_s(jj),'%.2f');
end
xticklabels(x_labels);

yticks(1:freq_step:numel(params.f));
freqs = params.f(1:freq_step:numel(params.f))/1000; %Convert the frequencies into kHz
for jj = 1:numel(freqs)
    y_labels{jj} = num2str(freqs(jj),'%.2f');
end
spacing = 2;
yticklabels(fliplr(y_labels));
set(gca,'FontName','Arial','FontSize',14,'FontWeight','Bold');
title('Cochleagram of Stimulus','FontSize',24,'FontWeight','bold');
xlabel('Time [s]','FontSize',22,'FontWeight','bold');
ylabel('Frequency [kHz]','FontSize',22,'FontWeight','bold');
colorbar;
colormap('jet');