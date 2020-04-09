
low_th = 0.49;
high_th = 0.524;
ix = results.spike_width_ms>low_th & results.spike_width_ms<high_th;
figure;
hold on;
plot(mean_wfs(ix,:)','b');
