% sort_dirs{1} = '/mnt/40086D4C086D41D0/Reverb_data/Ronnie/P05_Ronnie/P05-alex_reverb_wtih_noise_same';
sort_dirs{1} = '/mnt/40086D4C086D41D0/Reverb_data/Ronnie/P06_Ronnie/P06-alex_reverb_wtih_noise_same';
sort_dirs{2} = '/mnt/40086D4C086D41D0/Reverb_data/Ronnie/P13_Ronnie/P13-alex_reverb_wtih_noise_same';
sort_dirs{3} = '/mnt/40086D4C086D41D0/Reverb_data/Derry/P02_Derry/P02-alex_reverb_wtih_noise_same';
sort_dirs{4} = '/mnt/40086D4C086D41D0/Reverb_data/Derry/P03_Derry/P03-alex_reverb_wtih_noise_same';
sort_dirs{5} = '/mnt/40086D4C086D41D0/Reverb_data/Derry/P05_Derry/P05-alex_reverb_wtih_noise_same';
sort_dirs{6} = '/mnt/40086D4C086D41D0/Reverb_data/Derry/P08_Derry/P08-alex_reverb_wtih_noise_same';
sort_dirs{7} = '/mnt/40086D4C086D41D0/Reverb_data/Kilkenny/P05_Kilkenny/P05-alex_reverb_wtih_noise_same';
sort_dirs{8} = '/mnt/40086D4C086D41D0/Reverb_data/Kilkenny/P06_Kilkenny/P06-alex_reverb_wtih_noise_same';
num_pen = length(sort_dirs);
rec_dur = 40;
font_axis = 20;

mean_wfs = [];
for pen = 1:num_pen
    fprintf('== Processing pen %0.f/%0.f ==\n',pen,num_pen);
    file_path = [sort_dirs{pen},'/CRA/clust_info.mat'];
    load(file_path);
    mean_wfs = [mean_wfs;clust_info.mean_wfs(clust_info.good_units(:,2),:)];
    clear('clust_info');
    fprintf('== Done! ==\n');
end

[results] = fast_reg(mean_wfs);
[idx,C,X,inh_flag] = run_kmeans(results,false);
% inh_rate = mean(clust_info.no_spikes(idx==inh_flag))/(rec_dur*60);
% exc_rate = mean(clust_info.no_spikes(idx~=inh_flag))/(rec_dur*60);
figure('units','normalized','outerposition',[0 0 1 1]);
plot(X(idx==inh_flag,1),X(idx==inh_flag,2),'b.','MarkerSize',20)
hold on;
plot(X(idx~=inh_flag,1),X(idx~=inh_flag,2),'r.','MarkerSize',20)
plot(C(:,1),C(:,2),'kx',...
    'MarkerSize',15,'LineWidth',3)
legend('Putative inhibitory','Putative excitatory','Centroids',...
    'Location','NE')
title('Inhibitory/Excitatory cluster assignment and centroids');
xlabel('Spike width [ms]');
ylabel('$\frac{Peak}{Trough}$ Ratio','Interpreter','latex');
dim = [.68 .5 .3 .3];
% anot = annotation('textbox',dim,'String',sprintf('Putative inhibitory avg rate=%.1f (sp/s)\nPutative excitatory avg rate=%.1f (sp/s)',inh_rate,exc_rate),'FitBoxToText','on');
anot.FontSize = 12;
anot.FontWeight = 'bold';
set(gca,'FontName','Arial','FontSize',font_axis,'FontWeight','Bold');
hold off;
set(gcf,'color','w');
save_dir = ['/home/alex/Desktop/Figures/IR_investigate/I_E_plots/'];
save_name = [save_dir,'_Alldata_IE_Clustering.jpg'];
export_fig(save_name);
close all;

