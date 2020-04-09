function plot_extr_ie(sorted_dir)

%% Define params
feature1 = 1;
feature2 = 2;
features = [feature1 feature2];
%1 - Spike width (peak-to-trough width), this is the distance form peak (depolarization) to
%trough (hyperpolarization) in the wf, var: spike_width_ms
%2 - Width at 1/2 height trough, this is the width of the trough part at 1/2
%max value, var: width_hh_trough_ms
%3 - End slope, this is the slope of the decaying part of the trough, var: end_slope
%4 - Peak-to-trough ratio, this is the ratio of the amplitudes of the peak
%to the trough PA/TA, var: peak_trough_ratio
%5 - Width at 1/2 height peak, this is the width of the peak part at 1/2
%max value, var: width_hh_peak_m
font_axis = 20; %Font size for plotting
%% Load and name
file_path = fullfile(sorted_dir,'clust_info.mat');
load(file_path);
slash_ix1 = find(file_path == '/', 1, 'last');
var_name = file_path(slash_ix1+1:end);
var_name = strrep(var_name,'_',' ');
%% Run k-means
[clust_info] = run_kmeans(clust_info,features);
inh_flag = clust_info.clustering.inh_flag;
X_norm = clust_info.clustering.X_norm;
C_norm = clust_info.clustering.C_norm;
x_lab_norm = clust_info.clustering.x_lab_norm;
y_lab_norm = clust_info.clustering.y_lab_norm;
idx = clust_info.clustering.idx;
figure('units','normalized','outerposition',[0 0 1 1]);
plot(X_norm(idx==inh_flag,1),X_norm(idx==inh_flag,2),'b.','MarkerSize',20)
hold on;
plot(X_norm(idx~=inh_flag,1),X_norm(idx~=inh_flag,2),'r.','MarkerSize',20)
plot(C_norm(:,1),C_norm(:,2),'kx',...
    'MarkerSize',15,'LineWidth',3)
legend('Putative inhibitory','Putative excitatory','Centroids',...
    'Location','NE')
title('Inhibitory/Excitatory cluster assignment and centroids [Normalized]');
xlabel(x_lab_norm);
ylabel(y_lab_norm);
set(gca,'FontName','Arial','FontSize',font_axis,'FontWeight','Bold');
axis equal;
hold off;
set(gcf,'color','w');
save_dir = fullfile(sorted_dir,'/Plots/Clustering/');
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end
save_name = fullfile(save_dir,[var_name,'_IE_Clustering.png']);
export_fig(save_name);
close all;
save_name = fullfile(sorted_dir,'clust_info.mat');
save(save_name,'clust_info');


