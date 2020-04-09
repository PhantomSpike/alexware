function [clust_info] = run_kmeans(clust_info,features,plot_on)
%[clust_info] = run_kmeans(clust_info,plot_on)
%This function run k-means on the spike features for every cluster
%extracted with the featire_extract function.
%>>INPUT>>
%clust_info - clust_info var that cotnains all data and metadata
%features - vector with 2 values which determines which features to use for
%the cluster
%1 - Spike width (peak-to-trough width), this is the distance form peak (depolarization) to
%trough (hyperpolarization) in the wf, var: spike_width_ms
%2 - Width at 1/2 height trough, this is the width of the trough part at 1/2
%max value, var: width_hh_trough_ms
%3 - End slope, this is the slope of the decaying part of the trough, var: end_slope
%4 - Peak-to-trough ratio, this is the ratio of the amplitudes of the peak
%to the trough PA/TA, var: peak_trough_ratio
%5 - Width at 1/2 height peak, this is the width of the peak part at 1/2
%max value, var: width_hh_peak_ms

switch features(1)
    case 1
        x_feature = clust_info.results.spike_width_ms;
        x_lab = ('Spike width [ms]');
    case 2
        x_feature = clust_info.results.width_hh_trough_ms;
        x_lab = ('Width at ^{1}/_{2}height trough [ms]');
    case 3
        x_feature = clust_info.results.end_slope;
        x_lab = ('End slope [^{uV}/_{ms}]');
    case 4
        x_feature = clust_info.results.peak_trough_ratio;
        x_lab = ('^{Peak}/_{Trough} [AU]');
    case 5
        x_feature = clust_info.results.width_hh_peak_ms;
        x_lab = ('Width at ^{1}/_{2}height peak [ms]');
end

switch features(2)
    case 1
        y_feature = clust_info.results.spike_width_ms;
        y_lab = ('Spike width [ms]');
    case 2
        y_feature = clust_info.results.width_hh_trough_ms;
        y_lab = ('Width at ^{1}/_{2}height trough [ms]');
    case 3
        y_feature = clust_info.results.end_slope;
        y_lab = ('End slope [^{uV}/_{ms}]');
    case 4
        y_feature = clust_info.results.peak_trough_ratio;
        y_lab = ('^{Peak}/_{Trough} [AU]');
    case 5
        y_feature = clust_info.results.width_hh_peak_ms;
        y_lab = ('Width at ^{1}/_{2}height peak [ms]');
end

font_axis = 20;

if ~exist('plot_on','var')
    plot_on = false;
end

X = [x_feature,y_feature];
mean_X = mean(X);
std_X = std(X);
X_norm = (X - mean_X)./std_X;
[idx,C_norm] = kmeans(X_norm,2,'Replicates',100,'Distance','sqeuclidean');
C = (C_norm.*std_X) + mean_X;
[~,ix_min_width] = min(X(:,1));
inh_flag = idx(ix_min_width);

clust_info.clustering.X = X;
clust_info.clustering.X_norm = X_norm;
clust_info.clustering.C = C;
clust_info.clustering.C_norm = C_norm;
clust_info.clustering.X = X;
clust_info.clustering.idx(:,1) = clust_info.results.normal_cluster_id;
idx_copy = idx; %Here we want to assign ix to exc and inh and always the same
idx(idx_copy~=inh_flag) = 1; %These are excitatory
idx(idx_copy==inh_flag) = 2; %These are inhibitory
inh_flag = 2; %Inhibitory neurons are always 2
clust_info.clustering.idx(:,2) = idx;
clust_info.clustering.inh_flag = 2;
clust_info.clustering.x_lab = x_lab;
clust_info.clustering.y_lab = y_lab;
ix1 = find(x_lab == '[', 1, 'last');
ix2 = find(y_lab == '[', 1, 'last');
x_lab_norm = [x_lab(1:ix1-1),' [Z-score]'];
y_lab_norm = [y_lab(1:ix2-1),' [Z-score]'];
clust_info.clustering.x_lab_norm = x_lab_norm;
clust_info.clustering.y_lab_norm = y_lab_norm;


if plot_on
    figure('units','normalized','outerposition',[0 0 1 1]);
    plot(X(idx==inh_flag,1),X(idx==inh_flag,2),'b.','MarkerSize',20)
    hold on;
    plot(X(idx~=inh_flag,1),X(idx~=inh_flag,2),'r.','MarkerSize',20)
    plot(C(:,1),C(:,2),'kx',...
        'MarkerSize',15,'LineWidth',3)
    legend('Putative inhibitory','Putative excitatory','Centroids',...
        'Location','NE')
    title('Inhibitory/Excitatory cluster assignment and centroids');
    xlabel(x_lab);
    ylabel(y_lab);
    set(gca,'FontName','Arial','FontSize',font_axis,'FontWeight','Bold');
%     axis equal;
    hold off;
    set(gcf,'color','w');
    
    %Normalized plot (Z-scored)
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
end

