function plot_corr(sorted_dir,save_dir)
%This function plots the histograms for the different features
%% Define params
font_axis = 10; %Font size for plotting
font_axis_y = 10;
num_bins = 40;
mrk_size = 4;
if ~exist('save_dir','var')
    save_dir = fullfile(sorted_dir,'/Plots/Histogram_Plots');
    if ~exist(save_dir, 'dir')
        mkdir(save_dir);
    end
end

%% Load the clust_info
file_path = fullfile(sorted_dir,'clust_info.mat');
load(file_path);
%% Get and plot correlation matrix
X = [clust_info.results.spike_width_ms,clust_info.results.width_hh_trough_ms,clust_info.results.end_slope,clust_info.results.peak_trough_ratio,clust_info.results.width_hh_peak_ms];
norm_X = (X-mean(X))./std(X,[],1); %Z-score X to normalize for the different nature of the features 
[C,P] = corrcoef(X); %Compute Correlation matrix

%Plot the correlation matrix
figure('units','normalized','outerposition',[0 0 1 1]);
imagesc(C);
colormap('jet');
colorbar;
caxis([-1 1])
labelNames = {'Peak-to-trough width','Width 1/2 h trough','End slope','Peak-to-trough ratio','Width 1/2 h peak'};
n = length(labelNames);
set(gca, 'XTick', 1:n); % center x-axis ticks on bins
set(gca, 'YTick', 1:n); % center y-axis ticks on bins
set(gca,'XTickLabel',labelNames,'FontSize',font_axis,'fontweight','bold');   % gca gets the current axis
set(gca,'YTickLabel',labelNames,'FontSize',font_axis,'fontweight','bold');   % gca gets the current axis
save_name = fullfile(save_dir,['Correlogram.png']);
export_fig(save_name);
close;

%Plot histogram and scatter plots with correlations for all features for raw data
figure('units','normalized','outerposition',[0 0 1 1]);
[R,PValue,H1] = corrplot(X,'varNames',{'P-t-t w','W1/2ht','Endslope','P-t-t r','W1/2h p'},'type','Spearman','testR','on');
for h = 1:size(H1,2)
    H1(h,h).NumBins = num_bins;
end
save_name = fullfile(save_dir,['Correlogram_scatter.png']);
export_fig(save_name);
close;

%Plot histogram and scatter plots for all features for raw data
figure('units','normalized','outerposition',[0 0 1 1]);
[S,ax,BigAx,H2,HAx] = plotmatrix(X);

for h = 1:size(H2,2)
    H2(h).NumBins = num_bins;
end

for s = 1:numel(S)
    S(s).MarkerSize = mrk_size;
end

ax(1,1).YLabel.String='Peak-trough width'; ax(1,1).FontWeight = 'bold';ax(1,1).FontSize = font_axis_y;
ax(2,1).YLabel.String='Width ^{1}/_{2}h trough'; ax(2,1).FontWeight = 'bold';ax(2,1).FontSize = font_axis_y;
ax(3,1).YLabel.String='End slope'; ax(3,1).FontWeight = 'bold';ax(3,1).FontSize = font_axis_y;
ax(4,1).YLabel.String='^{Peak}/_{Trough}';ax(4,1).FontWeight = 'bold';ax(4,1).FontSize = font_axis_y;
ax(5,1).YLabel.String='Width ^{1}/_{2}h peak';ax(5,1).FontWeight = 'bold';ax(5,1).YLabel.FontSize = font_axis_y;
ax(5,1).XLabel.String='Peak-trough width'; ax(5,1).FontWeight = 'bold';ax(5,1).XLabel.FontSize = font_axis;
ax(5,2).XLabel.String='Width ^{1}/_{2}h trough'; ax(5,2).FontWeight = 'bold';ax(5,2).FontSize = font_axis;
ax(5,3).XLabel.String='End slope'; ax(5,3).FontWeight = 'bold';ax(5,3).FontSize = font_axis;
ax(5,4).XLabel.String='^{Peak}/_{Trough}';ax(5,4).FontWeight = 'bold';ax(5,4).FontSize = font_axis;
ax(5,5).XLabel.String='Width ^{1}/_{2} height peak';ax(5,5).FontWeight = 'bold';ax(5,5).FontSize = font_axis;

save_name = fullfile(save_dir,['Scatter_plots.png']);
export_fig(save_name);
close;

%Plot histogram and scatter plots for all features for z-scored data
figure('units','normalized','outerposition',[0 0 1 1]);
[S,ax,BigAx,H3,HAx] = plotmatrix(norm_X);
for h = 1:size(H3,2)
    H3(h).NumBins = num_bins;
end

for s = 1:numel(S)
    S(s).MarkerSize = mrk_size;
end

ax(1,1).YLabel.String='Peak-trough width'; ax(1,1).FontWeight = 'bold';ax(1,1).FontSize = font_axis_y;
ax(2,1).YLabel.String='Width ^{1}/_{2}h trough'; ax(2,1).FontWeight = 'bold';ax(2,1).FontSize = font_axis_y;
ax(3,1).YLabel.String='End slope'; ax(3,1).FontWeight = 'bold';ax(3,1).FontSize = font_axis_y;
ax(4,1).YLabel.String='^{Peak}/_{Trough}';ax(4,1).FontWeight = 'bold';ax(4,1).FontSize = font_axis_y;
ax(5,1).YLabel.String='Width ^{1}/_{2}h peak';ax(5,1).FontWeight = 'bold';ax(5,1).YLabel.FontSize = font_axis_y;
ax(5,1).XLabel.String='Peak-trough width'; ax(5,1).FontWeight = 'bold';ax(5,1).XLabel.FontSize = font_axis;
ax(5,2).XLabel.String='Width ^{1}/_{2}h trough'; ax(5,2).FontWeight = 'bold';ax(5,2).FontSize = font_axis;
ax(5,3).XLabel.String='End slope'; ax(5,3).FontWeight = 'bold';ax(5,3).FontSize = font_axis;
ax(5,4).XLabel.String='^{Peak}/_{Trough}';ax(5,4).FontWeight = 'bold';ax(5,4).FontSize = font_axis;
ax(5,5).XLabel.String='Width ^{1}/_{2} height peak';ax(5,5).FontWeight = 'bold';ax(5,5).FontSize = font_axis;
set(gca,'FontName','Arial','FontSize',font_axis,'FontWeight','Bold');
save_name = fullfile(save_dir,['Scatter_plots_normalized.png']);
export_fig(save_name);
close;
