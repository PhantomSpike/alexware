function plot_histograms(sorted_dir,save_dir)
%This function plots the histograms for the different features
%% Define params
font_axis = 12; %Font size for plotting

if ~exist('save_dir','var')
    save_dir = fullfile(sorted_dir,'/Plots/Histogram_Plots');
    if ~exist(save_dir, 'dir')
        mkdir(save_dir);
    end
end

%% Load the clust_info
file_path = fullfile(sorted_dir,'clust_info.mat');
load(file_path);
%% Plot the histograms
row = 2;
col = 3;
no_bins = 40;
per = 0.01;
edgel = 0.047; edger = per; edgeh = 0.01; edgeb = 0.07; space_h = 0.02; space_v = 0.07;
[pos]=subplot_pos(row,col,edgel,edger,edgeh,edgeb,space_h,space_v);

figure('units','normalized','outerposition',[0 0 1 1]);
subplot('position',pos{1});
set(gca,'FontName','Arial','FontSize',font_axis);
histogram(clust_info.results.spike_width_ms,no_bins,'Normalization','probability');
xlabel('Peak-to-trough time [ms]','FontSize',font_axis,'fontweight','bold');
ylabel('Probability');
set(gca,'FontName','Arial','FontSize',font_axis,'fontweight','bold');

subplot('position',pos{2});
set(gca,'FontName','Arial','FontSize',font_axis);
histogram(clust_info.results.width_hh_trough_ms,no_bins,'Normalization','probability');
xlabel('Width @ ^{1}/_{2}height trough [ms]','FontSize',font_axis,'fontweight','bold');
set(gca,'FontName','Arial','FontSize',font_axis,'fontweight','bold');
set(gca,'ytick',[]);

subplot('position',pos{3});
histogram(clust_info.results.end_slope,no_bins,'Normalization','probability');
xlabel('End slope [uV/ms]','FontSize',font_axis,'fontweight','bold');
set(gca,'FontName','Arial','FontSize',font_axis,'fontweight','bold');
set(gca,'ytick',[]);

subplot('position',pos{4});
histogram(clust_info.results.peak_trough_ratio,no_bins,'Normalization','probability');
xlabel('Peak-to-trough ratio','FontSize',font_axis,'fontweight','bold');
ylabel('Probability');
set(gca,'FontName','Arial','FontSize',font_axis,'fontweight','bold');

subplot('position',pos{5});
histogram(clust_info.results.width_hh_peak_ms,no_bins,'Normalization','probability');
xlabel('Width @ ^{1}/_{2}height peak [ms]','FontSize',font_axis,'fontweight','bold');
set(gca,'FontName','Arial','FontSize',font_axis,'fontweight','bold');
set(gca,'ytick',[]);

set(gcf,'color','w');
save_name = fullfile(save_dir,'Histograms_spike_features.png');
export_fig(save_name);
set(gca,'ytick',[]);
close all;
