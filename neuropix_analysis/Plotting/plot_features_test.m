function plot_features_test(sorted_dir,save_dir)
%This function plots the raw waveforms together with the extracted features
%so you can visually inspect if there is an error somewhere

%% Define params
skip = 3; %How many clusters to skip for the plottign i.e. every n-th cluster will be plotted
font_axis = 20; %Font size for plotting
lw = 3; %Specify linewidth


if ~exist('save_dir','var')
    save_dir = fullfile(sorted_dir,'/Plots/Test_Plots');
    if ~exist(save_dir, 'dir')
        mkdir(save_dir);
    end
end

%% Load the clust_info
file_path = fullfile(sorted_dir,'clust_info.mat');
load(file_path);
%% Test plotting
num_good_clust = size(clust_info.results.normal_int_wf,1);

for cl = 1:skip:num_good_clust
    figure('units','normalized','outerposition',[0 0 1 1]);
    clust_id = clust_info.results.normal_cluster_id(cl,1);
    t = clust_info.int_time_ms;
    wf = clust_info.results.normal_int_wf(cl,:);
    sp_width_ms = clust_info.results.spike_width_ms(cl);
    width_hh_peak_ms = clust_info.results.width_hh_peak_ms(cl);
    width_hh_trough_ms = clust_info.results.width_hh_trough_ms(cl);
    slope = clust_info.results.end_slope(cl);
    p_t_ratio = clust_info.results.peak_trough_ratio(cl);
    ix_peak = clust_info.results.ix_peak(cl);
    ix_trough = clust_info.results.ix_trough(cl);
    left_p_ix = clust_info.results.left_bound_ix_p(cl);
    right_p_ix = clust_info.results.right_bound_ix_p(cl);
    left_t_ix = clust_info.results.left_bound_ix_t(cl);
    right_t_ix = clust_info.results.right_bound_ix_t(cl);
    m = clust_info.results.end_slope(cl);
    c = clust_info.results.y_int(cl);
    slope_st = clust_info.results.ix_endslope_start(cl);
    slope_end = clust_info.results.ix_endslope_end(cl);
    y_hat = m*t(slope_st:slope_end) + c;
    plot(t,wf,'LineWidth',lw);
    hold on;
    xlabel('Time [ms]');
    ylabel('Amplitude [uV]');
    title(['Check plot cluster # ',num2str(clust_id)]);
    plot(t(ix_peak:ix_trough),ones(length(ix_peak:ix_trough),1)*wf(ix_peak),'k-','LineWidth',lw); %Plot the  line showing spike width (peak-to-trough distance) in black
    plot(t(left_p_ix:right_p_ix),ones(length(left_p_ix:right_p_ix),1)*wf(left_p_ix),'r-','LineWidth',lw); %Plot the  line showing peak width at 1/2 height in red
    plot(t(left_t_ix:right_t_ix),ones(length(left_t_ix:right_t_ix),1)*wf(left_t_ix),'g-','LineWidth',lw); %Plot the  line showing peak width at 1/2 height in green
    plot(t(slope_st:slope_end),y_hat,'m-','LineWidth',lw); %Plot the predicted end-slope line in magenta
    xline(t(ix_trough),'k--','LineWidth',lw-0.5); %Plot the vertical line from the trough down
    legend('Mean wf',['p-t dist ',num2str(sp_width_ms,2),' ms'],...
        ['W @ 1/2 h peak ',num2str(width_hh_peak_ms,2),' ms'],...
        ['W @ 1/2 h trough ',num2str(width_hh_trough_ms,2),' ms'],...
        ['End slope ',num2str(slope,2),' uV/ms'],...
        ['p-t ratio ',num2str(p_t_ratio,2)],'Location','southwest');
    set(gca,'FontName','Arial','FontSize',font_axis,'FontWeight','Bold');
    set(gcf,'color','w');
    save_name = fullfile(save_dir,['Cluster_',num2str(clust_id),'.png']);
    export_fig(save_name);
    close;
end