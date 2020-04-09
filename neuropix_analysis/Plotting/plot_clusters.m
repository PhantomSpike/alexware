function plot_clusters(clust_info,clusters,save_dir)

num_clust = numel(clusters);

row = 5;
col = 4;
per = 0.005;
edgel = 0.035; edger = per; edgeh = per; edgeb = 0.05; space_h = per; space_v = 0.01;
[pos]=subplot_pos(row,col,edgel,edger,edgeh,edgeb,space_h,space_v);

%% Plot waveforms
clust_spacing = row*col - 1;
num_groups = ceil(num_clust/(clust_spacing+1));
last_spacing = num_clust - (num_groups-1)*(clust_spacing+1) - 1;
first_clust_ix = [1:clust_spacing+1:(num_groups)*(clust_spacing+1)];
for group = 1:num_groups
    
    first_clust = first_clust_ix(group);
    last_clust = first_clust + clust_spacing;
    if group == num_groups
        last_clust = first_clust + last_spacing;
    end
    num_subplots = numel(first_clust:last_clust);
    plot_ix = [first_clust:last_clust];
    min_val = min(min(clust_info.mean_wfs(plot_ix,:),[],2));
    max_val = max(max(clust_info.mean_wfs(plot_ix,:),[],2));
    figure('units','normalized','outerposition',[0 0 1 1]);
    for ii = 1:num_subplots
        subplot('position',pos{ii});
        ix_clust = plot_ix(ii);
        clust = clusters(ix_clust);
        ix = find(clust_info.clust_id == clust);
        shadedErrorBar(clust_info.wf_time,clust_info.mean_wfs(ix,:),clust_info.std_wfs(ix,:),'b');
        legend([num2str(clust), ' ',clust_info.clust_quality{ix,1},' ',num2str(clust_info.clust_depth(ix),'%.0f'),'um ',num2str(clust_info.no_spikes(ix)), 'sp ',num2str(clust_info.clust_spikedur(ix),'%.2f'),'ms'])
        if ii == (row - 1)*col + 1
            xlabel('Time [ms]');
        elseif ii > (row - 1)*col + 1 && ii <= row*col
            xlabel('Time [ms]');
            set(gca,'ytick',[]);
        else
            set(gca,'ytick',[]);
            set(gca,'xtick',[]);
        end
        set(gca,'FontName','Arial','FontSize',10,'FontWeight','Bold');
        xlim([clust_info.wf_time(1) clust_info.wf_time(end)]);
        ylim([min_val max_val]);
    end
    if exist('save_dir','var')
        file_name = ['Spike_wfs_Clusters_',num2str(clusters(plot_ix(1))),'_to_',num2str(clusters(plot_ix(end)))];
        save_name = fullfile(save_dir,[file_name,'.png']);
        export_fig(save_name);
        close;
    end
end
