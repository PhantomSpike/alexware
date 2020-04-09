%% Find good clusters
NPSP_threshold = 60;
for ii = 1:numel(psth)
ix(ii) = psth(ii).NPSP < NPSP_threshold;
end
good_clusters = find(ix ~= 0);
%% Plot them
num_rows = numel(psth(1).params.dB_levels);
num_columns = numel(psth(1).params.freqs);
num_stim = num_rows*num_columns;
index = reshape(1:num_stim, num_columns, num_rows).';
index = index(:);
scaling_factor = (1000/psth(1).params.t_bin_ms);
for cluster = 1:numel(good_clusters)
    figure('units','normalized','outerposition',[0 0 1 1]);
    cluster_ix = good_clusters(cluster);
    psth_plot = psth(cluster_ix).X_dbft;
    psth_plot = reshape(psth_plot,num_rows*num_columns,numel(psth(1).params.t_ms));
    psth_plot = scaling_factor*psth_plot;
    per = 0.02;
    edgel = 0.03; edger = per; edgeh = per; edgeb = 0.05; space_h = 0.01; space_v =0.01;
    [pos]=subplot_pos(num_rows,num_columns,edgel,edger,edgeh,edgeb,space_h,space_v);
    for stim = 1:num_stim
        ix = index(stim);
        subplot('position',pos{ix});
        bar(psth(1).params.t_ms, psth_plot(stim,:), 'hist');
        xlim([psth(1).params.t_ms(1) psth(1).params.t_ms(end)]);
        psth_max = max(psth_plot(:));
        ylim([0 psth_max]);
        start_line = line(zeros(2,1),[0,psth_max],'Color','r');
        end_line = line(100*ones(2,1),[0,psth_max],'Color','r');
        axis off;
        if stim == num_rows
            axis on;
            xlabel('Time [ms]');
            ylabel('Spike rate [Hz]');
            legend(['Cluster ID = ',num2str(psth(cluster_ix).cluster_id)]);
            set(gca,'FontName','Arial','FontSize',10,'FontWeight','Bold');
        end
    end
    drawnow;
    pause;
    close all;
end