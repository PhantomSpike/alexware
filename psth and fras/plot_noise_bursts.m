%% Find good clusters
clear good_clusters;
count = 0;
for ii = 1:numel(psth_sahani)
    if psth_sahani(ii).NPSP < 40
        count = count + 1;
        good_clusters(count,1) = ii;
        good_clusters(count,2) = psth_sahani(ii).NPSP;
    end
end
good_clusters = sortrows(good_clusters,2);
%% Plot the results
row = 5;
col = 10;
step = 10;
per = 0.01;
edgel = 0.03; edger = per; edgeh = per; edgeb = 0.05; space_h = 0.005; space_v =0.005;
[pos]=subplot_pos(row,col,edgel,edger,edgeh,edgeb,space_h,space_v);
for ii = 1:row*col
    subplot('position',pos{ii});
    cluster = good_clusters(ii,1);
    imagesc(psth_sahani(cluster).data);
    axis off;
    if ii == col*(row - 1) + 1
        axis on;
        time_ms = psth_sahani(1).params.edges_ms([1:step:numel(psth_sahani(1).params.edges_ms)-1]);
        xticks([1:step:numel(psth_sahani(1).params.edges_ms)-1]);
        for jj = 1:numel(time_ms)
            x_labels{jj} = num2str(time_ms(jj),'%.0f');
        end
        xticklabels(x_labels);
        set(gca,'FontName','Arial','FontSize',10,'FontWeight','Bold');
        xlabel('Time [ms]','FontSize',16,'FontWeight','bold');
        ylabel('Trial number','FontSize',16,'FontWeight','bold');
    end
end

%% Plot the results
per = 0.02;
edgel = 0.03; edger = per; edgeh = per; edgeb = 0.05; space_h = 0.01; space_v =0.01;
[pos]=subplot_pos(10,8,edgel,edger,edgeh,edgeb,space_h,space_v);
all_psth = mean([psth_sahani.data]);
max_val = max(all_psth(:));
for ii = 1:80
    subplot('position',pos{ii});
    cluster = good_clusters(ii);
    plot(mean(psth_sahani(cluster).data));
    ylim([0 max_val]);
end

%% Plot the results for all
ix = 139:238;
per = 0.02;
edgel = 0.03; edger = per; edgeh = per; edgeb = 0.05; space_h = 0.01; space_v =0.01;
[pos]=subplot_pos(10,10,edgel,edger,edgeh,edgeb,space_h,space_v);
for ii = 1:100
    subplot('position',pos{ii});
    cluster = ix(ii);
    imagesc(psth_sahani(cluster).data);
end