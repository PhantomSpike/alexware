function [psth_final,t_ms] = explore_psth(data_root,sorted_root)
plot_name = 'repeated-noise-6';
save_dir = '/home/phant0msp1ke/Desktop/PSTHs_Ben/';
% plot_name = [plot_name,'_',num2str(start_time),'s to ',num2str(end_time),'s'];

start_time_s = 2;
end_time_s = 252;
chunk_s = 50;
synch_ch = load_synch(data_root,start_time_s,end_time_s,chunk_s);

trig_min_length = 30000; %The minimum langth of one trigger in samples

qualityOpt = 'nonoise';
Y = get_spike_times(sorted_root,qualityOpt);


sweep_duration_ms = 1000; %The duration of one time interval of interest in ms. Actual Benware sweep duration is around 105.2ms
t_bin_ms = 5; %Time bin in ms
fs_s = 30000; %Sampling rate in samples/s
t_ms = [-200:t_bin_ms:sweep_duration_ms]; %The edges of the histogram

diff_sig = diff(synch_ch); %Find the difference between every n+1 sample - n sample. This tells us the the beginning/end of each sweep

%Find the sample no for the beginning of each sweep
start_ix = find(diff_sig==1); %Note that this diff will differ depending on which sync channel we use so always check it before analyzing
end_ix = find(diff_sig==-1); %Note that this diff will differ depending on which sync channel we use so always check it before analyzing
start_ix_ms = (start_ix/fs_s).*1000; %Convert the starting sample numbers to times in ms

if start_ix_ms(end) + sweep_duration_ms > end_time_s*1000
    start_ix(end) = [];
    start_ix_ms(end) = [];
end
diff_ix = end_ix - start_ix; %Find the length of each sweep in samples
start_ix = start_ix(diff_ix >= trig_min_length); %Keep only the triggers which have length >= minimum triger length


num_triggers = length(start_ix);

spike_times_ms = Y(:,1).*1000; % Get the spike times in ms
clusters = Y(:,2);
cluster_id = unique(clusters); %Sort the clusters which are good in ascending order
total_no_clusters = length(cluster_id); %The total number of unqiue clusters

%Find the psths for every cluster, stimulus and repetition and store
%them in a cell array

for cluster = 1:total_no_clusters
    current_cluster_id = cluster_id(cluster);
    fprintf('== Processing cluster %.0f/%.0f ==\n',cluster,total_no_clusters);
    for trigger = 1:num_triggers
        psth{cluster}(trigger,:) = histc(spike_times_ms(clusters == current_cluster_id),start_ix_ms(trigger) + t_ms);
    end
    psth{cluster} = psth{cluster}(:,1:end-1); %Delete the last bin which is weird
end
avg_psth = cellfun(@(x)(mean(x,1)),psth,'UniformOutput',false); %Find the average across repetitions
%Convert the cell array into 3D array with dimesnions stimuli x
%clusters x time bins
int_mat = cellfun(@(x)reshape(x,1,1,[]),avg_psth,'un',0);
psth_final(:,:) = cell2mat(int_mat);
psth_final = psth_final.*(1000/5);
%% Plotting 2
% channel_no = 304;
% figure('units','normalized','outerposition',[0 0 1 1]);
% bar(t_ms(1:end-1), psth_final(:,channel_no), 'hist');
% title(['Histogram of channel no ',num2str(channel_no)]);
% xlabel('Time [ms]');
% ylabel('Spikes/s');
% fig=gcf;
% set(findall(fig,'-property','FontSize'),'FontSize',20);
% save_name = [save_dir,plot_name,' channel ',num2str(channel_no),'.png'];
% export_fig(save_name);
% %% Plotting 3
% per = 0.01;
% edgel = 0.04; edger = per; edgeh = per; edgeb = 0.04; space_h = 0.02; space_v = 0.02;
% rows = 5;
% columns = 3;
% cluster_list = [0:15:90];
% [pos]=subplot_pos(rows,columns,edgel,edger,edgeh,edgeb,space_h,space_v);
% %     figure('units','normalized','outerposition',[0 0 1 1]);
% for chunk = 1:6
%     figure('units','normalized','outerposition',[0 0 1 1]);
%     cluster = cluster_list(chunk);
%     for ii = 1:15
%         subplot('position',pos{ii});
%         bar(t_ms(1:end-1), psth_final(cluster + ii,:), 'hist');
%         xlim([t_ms(1) t_ms(end-1)]);
%         legend(['Cluster ', num2str(cluster + ii)],'Location','northeast');
%         if ii == ((rows-1)*columns + 1)
%             xlabel('Time [ms]');
%             ylabel('Rate [spikes/s]')
%             set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold');
%         end
%     end
%     save_name = [save_dir,plot_name,' Chunk #',num2str(chunk),'.png'];
%     export_fig(save_name);
%     close(gcf);
% end

%% Plotting 3
per = 0.01;
edgel = 0.04; edger = per; edgeh = per; edgeb = 0.04; space_h = 0.02; space_v = 0.02;
rows = 6;
columns = 3;
cluster_list = [0:15:90];
[pos]=subplot_pos(rows,columns,edgel,edger,edgeh,edgeb,space_h,space_v);

for ii = 1:17
    subplot('position',pos{ii});
    bar(t_ms(1:end-1), psth_final(ii,:), 'hist');
    xlim([t_ms(1) t_ms(end-1)]);
    legend(['Cluster ', num2str(ii)],'Location','northeast');
    xlabel('Time [ms]');
    ylabel('Rate [spikes/s]')
    set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold');
end


