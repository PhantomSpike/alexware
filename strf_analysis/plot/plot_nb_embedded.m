function [psth_noise, psth_all_anech_avg, psth_all_reverb1_avg, psth_all_reverb2_avg] = plot_nb_embedded(pen,plot,store)

%Give the correct directory with Kilosorted files
data_root = ['//media/alex/4TB SSD/DATA/Ronnie_23_01_2018/bin_files/Preprocessed/CRA_done/',pen,'-reverb_with_noise_same_cra/',pen,'-reverb_with_noise_same_cra_sorted'];
%Give the correct Benware Grid folder
grid_root = ['/media/alex/4TB SSD/DATA/Ronnie_23_01_2018/benware_gridinfo/',pen,'-reverb_with_noise_same'];
sahani_path = ['/media/alex/4TB SSD/DATA/Ronnie_23_01_2018/metadata/',pen,'-reverb_with_noise_same_psth_sahani'];
synch_path = ['/media/alex/4TB SSD/DATA/Ronnie_23_01_2018/metadata/synch_ch_reverbwithnoisesame_',pen];

folder_name = [pen,'/noise_bursts/'];

NPSP_val = 40;
t_bin_ms = 10;
first_bin_ms = -100;
last_bin_ms = 600;
min_trig_length_s = 20;
fs = 30000;
qualityOpt = 'nonoise';
load('/media/alex/4TB SSD/DATA/Ronnie_23_01_2018/metadata/reverb_ordered_metadata_same');
%% Get the noise position in the reverb stim
load(sahani_path);
for stim_no = 1:numel(reverb_metadata)
    noisepos(stim_no) = reverb_metadata(stim_no).noisepos;
    noisepos = noisepos';
end

%% Get the synch channel
load(synch_path);

%% Get the triggers
if strcmp(pen,'P06')
    [start_time_ms] = get_triggers2(synch_ch,min_trig_length_s,fs);
else
    [start_time_ms] = get_triggers(synch_ch,min_trig_length_s,fs);
end

%% Get the spike times
Y = get_spike_times(data_root,qualityOpt);

%% Load the psth_sahani and sort clusters
load(sahani_path);
count = 0;
for cluster_no = 1:numel(psth_sahani)
    if psth_sahani(cluster_no).NPSP < NPSP_val
        count = count + 1;
        good_clusters(count,1) = psth_sahani(cluster_no).cluster_id;
        good_clusters(count,2) = psth_sahani(cluster_no).NPSP;
    end
end
good_clusters = sortrows(good_clusters,2);
good_cluster_id = good_clusters(:,1);
no_good_clusters = numel(good_cluster_id);

%% Get the psth
psth_noise = psth_func_noise2(Y,good_cluster_id,start_time_ms,first_bin_ms,last_bin_ms,t_bin_ms,noisepos,grid_root);

%% Plot the noise burts according to condition
if plot
    gen_name = '/home/phant0msp1ke/Desktop/Kernels/';
    clust_spacing = 19;
    num_groups = ceil(no_good_clusters/(clust_spacing+1));
    last_spacing = no_good_clusters - (num_groups-1)*(clust_spacing+1) - 1;
    first_clust_ix = [1:clust_spacing+1:(num_groups)*(clust_spacing+1)];
    anech_ix = [1:14];
    reverb1_ix = [15:28];
    reverb2_ix = [29:42];
    edges_ms = psth_noise(1).params.edges_ms;
    
    for group = 1:num_groups
        num_columns = 5;
        num_rows = 4;
        per = 0.02;
        edgel = 0.04; edger = per; edgeh = per; edgeb = 0.05; space_h = 0.01; space_v =0.01;
        [pos]=subplot_pos(num_rows,num_columns,edgel,edger,edgeh,edgeb,space_h,space_v);
        
        first_clust = first_clust_ix(group);
        last_clust = first_clust + clust_spacing;
        if group == num_groups
            last_clust = first_clust + last_spacing;
        end
        cluster_ix = [first_clust:last_clust];
        num_clust = numel(cluster_ix);
        if store
            figure('units','normalized','outerposition',[0 0 1 1]);
        else
            figure;
        end
        for cluster = 1:num_clust
            current_clust = cluster_ix(cluster);
            cluster_id = psth_noise(current_clust).cluster_id;
            psth_anech = psth_noise(current_clust).data(anech_ix,:);
            psth_anech = mean(psth_anech)*(1000/t_bin_ms);
            psth_reverb1 = psth_noise(current_clust).data(reverb1_ix,:);
            psth_reverb1 = mean(psth_reverb1)*(1000/t_bin_ms);
            psth_reverb2 = psth_noise(current_clust).data(reverb2_ix,:);
            psth_reverb2 = mean(psth_reverb2)*(1000/t_bin_ms);
            psth_all_anech(current_clust,:) = psth_anech;
            psth_all_reverb1(current_clust,:) = psth_reverb1;
            psth_all_reverb2(current_clust,:) = psth_reverb2;
            subplot('position',pos{cluster});
            stairs(edges_ms,psth_anech,'LineWidth',2);
            hold on;
            stairs(edges_ms,psth_reverb1,'k','LineWidth',2);
            stairs(edges_ms,psth_reverb2,'r','LineWidth',2);
            xlim([-50 300]);
            %     set(findall(gca, 'Type', 'Line'),'LineWidth',2);
            if cluster == num_columns*(num_rows - 1) + 1
                axis on;
                set(gca,'FontName','Arial','FontSize',10,'FontWeight','Bold');
                xlabel('Time [ms]','FontSize',14,'FontWeight','bold');
                ylabel('Firing rate [spikes/s]','FontSize',16,'FontWeight','bold');
                legend('Anech','Small','Big','Location','northeast');
            end
        end
        hold off;
        file_name = [gen_name,folder_name,'Clusters_',num2str(first_clust),'_',num2str(last_clust),'.png'];
        if store
            export_fig(file_name);
            close;
        end
    end
    psth_all_anech_avg = mean(psth_all_anech);
    psth_all_reverb1_avg = mean(psth_all_reverb1);
    psth_all_reverb2_avg = mean(psth_all_reverb2);
    if store
        figure('units','normalized','outerposition',[0 0 1 1]);
    else
        figure;
    end
    stairs(edges_ms,psth_all_anech_avg,'LineWidth',2);
    hold on;
    stairs(edges_ms,psth_all_reverb1_avg,'k','LineWidth',2);
    stairs(edges_ms,psth_all_reverb2_avg,'r','LineWidth',2);
    set(gca,'FontName','Arial','FontSize',10,'FontWeight','Bold');
    xlabel('Time [ms]','FontSize',27,'FontWeight','bold');
    ylabel('Firing rate [spikes/s]','FontSize',27,'FontWeight','bold');
    title(['Embeddd noise burst response over all neurons for ',pen],'FontSize',27,'FontWeight','Bold');
    legend('Anech','Small','Big','Location','northeast');
    file_name2 = [gen_name,folder_name,pen,'_Avg_resp_all_neurons','.png'];
    if store
        export_fig(file_name2);
        close;
    end
end
end
