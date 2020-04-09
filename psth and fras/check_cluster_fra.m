synch_path = '/media/alex/Alex SSD/Kilosort_optimize/metadata/P02_quning_synch.mat';
grid_root = '/media/alex/Alex SSD/Derry/benware_gridinfo/P02-quning';
sorted_root = '/media/alex/Alex SSD/Kilosort_optimize/Sorted/1_Normal/P02-quning/CRA';
load([sorted_root,'/clust_info']);
save_dir = '/media/alex/Alex SSD/Kilosort_optimize/Figures/';
%% Generate the fras
fra_psth = fra_pixels2(synch_path,grid_root,sorted_root);
%% Find clusters with p-val from the anovan test lower than a certain threshold
p_t = 0.05;
pval = sortrows(fra_psth(1).pval,1);
ix_t = pval(pval(:,1) < p_t,3);
clusters = fra_psth(1).cluster_id(ix_t); 

%% Plot the FRAs
plot_fra_pixels(fra_psth,clusters,save_dir);

%% Plot the corresponding spike clusters
plot_clusters(clust_info,clusters,save_dir);