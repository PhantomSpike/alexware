function clust_times = gets_clust(spike_times,clust_id)

clust_ix = spike_times(:,2) == clust_id;
clust_times = spike_times(clust_ix,1);