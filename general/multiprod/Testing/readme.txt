Description of the clust_info file:

There are several useful variable that are:

*clust_id - Vector which has all the unique ids for the good and bad clusters as given by KS (Kilosort)
*clust_quality - Cell character array which has all the labels for the clusters i.e. good or bad. It has the same order as clust_id
*good_units - Matrix with three columns. First column is the cluster id. Second column is the index of cluster id. 
Third column is '0' for good units that were not used in the clsutering analysis because they were shit. The '1' are for good units that were used for clustering.
*clustering.idx - Matrix with two columns. The first column has '1' and '2' which tell you if the neuron was clustered as putative inhibitory or putative excitatory (Cfr clust_info.clustering.inh_flag).
The second column has the cluster id fo that neuron.
*clust_info.clustering.inh_flag = Either '1' or '2'. Tells you which flag means inhibitory neuron. Cfr clustering.idx
*npsp - Matrix with two columns. First column has the NPSP value for the reverb stimuli. The second column has the cluster id. The order is the same as in clust_id.
*spikeTimes - Matrix with two columns. The first column has every spike that ocurred from every cluster in seconds. The second one is the cluster id.