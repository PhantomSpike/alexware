%This scirpt will load the clust_info, extract spike times for every
%cluster and then compute the CCG between all good clusters and every other
%cluster
sorted_dir = '/mnt/40086D4C086D41D0/Inhibitory_excitatory_analysis/For_analysis/Ronnie/P08/P08-quning.2';
%% Params
t_bin_ms = 1; %The binning in ms
tau = 50; %How many shifts tau to do in the future and past
test_start_ms = 1; %Where to start looking for increase/decrease in firing in the CCG
test_end_ms = 5; %Where to stop looking for increase/decrease in firing in the CCG
base_start_ms = 25; %Where to start for the baseline window computation
base_end_ms = 50; %Where to finnish for the baseline window computation
z_th = 3; %How many std from the mean has the increase/decrease in firing rate to be sginificant 
sig_bins = 2; %How many bins need to be significant to count the neuron
%% Load data
file_path = fullfile(sorted_dir,'clust_info.mat');
load(file_path);
spike_times = clust_info.spikeTimes;
%% Set things
t_shift = [-tau:tau];
t_end_s = max(spike_times(:,1)) + 0.1;
t_bin_s = t_bin_ms/1000;
t_s = [0:t_bin_s:t_end_s];
num_t_points = length(t_s)-1;
num_shifts = length(t_shift);
zero_bin = ceil(num_shifts/2);
test_start_bin = zero_bin + ceil(test_start_ms/t_bin_ms);
test_end_bin = zero_bin + ceil(test_end_ms/t_bin_ms);
base_start_bin_pos = zero_bin + ceil(base_start_ms/t_bin_ms);
base_end_bin_pos = zero_bin + ceil(base_end_ms/t_bin_ms);
base_start_bin_neg = zero_bin - ceil(base_start_ms/t_bin_ms);
base_end_bin_neg = zero_bin - ceil(base_end_ms/t_bin_ms);
%% Get spike times for every cluster and store in matrix
clust_id = clust_info.clustering.idx(:,1); %Get the cluster id of only the good clusters that where used for I/E classification
num_clust = length(clust_id); %Find the total number of clusters
%Initialize matrices that will store the data
X_gt = zeros(num_clust,num_t_points); %g - all good clsuters, t - time points for histogram
X_ggs = zeros(num_clust,num_clust,num_shifts); %g - all good clsuters, t - time points for histogram

for c = 1:num_clust
    curr_id = clust_id(c);
    clust_times = gets_clust(spike_times,curr_id);
    spike_vec = histcounts(clust_times,t_s);
    spike_vec = spike_vec - mean(spike_vec); %Remove the mean for every cluster
    X_gt(c,:) = spike_vec;
end
%% Compute Correlograms for every good normal neuron with every other neuron

tic;
for s = 1:num_shifts
    fprintf('== Processing shift %.0f/%0.f ==\n',s,num_shifts);
    tau_current = t_shift(s);
    X_gt_shift = circshift(X_gt,tau_current,2); %Shift the response matrix X_ct with value = tau along time (second  dimension,rows)
    X_ggs(:,:,s) = (X_gt*X_gt_shift'); %Compute the cross-correlogram
end
X_ggs = X_ggs/num_t_points;
fprintf('== Done! Processing took %0.fs ==\n',toc);

%% Compute significance using Z-scores
base_X_ggs = cat(3,X_ggs(:,:,base_start_bin_neg:base_end_bin_neg),X_ggs(:,:,base_start_bin_pos:base_end_bin_pos)); %Get the data necessary for the base period 
test_X_ggs = X_ggs(:,:,test_start_bin:test_end_bin); %Get the data necessary for the test period 
mean_base_X_gg = mean(base_X_ggs,3); %Compute the mean of the base period for every pair
std_base_X_gg = std(base_X_ggs,[],3); %Compute the std of the base period for every pair
z_score_X_ggs = (test_X_ggs - mean_base_X_gg)./std_base_X_gg; %Z-score the data

for c = 1:num_clust
    z_curr_clust{c} = squeeze(z_score_X_ggs(c,:,:));
    z_curr_clust{c}(c,:) = 0; %Remove the result of the neuron with itself
    ix_exc = find(z_curr_clust{c}>=z_th); %Find the CCG values that are significantly bigger than the baseline
    ix_inh = find(z_curr_clust{c}<=-z_th); %Find the CCG values that are significantly smaller than the baseline
    %Find which neurons excite/inhibit which other neurons
    [clust_ix_excited,t_ix_exc] = ind2sub(size(z_curr_clust{c}),ix_exc); %Convert lin ix to 2D
    [clust_ix_inhibited,t_ix_inh] = ind2sub(size(z_curr_clust{c}),ix_inh);
    
    excited_ix = unique(clust_ix_excited); % Find all the neurons that have exc interaction with this neuron
    N_exc = histc(clust_ix_excited, excited_ix); %Find the ones which have at least x significant bins
    excited_ix = excited_ix(N_exc>=sig_bins);
    clust_id_excited = clust_id(excited_ix);
    
    inhibited_ix = unique(clust_ix_inhibited); % Find all the neurons that have interaction
    N_inh = histc(clust_ix_inhibited, inhibited_ix); %
    inhibited_ix = inhibited_ix(N_inh>=sig_bins);
    clust_id_inhibited = clust_id(inhibited_ix);
    
    %Classify the neurons based on the results
    if ~isempty(clust_id_excited) || ~isempty(clust_id_inhibited)
       if isempty(clust_id_inhibited) || length(clust_id_excited)>length(clust_id_inhibited)
           ei_label(c) = 1;
           id_excited{c} = clust_id_excited;
           id_inhibited{c} = [];
       elseif length(clust_id_excited)==length(clust_id_inhibited)
           ei_label(c) = 3;
           id_excited{c} = clust_id_excited;
           id_inhibited{c} = clust_id_inhibited;
       else
           ei_label(c) = 2;
           id_excited{c} = [];
           id_inhibited{c} = clust_id_inhibited;
       end
    else
        ei_label(c) = 0;
        id_excited{c} = [];
        id_inhibited{c} = [];
    end
    
end
clust_info.clustering.idx(:,3) = ei_label;
clust_info.clustering.id_excited = id_excited;
clust_info.clustering.id_inhibited = id_inhibited;
clust_info.clustering.xcorr_ccs = X_ggs;

clust_info.clustering.neuron_labels(1).type = 'No interaction';
clust_info.clustering.neuron_labels(2).type = 'Excitatory';
clust_info.clustering.neuron_labels(3).type = 'Inhibitory';
clust_info.clustering.neuron_labels(4).type = 'Both';
clust_info.clustering.neuron_labels(1).number = 0;
clust_info.clustering.neuron_labels(2).number = 1;
clust_info.clustering.neuron_labels(3).number = 2;
clust_info.clustering.neuron_labels(4).number = 3;

save_name = fullfile(sorted_dir,'clust_info.mat');
save(save_name,'clust_info');