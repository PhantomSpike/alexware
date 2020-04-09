function psth = psth_func_noise2(Y,good_cluster_id,start_time_ms,first_bin_ms,last_bin_ms,t_bin_ms,noisepos,grid_root)
%psth_final = psth_func(Y,start_time_ms,grid_root)
%Add a normal description


%% Initial parameters

cond_ix = find(noisepos~=0);
num_cond = numel(cond_ix);
noisepos_ms = noisepos*1000;
edges_ms = (first_bin_ms:t_bin_ms:last_bin_ms); %The edges of the histogram

%% Get PSTHs
spike_times_ms = Y(:,1).*1000; % Get the spike times in ms
clusters = Y(:,2);
no_clusters = length(good_cluster_id);

dir_info = dir([grid_root '/*Info.mat']); %Get the names of all files with this extension
grid_filename = fullfile(dir_info.folder, dir_info.name); %Form the name of the file to be opened
grid_ben = load(grid_filename);

%Initialize the struct array
psth(no_clusters).data = [];
psth(no_clusters).cluster_id = [];

fprintf('== Generating PSTHs ==\n');tic;
parfor cluster = 1:no_clusters
    fprintf('== Processing cluster %.0f/%.0f ==\n',cluster,no_clusters);
    cluster_id = good_cluster_id(cluster); %Load the current unit

    for stim = 1:num_cond
        ix = cond_ix(stim);
        ix_stim_ms = start_time_ms(grid_ben.grid.randomisedGridSetIdx==ix); %Find the starting point for this stimulus
        ix_stim_ms = ix_stim_ms + noisepos_ms(ix);
        psth(cluster).data(stim,:) = histc(spike_times_ms(clusters == cluster_id),ix_stim_ms + edges_ms);
    end
    psth(cluster).cluster_id = cluster_id;
end
psth(1).params.t_bin_ms = t_bin_ms;
psth(1).params.edges_ms = edges_ms;
fprintf('== Done! Processing took %0.fs ==\n',toc);
end