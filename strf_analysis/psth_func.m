function psth_final = psth_func(Y,start_time_ms,psth_end_time_s,t_bin_ms,grid_root)
%psth_final = psth_func(Y,start_time_ms,grid_root)
%Add a normal description


%% Initial parameters
no_conds = 6; %The number of sound stimuli
no_repeats = 10; %The number of repeats for each echoic condition
num_stim = no_conds*no_repeats; %Find the total number of stimuli
cond_ix = (1:no_repeats:num_stim); %Generate indices to extract the relevant conditions

psth_end_time_ms = psth_end_time_s*1000; %Convert to ms 
edges_ms = (0:t_bin_ms:psth_end_time_ms); %The edges of the histogram

%% Get PSTHs
spike_times_ms = Y(:,1).*1000; % Get the spike times in ms
clusters = Y(:,2);
ix_clusters = unique(clusters); %Find all the units in the input file
no_clusters = length(ix_clusters);

dir_info = dir([grid_root '/*Info.mat']); %Get the names of all files with this extension
grid_filename = fullfile(dir_info.folder, dir_info.name); %Form the name of the file to be opened
grid_ben = load(grid_filename);

%Initialize the struct array
psth_final(no_clusters).data = [];
psth_final(no_clusters).cluster_id = [];

fprintf('== Generating PSTHs ==\n');tic;
parfor cluster = 1:no_clusters
    fprintf('== Processing cluster %.0f/%.0f ==\n',cluster,no_clusters);
    cluster_id = ix_clusters(cluster); %Load the current unit
    psth = zeros(num_stim,numel(edges_ms));
    for stim = 1:num_stim
        ix_stim_ms = start_time_ms(grid_ben.grid.randomisedGridSetIdx==stim); %Find the starting point for this stimulus
        psth(stim,:) = histc(spike_times_ms(clusters == cluster_id),ix_stim_ms + edges_ms);
    end
    for stim_cond = 1:numel(cond_ix)
        if stim_cond == 1
            psth_final(cluster).data = psth((cond_ix(stim_cond):cond_ix(stim_cond) + no_repeats -1),:);
        else
            psth_final(cluster).data = [psth_final(cluster).data,psth(cond_ix(stim_cond):(cond_ix(stim_cond) + no_repeats-1),:)];
        end
    end
    psth_final(cluster).cluster_id = cluster_id;
end

fprintf('== Done! Processing took %0.fs ==\n',toc);
end