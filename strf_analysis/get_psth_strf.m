function [psth_strf] = get_psth_strf(Y,start_time_ms,t_s,grid_root,good_cluster_id,n_h)
%psth_final = psth_func(Y,start_time_ms,grid_root)
%Add a normal description


%% Initial parameters
no_conds = 6; %The number of sound stimuli
no_repeats = 10; %The number of repeats for each echoic condition
num_stim = no_conds*no_repeats; %Find the total number of stimuli
cond_ix = (1:no_repeats:num_stim); %Generate indices to extract the relevant conditions

edges_ms = 1000*t_s; %The edges of the histogram converted from s to ms

%% Get PSTHs
spike_times_ms = Y(:,1).*1000; % Get the spike times in ms
clusters = Y(:,2);
no_clusters = length(good_cluster_id);

dir_info = dir([grid_root '/*Info.mat']); %Get the names of all files with this extension
grid_filename = fullfile(dir_info.folder, dir_info.name); %Form the name of the file to be opened
grid_ben = load(grid_filename);

%Initialize the struct array
psth_strf(no_clusters).anech = [];
psth_strf(no_clusters).reverb1 = [];
psth_strf(no_clusters).reverb2 = [];
psth_strf(no_clusters).cluster_id = [];

fprintf('== Generating PSTHs ==\n');tic;
parfor cluster = 1:no_clusters
    fprintf('== Processing cluster %.0f/%.0f ==\n',cluster,no_clusters);
    cluster_id = good_cluster_id(cluster); %Load the current unit
    psth = zeros(num_stim,numel(edges_ms));
    for stim = 1:num_stim
        ix_stim_ms = start_time_ms(grid_ben.grid.randomisedGridSetIdx==stim); %Find the starting point for this stimulus
        psth(stim,:) = histc(spike_times_ms(clusters == cluster_id),ix_stim_ms + edges_ms);
    end
    
    psth = psth(:,n_h:end);
    
    psth_strf(cluster).anech = [psth(1:10,:),psth(11:20,:)];
    psth_strf(cluster).anech = mean(psth_strf(cluster).anech);
    psth_strf(cluster).reverb1 = [psth(21:30,:),psth(31:40,:)];
    psth_strf(cluster).reverb1 = mean(psth_strf(cluster).reverb1);
    psth_strf(cluster).reverb2 = [psth(41:50,:),psth(51:60,:)];
    psth_strf(cluster).reverb2 = mean(psth_strf(cluster).reverb2);
    psth_strf(cluster).cluster_id = cluster_id;
    
end
fprintf('== Done! Processing took %0.fs ==\n',toc);
end