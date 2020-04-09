%% Enter all the peentrations to be processed
extract = true;
plot_hist = false;
features = [1 2];
all_pen_dir = '/mnt/40086D4C086D41D0/Inhibitory_excitatory_analysis/For_analysis/All_pens';
pen{1} = '/mnt/40086D4C086D41D0/Inhibitory_excitatory_analysis/For_analysis/Ronnie/P04/P04-quning';
pen{2} = '/mnt/40086D4C086D41D0/Inhibitory_excitatory_analysis/For_analysis/Ronnie/P05/P05-quning';
pen{3} = '/mnt/40086D4C086D41D0/Inhibitory_excitatory_analysis/For_analysis/Ronnie/P06/P06-quning';
pen{4} = '/mnt/40086D4C086D41D0/Inhibitory_excitatory_analysis/For_analysis/Ronnie/P08/P08-quning_2';
pen{5} = '/mnt/40086D4C086D41D0/Inhibitory_excitatory_analysis/For_analysis/Ronnie/P09/P09-quning';
pen{6} = '/mnt/40086D4C086D41D0/Inhibitory_excitatory_analysis/For_analysis/Ronnie/P11/P11-quning_3';
pen{7} = '/mnt/40086D4C086D41D0/Inhibitory_excitatory_analysis/For_analysis/Ronnie/P13/P13-quning_2';
pen{8} = '/mnt/40086D4C086D41D0/Inhibitory_excitatory_analysis/For_analysis/Ronnie/P14/P14-quning';
pen{9} = '/mnt/40086D4C086D41D0/Inhibitory_excitatory_analysis/For_analysis/Dory/P00/P00-quning';
pen{10} = '/mnt/40086D4C086D41D0/Inhibitory_excitatory_analysis/For_analysis/Dory/P01/P01-quning_g1';
pen{11} = '/mnt/40086D4C086D41D0/Inhibitory_excitatory_analysis/For_analysis/Dory/P02/P02-quning_2';
pen{12} = '/mnt/40086D4C086D41D0/Inhibitory_excitatory_analysis/For_analysis/Dory/P04/P04-quning';
pen{13} = '/mnt/40086D4C086D41D0/Inhibitory_excitatory_analysis/For_analysis/Dory/P05/P05-quning';

%% Extract the waveforms, cluster, plot histograms and scatter plots
if extract
    for p = 1:length(pen)
        fprintf('== Processing pen %0.f/%0.f ==\n',p,length(pen));
        [clust_info] = process_wfs(pen{p},featur es);
        if plot_hist
            plot_histograms(pen{p});
            plot_corr(pen{p});
        end
    end
end
%% Combine the penetrations together
fprintf('== Getting the data ==\n');tic;
clust_info_temp.results = struct([]);
for p = 1:length(pen)
    fprintf('== Loading pen %0.f/%0.f ==\n',p,length(pen));
    load(fullfile(pen{p},'clust_info'));
    clust_info_temp.results = [clust_info_temp.results;clust_info.results];
end
clear clust_info;
names = fieldnames(clust_info_temp.results);
for n = 1:length(names)
    clust_info.results.(names{n}) = reach(clust_info_temp.results,names{n});
end
fprintf('== Done! This took %0.fs ==\n',toc);
save(fullfile(all_pen_dir,'clust_info.mat'),'clust_info');
%% Plot histograms and scatter plots for all pens together
if plot_hist
    plot_histograms(all_pen_dir);
    plot_corr(all_pen_dir);
end
%% Cluster and record metrics for all of the pens together
plot_on = 1;

file_path = fullfile(all_pen_dir,'clust_info.mat');
% load(file_path);
slash_ix1 = find(file_path == '/', 1, 'last');
var_name = file_path(slash_ix1+1:end);
var_name = strrep(var_name,'_',' ');

[clust_info] = run_kmeans(clust_info,features,plot_on);
save_dir = fullfile(all_pen_dir,['/Plots/Clustering/','feature_',num2str(features(1)),'_',num2str(features(2))]);
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end
clust_info_name = fullfile(save_dir,'clust_info.mat');
save(clust_info_name,'clust_info'); %Save the updates clust info
fig_handles = findobj('Type', 'figure');
save_name = fullfile(save_dir,[var_name,'_IE_Clustering_norm.png']);
export_fig(save_name,fig_handles(1));
save_name = fullfile(save_dir,[var_name,'_IE_Clustering.png']);
export_fig(save_name,fig_handles(2));
close all;
