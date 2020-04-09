function compute_plot_fra(sorted_dir)

%% Set up variables
qualityOpt = 'Good'; %Which clusters to use for the FRAs: 'Good' - Only the Good units; 'MUA' - Only the MUA units; 'Both' - Both Good and MUA (no noise)
p_val_thr = 0.05; %The pvalue threshold in the ANOVA test that we use to select FRAs

[temp,stim_name] = fileparts(sorted_dir);
[temp2,pen_name] = fileparts(temp);
[~,animal_name] = fileparts(temp2);
synch_dir = fullfile('/mnt/40086D4C086D41D0/Inhibitory_excitatory_analysis/Meta/synch_ch/',animal_name,pen_name); %This is the directory of synchronization channel that we use to know when stimuli were presented relative to the spikes
grid_dir = fullfile('/mnt/40086D4C086D41D0/Inhibitory_excitatory_analysis/Meta/benware_grids/',animal_name,stim_name); %This is the directory of the benware file that tells us what was the order of stimulus presentation for the FRAs (pure tones)
%% Get the clust_info
fprintf('== Getting the spike times ==\n');tic;
load(fullfile(sorted_dir,'clust_info'));
fprintf('== Done! This took %0.fs ==\n',toc);
%% Make the FRA psth
fprintf('== Generating the FRA PSTHs ==\n');tic;
plot_folder = fullfile(sorted_dir,'Plots'); %Points to where we wanna save the plots
fname = fullfile(plot_folder,'fra_psth.mat');

if ~isfile(fname)
    fra_psth = fra_pixels2(synch_dir,grid_dir,sorted_dir,qualityOpt);
else
    load(fname);
end

save_dir_exc = fullfile(plot_folder,'/FRA_excitatory');
save_dir_inh = fullfile(plot_folder,'/FRA_inhibitory');
if ~exist(save_dir_exc, 'dir') || ~exist(save_dir_inh, 'dir')
    mkdir(save_dir_exc);
    mkdir(save_dir_inh);
end

fprintf('== Done! This took %0.fs ==\n',toc);
%% Sort neurons by p-val
temp_p = fra_psth(1).pval; %Create temp var to work on 
temp_p_sort = sortrows(temp_p,1,'ascend'); %Sort based on the pval for freqeuncy from anova test 
cluster_id_sort = temp_p_sort(:,5); %Get the sorted ids of all clusters
p_val_freq_sort = temp_p_sort(:,1); %Get the sorted freq p vals for all clsuters
p_val_int_sort = temp_p_sort(:,3); %Get the sorted freq-level interaction p vals for all clsuters
clusters_id_final = cluster_id_sort(p_val_freq_sort<p_val_thr | p_val_int_sort<p_val_thr); %Select only clusters that pass the freq p val threshold
%% Split into excitatory and inhibitory
inh_flag = clust_info.clustering.inh_flag;
exc_id = clust_info.clustering.idx(clust_info.clustering.idx(:,2)~=inh_flag,1);
inh_id = clust_info.clustering.idx(clust_info.clustering.idx(:,2)==inh_flag,1);

exc_id_sig = intersect(clusters_id_final,exc_id,'stable');
inh_id_sig = intersect(clusters_id_final,inh_id,'stable');

per_tuned_exc = numel(exc_id_sig)/numel(exc_id);
per_tuned_inh = numel(inh_id_sig)/numel(inh_id);
%% Estimate BF, CF, Q10 & Q30 for excitatory and inhibitory separately
fprintf('== Computing FRA properties ==\n');tic;
params = fra_psth(1).params;
cluster_ids = fra_psth(1).cluster_id;
fra_psth(1).exc_id_sig = exc_id_sig;
fra_psth(1).per_tuned_exc = per_tuned_exc;
for c = 1:numel(exc_id_sig)
    cluster = exc_id_sig(c);
    cluster_ix = find(cluster_ids==cluster);
    [fra_psth(c).fra_prop_exc] = compute_fra_prop(fra_psth(cluster_ix).X_dbft,fra_psth(cluster_ix).base_rate,params);
    if ~isempty(fra_psth(c).fra_prop_exc.q10)
        fra_params.q10_exc(c) = fra_psth(c).fra_prop_exc.q10;
        fra_params.bw_10_exc(c) = fra_psth(c).fra_prop_exc.bw_10;
    end
    if ~isempty(fra_psth(c).fra_prop_exc.q30)
        fra_params.q30_exc(c) = fra_psth(c).fra_prop_exc.q30;
        fra_params.bw_30_exc(c) = fra_psth(c).fra_prop_exc.bw_30;
    end
    if ~isempty(fra_psth(c).fra_prop_exc.db_th)
        fra_params.db_th_exc(c) = fra_psth(c).fra_prop_exc.db_th;
        fra_params.bf_exc(c) = fra_psth(c).fra_prop_exc.bf;
        fra_params.cf_exc(c) = fra_psth(c).fra_prop_exc.cf;
    end
end

fra_psth(1).inh_id_sig = inh_id_sig;
fra_psth(1).per_tuned_inh= per_tuned_inh;
for c = 1:numel(inh_id_sig)
    cluster = inh_id_sig(c);
    cluster_ix = find(cluster_ids==cluster);
    [fra_psth(c).fra_prop_inh] = compute_fra_prop(fra_psth(cluster_ix).X_dbft,fra_psth(cluster_ix).base_rate,params);
    if ~isempty(fra_psth(c).fra_prop_inh.q10)
        fra_params.q10_inh(c) = fra_psth(c).fra_prop_inh.q10;
        fra_params.bw_10_inh(c) = fra_psth(c).fra_prop_inh.bw_10;
    end
    if ~isempty(fra_psth(c).fra_prop_inh.q30)
        fra_params.q30_inh(c) = fra_psth(c).fra_prop_inh.q30;
        fra_params.bw_30_inh(c) = fra_psth(c).fra_prop_inh.bw_30;
    end
    
    if ~isempty(fra_psth(c).fra_prop_inh.db_th)
        fra_params.db_th_inh(c) = fra_psth(c).fra_prop_inh.db_th;
        fra_params.bf_inh(c) = fra_psth(c).fra_prop_inh.bf;
        fra_params.cf_inh(c) = fra_psth(c).fra_prop_inh.cf;
    end
end
fprintf('== Done! This took %0.fs ==\n',toc);
%% Plot the FRAs
fprintf('== Plotting the FRAs ==\n');
plot_fra_properties(fra_psth,'exc',save_dir_exc);
plot_fra_properties(fra_psth,'inh',save_dir_inh);
%% Compute stats on FRA properties and plot box plots
fontsize = 20;
fprintf('== Computing FRA stats and making figures ==\n');
save_dir_stats = fullfile(plot_folder,'/Stats');
if ~exist(save_dir_stats, 'dir')
    mkdir(save_dir_stats);
end

%Perform Wilcoxon rank sum test =  Mann-Whitney U-test
[stats.q10_p_wilcox,hw] = ranksum(fra_params.q10_exc,fra_params.q10_inh);
%Perform Kolmogorov-Smirnov test to see if the two distributions are
%different
% [hks,pks] = kstest2(fra_params.q10_exc,fra_params.q10_inh);
%Plot the Q10
g1 = repmat({'Excitatory'},numel(fra_params.q10_exc),1);
g2 = repmat({'Inhibitory'},numel(fra_params.q10_inh),1);
g = [g1; g2];
figure('units','normalized','outerposition',[0 0 1 1]);
boxplot([fra_params.q10_exc';fra_params.q10_inh'],g,'Notch','on');
title('Q10 for excitatory and inhibitory neurons');
ylabel('Q10 [AU]');
ylim([0 5]);
annotation('textbox', [0.75, 0.8, 0.1, 0.1], 'String', ['Mann-Whitney U-test p=',num2str(stats.q10_p_wilcox,'%.2f')]);
set(gcf,'color','w');
set(gca,'FontName','Arial','FontSize',fontsize,'FontWeight','Bold');
file_name = ['Q10_boxplots'];
save_name = fullfile(save_dir_stats,[file_name,'.png']);
export_fig(save_name);
close;


%Plot a histogram fo the two distributions
% figure('units','normalized','outerposition',[0 0 1 1]);
% edges = [0:0.2:5];
% h_exc = histogram(fra_params.q10_exc,edges);
% h_exc.FaceColor = 'r';
% h_exc.Normalization = 'probability';
% hold on;
% h_inh = histogram(fra_params.q10_inh,edges);
% h_inh.FaceColor = 'b';
% h_inh.Normalization = 'probability';
% title('Distribution of Q10 for excitatory and inhibitory neurons');
% ylabel('Probability');
% xlabel('Q10');
% legend('Putative Excitatory','Putative Inhibitory')
% annotation('textbox', [0.665, 0.81, 0.1, 0.1], 'String', ['Kolmogorov-Smirnov test p=',num2str(pks,'%.2f')]);


%Perform Wilcoxon rank sum test =  Mann-Whitney U-test
[stats.q30_p_wilcox,hw] = ranksum(fra_params.q30_exc,fra_params.q30_inh);
%Plot the Q30
g1 = repmat({'Excitatory'},numel(fra_params.q30_exc),1);
g2 = repmat({'Inhibitory'},numel(fra_params.q30_inh),1);
g = [g1; g2];
figure('units','normalized','outerposition',[0 0 1 1]);
boxplot([fra_params.q30_exc';fra_params.q30_inh'],g,'Notch','on');
title('Q30 for excitatory and inhibitory neurons');
ylabel('Q30 [AU]');
ylim([0 5]);
annotation('textbox', [0.75, 0.8, 0.1, 0.1], 'String', ['Mann-Whitney U-test p=',num2str(stats.q30_p_wilcox,'%.2f')]);
set(gcf,'color','w');
set(gca,'FontName','Arial','FontSize',fontsize,'FontWeight','Bold');
file_name = ['Q30_boxplots'];
save_name = fullfile(save_dir_stats,[file_name,'.png']);
export_fig(save_name);
close;


%Perform Wilcoxon rank sum test =  Mann-Whitney U-test
[stats.bw10_p_wilcox,hw] = ranksum(fra_params.bw_10_exc,fra_params.bw_10_inh);
%Plot the BW10
g1 = repmat({'Excitatory'},numel(fra_params.bw_10_exc),1);
g2 = repmat({'Inhibitory'},numel(fra_params.bw_10_inh),1);
g = [g1; g2];
figure('units','normalized','outerposition',[0 0 1 1]);
boxplot([fra_params.bw_10_exc';fra_params.bw_10_inh'],g,'Notch','on','Labels',{'Excitatory','Inhibitory'});
title('BW10 for excitatory and inhibitory neurons');
ylabel('BW10 [kHz]');
annotation('textbox', [0.75, 0.8, 0.1, 0.1], 'String', ['Mann-Whitney U-test p=',num2str(stats.bw10_p_wilcox,'%.2f')]);
set(gcf,'color','w');
set(gca,'FontName','Arial','FontSize',fontsize,'FontWeight','Bold');
file_name = ['BW10_boxplots'];
save_name = fullfile(save_dir_stats,[file_name,'.png']);
export_fig(save_name);
close;


%Perform Wilcoxon rank sum test =  Mann-Whitney U-test
[stats.bw30_p_wilcox,hw] = ranksum(fra_params.bw_30_exc,fra_params.bw_30_inh);
%Plot the BW30
g1 = repmat({'Excitatory'},numel(fra_params.bw_30_exc),1);
g2 = repmat({'Inhibitory'},numel(fra_params.bw_30_inh),1);
g = [g1; g2];
figure('units','normalized','outerposition',[0 0 1 1]);
boxplot([fra_params.bw_30_exc';fra_params.bw_30_inh'],g,'Notch','on','Labels',{'Excitatory','Inhibitory'});
title('BW30 for excitatory and inhibitory neurons');
ylabel('BW30 [kHz]');
annotation('textbox', [0.75, 0.8, 0.1, 0.1], 'String', ['Mann-Whitney U-test p=',num2str(stats.bw30_p_wilcox,'%.2f')]);
set(gcf,'color','w');
set(gca,'FontName','Arial','FontSize',fontsize,'FontWeight','Bold');
file_name = ['BW30_boxplots'];
save_name = fullfile(save_dir_stats,[file_name,'.png']);
export_fig(save_name);
close;


%Perform Wilcoxon rank sum test =  Mann-Whitney U-test
[stats.db_th_p_wilcox,hw] = ranksum(fra_params.db_th_exc,fra_params.db_th_inh);
%Plot the dB thresholds
g1 = repmat({'Excitatory'},numel(fra_params.db_th_exc),1);
g2 = repmat({'Inhibitory'},numel(fra_params.db_th_inh),1);
g = [g1; g2];
figure('units','normalized','outerposition',[0 0 1 1]);
boxplot([fra_params.db_th_exc';fra_params.db_th_inh'],g,'Notch','on','Labels',{'Excitatory','Inhibitory'});
title('dB Threshold for excitatory and inhibitory neurons');
ylabel('Loudness [dB]');
annotation('textbox', [0.75, 0.8, 0.1, 0.1], 'String', ['Mann-Whitney U-test p=',num2str(stats.db_th_p_wilcox,'%.2f')]);
set(gcf,'color','w');
set(gca,'FontName','Arial','FontSize',fontsize,'FontWeight','Bold');
file_name = ['dB_Thresholds_boxplots'];
save_name = fullfile(save_dir_stats,[file_name,'.png']);
export_fig(save_name);
close;


%Perform Wilcoxon rank sum test =  Mann-Whitney U-test
[stats.cf_p_wilcox,hw] = ranksum(fra_params.cf_exc,fra_params.cf_inh);
%Plot the CF
g1 = repmat({'Excitatory'},numel(fra_params.cf_exc),1);
g2 = repmat({'Inhibitory'},numel(fra_params.cf_inh),1);
g = [g1; g2];
figure('units','normalized','outerposition',[0 0 1 1]);
boxplot(log2([fra_params.cf_exc';fra_params.cf_inh']),g,'Notch','on','Labels',{'Excitatory','Inhibitory'});
title('CF for excitatory and inhibitory neurons');
ylabel('freq [log2(kHz)]');
annotation('textbox', [0.75, 0.8, 0.1, 0.1], 'String', ['Mann-Whitney U-test p=',num2str(stats.cf_p_wilcox,'%.2f')]);
set(gcf,'color','w');
set(gca,'FontName','Arial','FontSize',fontsize,'FontWeight','Bold');
file_name = ['CF_boxplots'];
save_name = fullfile(save_dir_stats,[file_name,'.png']);
export_fig(save_name);
close;


%Perform Wilcoxon rank sum test =  Mann-Whitney U-test
[stats.bf_p_wilcox,hw] = ranksum(fra_params.bf_exc,fra_params.bf_inh);
%Plot the BF
g1 = repmat({'Excitatory'},numel(fra_params.bf_exc),1);
g2 = repmat({'Inhibitory'},numel(fra_params.bf_inh),1);
g = [g1; g2];
figure('units','normalized','outerposition',[0 0 1 1]);
boxplot(log2([fra_params.bf_exc';fra_params.bf_inh']),g,'Notch','on','Labels',{'Excitatory','Inhibitory'});
title('BF for excitatory and inhibitory neurons');
ylabel('freq [log2(kHz)]');
annotation('textbox', [0.75, 0.8, 0.1, 0.1], 'String', ['Mann-Whitney U-test p=',num2str(stats.cf_p_wilcox,'%.2f')]);
set(gcf,'color','w');
set(gca,'FontName','Arial','FontSize',fontsize,'FontWeight','Bold');
file_name = ['BF_boxplots'];
save_name = fullfile(save_dir_stats,[file_name,'.png']);
export_fig(save_name);
close;

%% Save the results
fprintf('== Saving the results ==\n');
fra_psth(1).stats = stats;
fra_psth(1).fra_params = fra_params;
save(fname,'fra_psth');