%% Definte the penetrations
make_fras = 1;
min_num_neurons = 5; %How many neurons of each type there have to be to use the penetration
all_pen_dir = '/mnt/40086D4C086D41D0/Inhibitory_excitatory_analysis/For_analysis/All_pens';
pen{1} = '/mnt/40086D4C086D41D0/Inhibitory_excitatory_analysis/For_analysis/Ronnie/P04/P04-quning';
pen{2} = '/mnt/40086D4C086D41D0/Inhibitory_excitatory_analysis/For_analysis/Ronnie/P05/P05-quning';
pen{3} = '/mnt/40086D4C086D41D0/Inhibitory_excitatory_analysis/For_analysis/Ronnie/P06/P06-quning';
pen{4} = '/mnt/40086D4C086D41D0/Inhibitory_excitatory_analysis/For_analysis/Ronnie/P08/P08-quning_2';
pen{5} = '/mnt/40086D4C086D41D0/Inhibitory_excitatory_analysis/For_analysis/Ronnie/P09/P09-quning';
pen{6} = '/mnt/40086D4C086D41D0/Inhibitory_excitatory_analysis/For_analysis/Ronnie/P11/P11-quning_3';
pen{7} = '/mnt/40086D4C086D41D0/Inhibitory_excitatory_analysis/For_analysis/Ronnie/P13/P13-quning_2';
pen{8} = '/mnt/40086D4C086D41D0/Inhibitory_excitatory_analysis/For_analysis/Ronnie/P14/P14-quning';
pen{9} = '/mnt/40086D4C086D41D0/Inhibitory_excitatory_analysis/For_analysis/Dory/P01/P01-quning_g1';
pen{10} = '/mnt/40086D4C086D41D0/Inhibitory_excitatory_analysis/For_analysis/Dory/P02/P02-quning_2';
pen{11} = '/mnt/40086D4C086D41D0/Inhibitory_excitatory_analysis/For_analysis/Dory/P04/P04-quning';
pen{12} = '/mnt/40086D4C086D41D0/Inhibitory_excitatory_analysis/For_analysis/Dory/P05/P05-quning';

%% Make FRAs, compute FRA properties and make plots
if make_fras
    for p = 1:length(pen)
        compute_plot_fra(pen{p})
    end
end
%% Combine all penetrations
fprintf('== Getting the data ==\n');tic;
fra_params_temp = struct([]);
fra_params = struct([]);
per_tuned_exc = [];
per_tuned_inh = [];
for p = 1:length(pen)
    fprintf('== Loading pen %0.f/%0.f ==\n',p,length(pen));
    load(fullfile(pen{p},'Plots/','fra_psth'));
    
    if numel(fra_psth(1).inh_id_sig)>min_num_neurons && numel(fra_psth(1).exc_id_sig)>min_num_neurons %Only take penetrations that have at least x number of exc and inh tuned neurons
        fra_params_temp = [fra_params_temp;structfun(@transpose,fra_psth(1).fra_params,'UniformOutput',false)];
        per_tuned_exc = [per_tuned_exc;fra_psth(1).per_tuned_exc];
        per_tuned_inh =[per_tuned_inh;fra_psth(1).per_tuned_inh];
    end
    
end
clear fra_psth;
names = fieldnames(fra_params_temp);
for n = 1:length(names)
    fra_params(1).(names{n}) = reach(fra_params_temp,names{n});
end
fprintf('== Done! This took %0.fs ==\n',toc);
save(fullfile(all_pen_dir,'fra_params.mat'),'fra_params');
%% Compute stats
%Q10
%Perform Wilcoxon rank sum test =  Mann-Whitney U-test
[pw_q10,~] = ranksum(fra_params.q10_exc,fra_params.q10_inh);
%Perform Kolmogorov-Smirnov test to see if the two distributions are
%different
[~,pks_q10] = kstest2(fra_params.q10_exc,fra_params.q10_inh);

%Q30
%Perform Wilcoxon rank sum test =  Mann-Whitney U-test
[pw_q30,~] = ranksum(fra_params.q30_exc,fra_params.q30_inh);
%Perform Kolmogorov-Smirnov test to see if the two distributions are
%different
[~,pks_q30] = kstest2(fra_params.q30_exc,fra_params.q30_inh);

%BW10
%Perform Wilcoxon rank sum test =  Mann-Whitney U-test
[pw_bw10,~] = ranksum(fra_params.bw_10_exc,fra_params.bw_10_inh);
%Perform Kolmogorov-Smirnov test to see if the two distributions are
%different
[~,pks_bw10] = kstest2(fra_params.bw_10_exc,fra_params.bw_10_inh);

%BW30
%Perform Wilcoxon rank sum test =  Mann-Whitney U-test
[pw_bw30,~] = ranksum(fra_params.bw_30_exc,fra_params.bw_30_inh);
%Perform Kolmogorov-Smirnov test to see if the two distributions are
%different
[~,pks_bw30] = kstest2(fra_params.bw_30_exc,fra_params.bw_30_inh);

%dB Thresholds
%Perform Wilcoxon rank sum test =  Mann-Whitney U-test
[pw_db_th,~] = ranksum(fra_params.db_th_exc,fra_params.db_th_inh);
%Perform Kolmogorov-Smirnov test to see if the two distributions are
%different
[~,pks_db_th] = kstest2(fra_params.db_th_exc,fra_params.db_th_inh);

%Percentage of tuned neurons
%Perform Wilcoxon rank sum test =  Mann-Whitney U-test
[p_tune,~] = ranksum(per_tuned_exc,per_tuned_inh);
%% Make plots
fontsize = 20;
tbox_size = 12;
fprintf('== Computing FRA stats and making figures ==\n');
save_dir_stats = fullfile(all_pen_dir,'/Stats');
if ~exist(save_dir_stats, 'dir')
    mkdir(save_dir_stats);
end

%Q10 Box plots
g1 = repmat({'Excitatory'},numel(fra_params.q10_exc),1);
g2 = repmat({'Inhibitory'},numel(fra_params.q10_inh),1);
g = [g1; g2];
figure('units','normalized','outerposition',[0 0 1 1]);
boxplot([fra_params.q10_exc;fra_params.q10_inh],g,'Notch','on');
title('Q10 for excitatory and inhibitory neurons');
ylabel('Q10 [AU]');
ylim([0 5]);
annotation('textbox', [0.75, 0.8, 0.1, 0.1], 'String', ['Mann-Whitney U-test p=',num2str(pw_q10,'%.2f')],'FontWeight','Bold','FontSize',tbox_size,'Linestyle','none');
set(gcf,'color','w');
set(gca,'FontName','Arial','FontSize',fontsize,'FontWeight','Bold');
file_name = ['Q10_boxplots'];
save_name = fullfile(save_dir_stats,[file_name,'.png']);
export_fig(save_name);
close;


%Q10 Histograms
figure('units','normalized','outerposition',[0 0 1 1]);
edges = [0:0.2:3];
h_exc = histogram(fra_params.q10_exc,edges);
h_exc.FaceColor = 'r';
h_exc.Normalization = 'probability';
hold on;
h_inh = histogram(fra_params.q10_inh,edges);
h_inh.FaceColor = 'b';
h_inh.Normalization = 'probability';
title('Distribution of Q10 for excitatory and inhibitory neurons');
ylabel('Probability');
xlabel('Q10');
legend('Putative Excitatory','Putative Inhibitory')
annotation('textbox', [0.73, 0.74, 0.1, 0.1], 'String', ['Kolmogorov-Smirnov test p=',num2str(pks_q10,'%.2f')],'FontWeight','Bold','FontSize',tbox_size,'Linestyle','none');
set(gcf,'color','w');
set(gca,'FontName','Arial','FontSize',fontsize,'FontWeight','Bold');
file_name = ['Q10_histograms'];
save_name = fullfile(save_dir_stats,[file_name,'.png']);
export_fig(save_name);
close;


%Q30 Box plots
g1 = repmat({'Excitatory'},numel(fra_params.q30_exc),1);
g2 = repmat({'Inhibitory'},numel(fra_params.q30_inh),1);
g = [g1; g2];
figure('units','normalized','outerposition',[0 0 1 1]);
boxplot([fra_params.q30_exc;fra_params.q30_inh],g,'Notch','on');
title('Q30 for excitatory and inhibitory neurons');
ylabel('Q30 [AU]');
ylim([0 5]);
annotation('textbox', [0.75, 0.8, 0.1, 0.1], 'String', ['Mann-Whitney U-test p=',num2str(pw_q30,'%.2f')],'FontWeight','Bold','FontSize',tbox_size,'Linestyle','none');
set(gcf,'color','w');
set(gca,'FontName','Arial','FontSize',fontsize,'FontWeight','Bold');
file_name = ['Q30_boxplots'];
save_name = fullfile(save_dir_stats,[file_name,'.png']);
export_fig(save_name);
close;


%Q30 Histograms
figure('units','normalized','outerposition',[0 0 1 1]);
edges = [0:0.2:3];
h_exc = histogram(fra_params.q30_exc,edges);
h_exc.FaceColor = 'r';
h_exc.Normalization = 'probability';
hold on;
h_inh = histogram(fra_params.q30_inh,edges);
h_inh.FaceColor = 'b';
h_inh.Normalization = 'probability';
title('Distribution of Q30 for excitatory and inhibitory neurons');
ylabel('Probability');
xlabel('Q30');
legend('Putative Excitatory','Putative Inhibitory')
annotation('textbox', [0.73, 0.74, 0.1, 0.1], 'String', ['Kolmogorov-Smirnov test p=',num2str(pks_q30,'%.2f')],'FontWeight','Bold','FontSize',tbox_size,'Linestyle','none');
set(gcf,'color','w');
set(gca,'FontName','Arial','FontSize',fontsize,'FontWeight','Bold');
file_name = ['Q30_histograms'];
save_name = fullfile(save_dir_stats,[file_name,'.png']);
export_fig(save_name);
close;


%BW10 Box plots
g1 = repmat({'Excitatory'},numel(fra_params.bw_10_exc),1);
g2 = repmat({'Inhibitory'},numel(fra_params.bw_10_inh),1);
g = [g1; g2];
figure('units','normalized','outerposition',[0 0 1 1]);
boxplot([fra_params.bw_10_exc;fra_params.bw_10_inh],g,'Notch','on');
title('BW10 for excitatory and inhibitory neurons');
ylabel('BW10 [kHz]');
annotation('textbox', [0.75, 0.8, 0.1, 0.1], 'String', ['Mann-Whitney U-test p=',num2str(pw_bw10,'%.2f')],'FontWeight','Bold','FontSize',tbox_size,'Linestyle','none');
set(gcf,'color','w');
set(gca,'FontName','Arial','FontSize',fontsize,'FontWeight','Bold');
file_name = ['BW10_boxplots'];
save_name = fullfile(save_dir_stats,[file_name,'.png']);
export_fig(save_name);
close;


%BW10 Histograms
figure('units','normalized','outerposition',[0 0 1 1]);
edges = [0:2:30];
h_exc = histogram(fra_params.bw_10_exc,edges);
h_exc.FaceColor = 'r';
h_exc.Normalization = 'probability';
hold on;
h_inh = histogram(fra_params.bw_10_inh,edges);
h_inh.FaceColor = 'b';
h_inh.Normalization = 'probability';
title('Distribution of BW10 for excitatory and inhibitory neurons');
ylabel('Probability');
xlabel('BW10 [kHz]');
legend('Putative Excitatory','Putative Inhibitory')
annotation('textbox', [0.73, 0.74, 0.1, 0.1], 'String', ['Kolmogorov-Smirnov test p=',num2str(pks_bw10,'%.2f')],'FontWeight','Bold','FontSize',tbox_size,'Linestyle','none');
set(gcf,'color','w');
set(gca,'FontName','Arial','FontSize',fontsize,'FontWeight','Bold');
file_name = ['BW10_histograms'];
save_name = fullfile(save_dir_stats,[file_name,'.png']);
export_fig(save_name);
close;

%BW30 Box plots
g1 = repmat({'Excitatory'},numel(fra_params.bw_30_exc),1);
g2 = repmat({'Inhibitory'},numel(fra_params.bw_30_inh),1);
g = [g1; g2];
figure('units','normalized','outerposition',[0 0 1 1]);
boxplot([fra_params.bw_30_exc;fra_params.bw_30_inh],g,'Notch','on');
title('BW30 for excitatory and inhibitory neurons');
ylabel('BW30 [kHz]');
annotation('textbox', [0.75, 0.8, 0.1, 0.1], 'String', ['Mann-Whitney U-test p=',num2str(pw_bw30,'%.2f')],'FontWeight','Bold','FontSize',tbox_size,'Linestyle','none');
set(gcf,'color','w');
set(gca,'FontName','Arial','FontSize',fontsize,'FontWeight','Bold');
file_name = ['BW30_boxplots'];
save_name = fullfile(save_dir_stats,[file_name,'.png']);
export_fig(save_name);
close;


%BW30 Histograms
figure('units','normalized','outerposition',[0 0 1 1]);
edges = [0:2:30];
h_exc = histogram(fra_params.bw_30_exc,edges);
h_exc.FaceColor = 'r';
h_exc.Normalization = 'probability';
hold on;
h_inh = histogram(fra_params.bw_30_inh,edges);
h_inh.FaceColor = 'b';
h_inh.Normalization = 'probability';
title('Distribution of BW30 for excitatory and inhibitory neurons');
ylabel('Probability');
xlabel('BW30 [kHz]');
legend('Putative Excitatory','Putative Inhibitory')
annotation('textbox', [0.73, 0.74, 0.1, 0.1], 'String', ['Kolmogorov-Smirnov test p=',num2str(pks_bw30,'%.2f')],'FontWeight','Bold','FontSize',tbox_size,'Linestyle','none');
set(gcf,'color','w');
set(gca,'FontName','Arial','FontSize',fontsize,'FontWeight','Bold');
file_name = ['BW30_histograms'];
save_name = fullfile(save_dir_stats,[file_name,'.png']);
export_fig(save_name);
close;

%dB thresholds Box plots
g1 = repmat({'Excitatory'},numel(fra_params.db_th_exc),1);
g2 = repmat({'Inhibitory'},numel(fra_params.db_th_inh),1);
g = [g1; g2];
figure('units','normalized','outerposition',[0 0 1 1]);
boxplot([fra_params.db_th_exc;fra_params.db_th_inh],g,'Notch','on','Labels',{'Excitatory','Inhibitory'});
title('Threshold for excitatory and inhibitory neurons');
ylabel('Threshold [dB]');
ylim([20 80]);
annotation('textbox', [0.75, 0.8, 0.1, 0.1], 'String', ['Mann-Whitney U-test p=',num2str(pw_db_th,'%.2f')],'FontWeight','Bold','FontSize',tbox_size,'Linestyle','none');
set(gcf,'color','w');
set(gca,'FontName','Arial','FontSize',fontsize,'FontWeight','Bold');
file_name = ['dB_Thresholds_boxplots'];
save_name = fullfile(save_dir_stats,[file_name,'.png']);
export_fig(save_name);
close;

%dB thresholds Histograms
figure('units','normalized','outerposition',[0 0 1 1]);
edges = [35:10:80];
h_exc = histogram(fra_params.db_th_exc,edges);
h_exc.FaceColor = 'r';
h_exc.Normalization = 'probability';
hold on;
h_inh = histogram(fra_params.db_th_inh,edges);
h_inh.FaceColor = 'b';
h_inh.Normalization = 'probability';
title('Distribution of thresholds for excitatory and inhibitory neurons');
ylabel('Probability');
xlabel('Threshold [dB]');
legend('Putative Excitatory','Putative Inhibitory')
annotation('textbox', [0.73, 0.74, 0.1, 0.1], 'String', ['Kolmogorov-Smirnov test p=',num2str(pks_db_th,'%.2f')],'FontWeight','Bold','FontSize',tbox_size,'Linestyle','none');
set(gcf,'color','w');
set(gca,'FontName','Arial','FontSize',fontsize,'FontWeight','Bold');
file_name = ['Thresholds_histograms'];
save_name = fullfile(save_dir_stats,[file_name,'.png']);
export_fig(save_name);
close;
%% Save results