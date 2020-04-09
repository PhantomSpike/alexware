function plot_tconst_neurons_grande(kernel_dir,NPSP_th,get_neurons)
%% Params
chop_ms = 190;
hist_type = 'bar';
if ~exist('NPSP_th','var') || isempty(NPSP_th)
    NPSP_th = 40;
end

if ~exist('get_neurons','var') || isempty(get_neurons)
    get_neurons = 'all';
end

[save_dir,~] = fileparts(kernel_dir);
plot_dir = fullfile(save_dir,'Plots');
if ~exist(plot_dir,'dir')
    mkdir(plot_dir);
end

save_dir_full = fullfile(plot_dir,get_neurons);
if ~exist(save_dir_full,'dir')
    mkdir(save_dir_full);
end

r_type{1} = 'anech';
r_type{2} = 'small';
r_type{3} = 'big';
n_rooms = length(r_type);
%% Load the data
load(fullfile(kernel_dir,'info'),'info');
temp_files = dir([kernel_dir,'/*.mat']);
temp = load(fullfile(temp_files(1).folder,temp_files(1).name));
model = temp.kernel.model;
%% Select the data
switch get_neurons
    case 'all'
        ix_qualia = ones(length(info.cluster_id),1); %Get all the neurons
    case 'good'
        ix_qualia = cell2mat(cellfun(@(x) strcmp(x,'good'),info.quality,'UniformOutput',false)); %Find all the good units
    case 'mua'
        ix_qualia = cell2mat(cellfun(@(x) strcmp(x,'mua'),info.quality,'UniformOutput',false)); %Find all the mua units
end
ix_npsp = info.NPSP<NPSP_th; %Find all the neurons below certain NPSP
ix = ix_qualia & ix_npsp; %Find the intersection of the two

NPSPs = info.NPSP(ix);
cluster_ids = info.cluster_id(ix);
animal_names = info.animal_name(ix);
pen_names = info.pen_name(ix);
qualities = info.quality(ix);

%Sort in increasing NPSP
[NPSPs,ix_select] = sort(NPSPs,'ascend');
animal_names = animal_names(ix_select);
pen_names = pen_names(ix_select);
cluster_ids = cluster_ids(ix_select);
qualities = qualities(ix_select);
n_clust = length(cluster_ids);
%% First load all the selected kernels
fprintf('== Loading the data ==\n');tic;
for k = 1:n_clust
    c_name = fullfile(kernel_dir,strjoin({animal_names{k},pen_names{k},num2str(cluster_ids(k))},'_'));
    load(c_name,'kernel');
    kernels{k,1} = kernel;
end
fprintf('== Done! This took %0.fs ==\n',toc);
%% First compute the bf and tau for every cluster
freqs = fliplr(kernels{1}.freqs); %Get the freqeuncies but flip them because they are going low->high and cochleagram is high->low
n_h = kernel.n_h;
dt_ms = round(kernel.dt_ms);
chop_ix = round(chop_ms/dt_ms);
h = (1:1:chop_ix)';
h = dt_ms*h;
fprintf('== Calcuating tau and bf ==\n');tic;
for k = 1:n_clust
    sprintf('== Cluster %0.f/%0.f ==\n',k,n_clust);
    for r = 1:n_rooms
        room = r_type{r};
        switch model
            case {'sep','sep_kh'}      
                [~,ix] = max(kernels{k}.(room).k_f);
                bf(k).(room) = freqs(ix); %Find the corresponding frequency
                k_h = flipud(kernels{k}.(room).k_h); %Get the k_h
                
            case {'ridge','lasso','elastic'}
                k_fh = fliplr(kernels{k}.(room).main{11}.k_fh);
                k_fh = k_fh(:,1:chop_ix);
                k_fh_neg = abs(min(k_fh,0));
                k_fh_pos = abs(max(k_fh,0));
                k_h_neg = mean(k_fh_neg);
                k_h_pos = mean(k_fh_pos);
                k_f = mean(k_fh_pos,2); %Take the mean across history steps
                [~,ix] = max(k_f); %Find the index of the max frequency
                bf(k).(room) = freqs(ix); %Find the corresponding frequency
        end
        k_h_neg = k_h_neg./sum(k_h_neg(:)); %Scale the values to sum to 1 for inhibition
        k_h_pos = k_h_pos./sum(k_h_pos(:)); %Scale the values to sum to 1 for inhibition
        tau_neg(k).(room) = (k_h_neg*h); %Compute a weighted sum of all values
        tau_pos(k).(room) = (k_h_pos*h); %Compute a weighted sum of all values
    end 
    bf_mean(k) = mean([bf(k).small,bf(k).big]);
    [~,ix_bf] = min(abs(bf_mean(k) - freqs)); %Find the closest freqeuncy to the mean one  from the actual freqs
    bf_closest(k) = freqs(ix_bf);
end
fprintf('== Done! This took %0.fs ==\n',toc);

%% Stats
tau_anech_neg = reach(tau_neg,'anech');
tau_anech_pos = reach(tau_pos,'anech');
tau_small_neg = reach(tau_neg,'small');
tau_small_pos = reach(tau_pos,'small');
tau_big_neg = reach(tau_neg,'big');
tau_big_pos = reach(tau_pos,'big');
[p_val_neg,~,~] = signrank(tau_big_neg,tau_small_neg);
[p_val_pos,~,~] = signrank(tau_big_pos,tau_small_pos);
%% Plot the small vs big tau neg with NPSP 
font_type = 'Liberation Sans';
sz = 35;
y_font_sz = sz;
x_font_sz = sz;
all_font_sz = sz;

lim_val_ms = 150;
c = NPSPs;
figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);
scatter(tau_small_neg,tau_big_neg,sz,c,'filled');
xlabel('\tau_{small} [ms]','FontSize',x_font_sz,'FontWeight','bold');
ylabel('\tau_{large} [ms]','FontSize',y_font_sz,'FontWeight','bold');
title(['Inhibitory \tau_{big} vs \tau_{small} for all neurons NPSP<',num2str(NPSP_th),' ',model,' kernel']);
annotation('textbox',[0.65 0.07 0.1 0.1],'String', sprintf('n=%0.f',n_clust),'LineStyle','none','Color','k','FontSize',all_font_sz,'FontWeight','bold');
annotation('textbox',[0.73 0.86 0.1 0.1],'String', sprintf('NPSP'),'LineStyle','none','Color','k','FontSize',all_font_sz,'FontWeight','bold');
colormap('jet');
colorbar;
axis equal;
l = [0 lim_val_ms];
xlim(l);
ylim(l);
hline = refline(1,0);
hline.Color = 'k';
hline.LineWidth = 3;
hline.LineStyle = '--';
set(gca,'FontSize',all_font_sz,'FontWeight','Bold');
set(gcf,'color','w');
save_name = fullfile(save_dir_full,['Inhibitory tau for ',get_neurons,' neurons NPSP<',num2str(NPSP_th),' ',model,' kernel.png']);
export_fig(save_name);
close all;
%% Plot the small vs big pos with NPSP
sz = 30;
c = NPSPs;
figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);
scatter(tau_small_pos,tau_big_pos,sz,c,'filled');
xlabel('\tau_{small} [ms]','FontSize',x_font_sz,'FontWeight','bold');
ylabel('\tau_{large} [ms]','FontSize',y_font_sz,'FontWeight','bold');
title(['Excitatory \tau_{big} vs \tau_{small} for all neurons NPSP<',num2str(NPSP_th),' ',model,' kernel']);
annotation('textbox',[0.65 0.07 0.1 0.1],'String', sprintf('n=%0.f',n_clust),'LineStyle','none','Color','k','FontSize',all_font_sz,'FontWeight','bold');
annotation('textbox',[0.73 0.86 0.1 0.1],'String', sprintf('NPSP'),'LineStyle','none','Color','k','FontSize',all_font_sz,'FontWeight','bold');
colormap('jet');
colorbar;
axis equal;
l = [0 lim_val_ms];
xlim(l);
ylim(l);
hline = refline(1,0);
hline.Color = 'k';
hline.LineWidth = 3;
hline.LineStyle = '--';
set(gca,'FontSize',all_font_sz,'FontWeight','Bold');
set(gcf,'color','w');
save_name = fullfile(save_dir_full,['Excitatory tau for ',get_neurons,' neurons NPSP<',num2str(NPSP_th),' ',model,' kernel.png']);
export_fig(save_name);
close all;
%% Plot the small vs big inhibitory tau with bf as colour
sz = 30;
c = log2(bf_closest);
figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);
scatter(tau_small_neg,tau_big_neg,sz,c,'filled');
xlabel('\tau_{small} [ms]','FontSize',x_font_sz,'FontWeight','bold');
ylabel('\tau_{large} [ms]','FontSize',y_font_sz,'FontWeight','bold');
title(['Inhibitory \tau_{big} vs \tau_{small} for all neurons NPSP<',num2str(NPSP_th),' ',model,' kernel']);
annotation('textbox',[0.65 0.07 0.1 0.1],'String', sprintf('n=%0.f',n_clust),'LineStyle','none','Color','k','FontSize',all_font_sz,'FontWeight','bold');
annotation('textbox',[0.73 0.9 0.1 0.1],'String', sprintf('Frequency [kHz]'),'LineStyle','none','Color','k','FontSize',all_font_sz,'FontWeight','bold');
colormap('jet');
colorbar;
axis equal;
l = [0 lim_val_ms];
xlim(l);
ylim(l);
hline = refline(1,0);
hline.Color = 'k';
hline.LineWidth = 3;
hline.LineStyle = '--';
set(gca,'FontSize',all_font_sz,'FontWeight','Bold');
set(gcf,'color','w');
save_name = fullfile(save_dir_full,['Inhibitory tau vs bf for ',get_neurons,' neurons NPSP<',num2str(NPSP_th),' ',model,' kernel.png']);
export_fig(save_name);
close all;
%% Plot the small vs big excitatory tau with bf as colour
axis_sz = 20;
sz = 30;
c = log2(bf_closest);
figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);
scatter(tau_small_pos,tau_big_pos,sz,c,'filled');
xlabel('\tau_{small} [ms]','FontSize',x_font_sz,'FontWeight','bold');
ylabel('\tau_{large} [ms]','FontSize',y_font_sz,'FontWeight','bold');
title(['Excitatory \tau_{big} vs \tau_{small} for all neurons NPSP<',num2str(NPSP_th),' ',model,' kernel']);
annotation('textbox',[0.65 0.07 0.1 0.1],'String', sprintf('n=%0.f',n_clust),'LineStyle','none','Color','k','FontSize',all_font_sz,'FontWeight','bold');
annotation('textbox',[0.73 0.9 0.1 0.1],'String', sprintf('Frequency [kHz]'),'LineStyle','none','Color','k','FontSize',all_font_sz,'FontWeight','bold');
colormap('jet');
colorbar;
axis equal;
l = [0 lim_val_ms];
xlim(l);
ylim(l);
hline = refline(1,0);
hline.Color = 'k';
hline.LineWidth = 3;
hline.LineStyle = '--';
set(gca,'FontSize',all_font_sz,'FontWeight','Bold');
set(gcf,'color','w');
save_name = fullfile(save_dir_full,['Excitatory tau vs bf for ',get_neurons,' neurons NPSP<',num2str(NPSP_th),' ',model,' kernel.png']);
export_fig(save_name);
close all;
%% Plot the small vs big inhibitory and excitatory tau w/o NPSP
sz = 30;
figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);
scatter(tau_small_neg,tau_big_neg,sz,'filled','MarkerEdgeColor','b','MarkerFaceColor','b'); hold on;
scatter(tau_small_pos,tau_big_pos,sz,'filled','MarkerEdgeColor','r','MarkerFaceColor','r');
xlabel('\tau_{small} [ms]','FontSize',x_font_sz,'FontWeight','bold');
ylabel('\tau_{large} [ms]','FontSize',y_font_sz,'FontWeight','bold');
% title(['Inhibitory \tau_{big} vs \tau_{small}  for ',get_neurons,' neurons NPSP<',num2str(NPSP_th),' ',model,' kernel']);
% annotation('textbox',[0.65 0.2 0.1 0.1],'String', sprintf('n=%0.f',n_clust),'LineStyle','none','Color','k','FontSize',all_font_sz,'FontWeight','bold');
% annotation('textbox',[0.65 0.15 0.1 0.1],'String', sprintf('inhibitory p=%0.e',p_val_neg),'LineStyle','none','Color','b','FontSize',all_font_sz,'FontWeight','bold');
% annotation('textbox',[0.65 0.1 0.1 0.1],'String', sprintf('excitatory p=%0.e',p_val_pos),'LineStyle','none','Color','r','FontSize',all_font_sz,'FontWeight','bold');
axis equal;
l = [0 lim_val_ms];
xlim(l);
ylim(l);
hline = refline(1,0);
hline.Color = 'k';
hline.LineWidth = 3;
hline.LineStyle = '--';
set(gca,'FontSize',all_font_sz,'FontWeight','Bold');
set(gcf,'color','w');
save_name = fullfile(save_dir_full,['Inhibitory and excitatory tau no npsp for ',get_neurons,' neurons with NPSP<',num2str(NPSP_th),' ',model,' kernel.png']);
export_fig(save_name);
close all;

%% Plot a histogram of inhibitory and excitatory tau_big vs tau_small as COM change ix
lim_val_ms = 60;
d_val = 2.5;
lw = 3;
%Neagtive
tau_diff_neg = tau_big_neg-tau_small_neg;
tau_sum_neg = tau_big_neg + tau_small_neg;
com_change_ix_neg = 100*(tau_diff_neg./tau_sum_neg);
edges = [-lim_val_ms:d_val:lim_val_ms];
counts_neg = histcounts(com_change_ix_neg,edges);
tau_diff_pos = tau_big_pos-tau_small_pos;
tau_sum_pos = tau_big_pos + tau_small_pos;
com_change_ix = 100*(tau_diff_pos./tau_sum_pos);
counts_pos = histcounts(com_change_ix,edges);
figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);
switch hist_type
    case 'stairs'
        fc_neg = 'none';
        fc_pos = 'none';
    case 'bar'
        fc_neg = 'b';
        fc_pos = 'r';
end

hold on;
histogram('BinEdges', edges,'BinCounts', counts_neg,'FaceColor',fc_neg,'EdgeColor','b','DisplayStyle',hist_type,'LineWidth', lw);
histogram('BinEdges', edges,'BinCounts', counts_pos,'FaceColor',fc_pos,'EdgeColor','r','DisplayStyle',hist_type,'LineWidth', lw);
xlabel('100*(\tau_{large}- \tau_{small})/(\tau_{large} + \tau_{small}) [Cetner of mass change index]','FontSize',x_font_sz,'FontWeight','bold');
ylabel('Number of neurons','FontSize',y_font_sz,'FontWeight','bold');
% title(['COM change of \tau  for ',get_neurons,' neurons NPSP<',num2str(NPSP_th),' ',model,' kernel']);
% annotation('textbox',[0.8 0.8 0.1 0.1],'String', sprintf('n=%0.f',n_clust),'LineStyle','none','Color','k','FontSize',all_font_sz,'FontWeight','bold');
% annotation('textbox',[0.8 0.75 0.1 0.1],'String', sprintf('inhibitory p=%0.e',p_val_neg),'LineStyle','none','Color','b','FontSize',all_font_sz,'FontWeight','bold');
% annotation('textbox',[0.8 0.7 0.1 0.1],'String', sprintf('excitatory p=%0.e',p_val_pos),'LineStyle','none','Color','r','FontSize',all_font_sz,'FontWeight','bold');
set(gca,'FontSize',all_font_sz,'FontWeight','Bold');
xline(0,'k','LineWidth',lw);
set(gcf,'color','w');
xlim([-lim_val_ms, lim_val_ms]);
hold off;
save_name = fullfile(save_dir_full,['Histogram ',hist_type,' of com change ix for ',get_neurons,' neurons with NPSP<',num2str(NPSP_th),' ',model,' kernel.png']);
export_fig(save_name);
close all;

%% Plot a histogram of inhibitory and excitatory tau_big vs tau_small in ms
lw = 3;
ms_spacing = 2;
lim_val_ms = 60;
edges = [-lim_val_ms:ms_spacing:lim_val_ms];
counts_neg = histcounts(tau_diff_neg,edges);
counts_pos = histcounts(tau_diff_pos,edges);
figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);
hold on;
histogram('BinEdges', edges,'BinCounts', counts_neg,'FaceColor',fc_neg,'EdgeColor','b','DisplayStyle',hist_type,'LineWidth', lw);
histogram('BinEdges', edges,'BinCounts', counts_pos,'FaceColor',fc_pos,'EdgeColor','r','DisplayStyle',hist_type,'LineWidth', lw);
xlabel('\tau_{large}- \tau_{small} [ms]','FontSize',x_font_sz,'FontWeight','bold');
ylabel('Number of neurons','FontSize',y_font_sz,'FontWeight','bold');
% title(['Change of inhibitory and excitatory \tau in ms for',get_neurons,' neurons NPSP<',num2str(NPSP_th),' ',model,' kernel']);
% annotation('textbox',[0.8 0.8 0.1 0.1],'String', sprintf('n=%0.f',n_clust),'LineStyle','none','Color','k','FontSize',axis_sz,'FontWeight','bold');
% annotation('textbox',[0.8 0.75 0.1 0.1],'String', sprintf('inhibitory p=%0.e',p_val_neg),'LineStyle','none','Color','b','FontSize',axis_sz,'FontWeight','bold');
% annotation('textbox',[0.8 0.7 0.1 0.1],'String', sprintf('excitatory p=%0.e',p_val_pos),'LineStyle','none','Color','r','FontSize',axis_sz,'FontWeight','bold');
set(gca,'FontSize',all_font_sz,'FontWeight','Bold');
xline(0,'k','LineWidth',lw);
set(gcf,'color','w');
xlim([-lim_val_ms, lim_val_ms]);
hold off;
save_name = fullfile(save_dir_full,['Histogram ',hist_type,' of tau ms change for ',get_neurons,' neurons with NPSP<',num2str(NPSP_th),' ',model,' kernel.png']);
export_fig(save_name);
close all;

%% Plot inhibitory tau vs frequency for all rooms
lw = 3;
for f = 1:length(freqs)
    curr_freq = freqs(f);
    ix = bf_closest==curr_freq;
    tau_anech_mean(f) = nanmean(tau_anech_neg(ix));
    tau_anech_std(f) = nanstd(tau_anech_neg(ix));
    tau_small_mean(f) = nanmean(tau_small_neg(ix));
    tau_small_std(f) = nanstd(tau_small_neg(ix));
    tau_big_mean(f) = nanmean(tau_big_neg(ix));
    tau_big_std(f) = nanstd(tau_big_neg(ix));
end
figure('units','normalized','outerposition',[0 0 1 1]);hold on;
plot_freq = freqs;
shadedErrorBar(plot_freq,tau_small_mean,tau_small_std./sqrt(n_clust),{'Color',[0.31 0.75 0.405]});
shadedErrorBar(plot_freq,tau_big_mean,tau_big_std./sqrt(n_clust),{'Color',[0.642 0.456 0.924]});
hold off;
xlabel('Frequency [kHz]','FontSize',x_font_sz,'FontWeight','bold');
ylabel('\tau [ms]','FontSize',y_font_sz,'FontWeight','bold');
title(['Inhibitory \tau vs frequency  for ',get_neurons,' neurons NPSP<',num2str(NPSP_th),' ',model,' kernel']);
annotation('textbox',[0.8 0.8 0.1 0.1],'String', sprintf('Big Room'),'LineStyle','none','Color',[0.642 0.456 0.924],'FontSize',all_font_sz,'FontWeight','bold');
annotation('textbox',[0.8 0.75 0.1 0.1],'String', sprintf('Small room'),'LineStyle','none','Color',[0.31 0.75 0.405],'FontSize',all_font_sz,'FontWeight','bold');
annotation('textbox',[0.8 0.7 0.1 0.1],'String', sprintf('n=%0.f',n_clust),'LineStyle','none','Color','k','FontSize',all_font_sz,'FontWeight','bold');
set(gca,'FontSize',all_font_sz,'FontWeight','Bold');
set(gcf,'color','w');
set(gca, 'XScale', 'log');
save_name = fullfile(save_dir_full,['Inhibitory tau vs freq  for ',get_neurons,' neurons NPSP<',num2str(NPSP_th),' ',model,' kernel.png']);
export_fig(save_name);
close all;

%% Plot excitatory tau vs frequency for all rooms
lw = 3;
for f = 1:length(freqs)
    curr_freq = freqs(f);
    ix = bf_closest==curr_freq;
    tau_anech_mean(f) = nanmean(tau_anech_pos(ix));
    tau_anech_std(f) = nanstd(tau_anech_pos(ix));
    tau_small_mean(f) = nanmean(tau_small_pos(ix));
    tau_small_std(f) = nanstd(tau_small_pos(ix));
    tau_big_mean(f) = nanmean(tau_big_pos(ix));
    tau_big_std(f) = nanstd(tau_big_pos(ix));
end
figure('units','normalized','outerposition',[0 0 1 1]);hold on;
plot_freq = freqs;
shadedErrorBar(plot_freq,tau_small_mean,tau_small_std./sqrt(n_clust),{'Color',[0.31 0.75 0.405]});
shadedErrorBar(plot_freq,tau_big_mean,tau_big_std./sqrt(n_clust),{'Color',[0.642 0.456 0.924]});
hold off;
xlabel('Frequency [kHz]','FontSize',x_font_sz,'FontWeight','bold');
ylabel('\tau [ms]','FontSize',y_font_sz,'FontWeight','bold');
title(['Excitatory \tau vs frequency  for ',get_neurons,' neurons NPSP<',num2str(NPSP_th),' ',model,' kernel']);
annotation('textbox',[0.8 0.8 0.1 0.1],'String', sprintf('Large room'),'LineStyle','none','Color',[0.642 0.456 0.924],'FontSize',all_font_sz,'FontWeight','bold');
annotation('textbox',[0.8 0.75 0.1 0.1],'String', sprintf('Small room'),'LineStyle','none','Color',[0.31 0.75 0.405],'FontSize',all_font_sz,'FontWeight','bold');
annotation('textbox',[0.8 0.7 0.1 0.1],'String', sprintf('n=%0.f',n_clust),'LineStyle','none','Color','k','FontSize',all_font_sz,'FontWeight','bold');
set(gca,'FontSize',all_font_sz,'FontWeight','Bold');
set(gcf,'color','w');
set(gca, 'XScale', 'log');
save_name = fullfile(save_dir_full,['Excitatory tau vs freq  for ',get_neurons,' neurons NPSP<',num2str(NPSP_th),' ',model,' kernel.png']);
export_fig(save_name);
close all;
%% Plot bf consistency accross clusters with histograms
bf_small = reach(bf,'small')/1000;
bf_big = reach(bf,'big')/1000;
diff_bf = abs(bf_big - bf_small);
figure('units','normalized','outerposition',[0 0 1 1]);
histogram(diff_bf,'Normalization','probability');
xlabel('\DeltaBF [kHz]','FontSize',x_font_sz,'FontWeight','bold');
ylabel('Probability','FontSize',y_font_sz,'FontWeight','bold');
title(['Histogram of BF consistency small<->big  for ',get_neurons,' neurons NPSP<',num2str(NPSP_th),' ',model,' kernel']);
set(gca,'FontName','Arial','FontSize',axis_sz,'FontWeight','Bold');
set(gcf,'color','w');
save_name = fullfile(save_dir_full,['Histogram of BF consistency  for ',get_neurons,' neurons NPSP<',num2str(NPSP_th),' ',model,' kernel.png']);
export_fig(save_name);
close all;
%% Plot inhibitory tau vs bf no binning small big separate
sz = 30;
tau_small = reach(tau_neg,'small');
tau_big = reach(tau_neg,'big');
figure('units','normalized','outerposition',[0 0 1 1]);
plot(bf_mean,tau_small,'.','Color',[0.31 0.75 0.405],'MarkerSize',sz);hold on;
plot(bf_mean,tau_big,'.','Color',[0.642 0.456 0.924],'MarkerSize',sz);
xlabel('BF [kHz]','FontSize',x_font_sz,'FontWeight','bold');
ylabel('\tau [ms]','FontSize',y_font_sz,'FontWeight','bold');
title(['Inhibitory \tau vs BF individual rooms  for ',get_neurons,' neurons NPSP<',num2str(NPSP_th),' ',model,' kernel']);
annotation('textbox',[0.8 0.8 0.1 0.1],'String', sprintf('Big room'),'LineStyle','none','Color',[0.642 0.456 0.924],'FontSize',all_font_sz,'FontWeight','bold');
annotation('textbox',[0.8 0.75 0.1 0.1],'String', sprintf('Small room'),'LineStyle','none','Color',[0.31 0.75 0.405],'FontSize',all_font_sz,'FontWeight','bold');
annotation('textbox',[0.8 0.7 0.1 0.1],'String', sprintf('n=%0.f',n_clust),'LineStyle','none','Color','k','FontSize',all_font_sz,'FontWeight','bold');
set(gca,'FontSize',all_font_sz,'FontWeight','Bold');
set(gca, 'XScale', 'log');
set(gcf,'color','w');
save_name = fullfile(save_dir_full,['Inhibitory tau vs BF individual rooms  for ',get_neurons,' neurons NPSP<',num2str(NPSP_th),' ',model,' kernel.png']);
export_fig(save_name);
close all;
%% Plot excitatory tau vs bf no binning small big separate
tau_small = reach(tau_pos,'small');
tau_big = reach(tau_pos,'big');
figure('units','normalized','outerposition',[0 0 1 1]);
plot(bf_mean,tau_small,'.','Color',[0.31 0.75 0.405],'MarkerSize',sz);hold on;
plot(bf_mean,tau_big,'.','Color',[0.642 0.456 0.924],'MarkerSize',sz);
xlabel('BF [kHz]','FontSize',x_font_sz,'FontWeight','bold');
ylabel('\tau [ms]','FontSize',y_font_sz,'FontWeight','bold');
title(['Excitatory \tau vs BF individual rooms  for ',get_neurons,' neurons NPSP<',num2str(NPSP_th),' ',model,' kernel']);
annotation('textbox',[0.8 0.8 0.1 0.1],'String', sprintf('Big room'),'LineStyle','none','Color',[0.642 0.456 0.924],'FontSize',all_font_sz,'FontWeight','bold');
annotation('textbox',[0.8 0.75 0.1 0.1],'String', sprintf('Small room'),'LineStyle','none','Color',[0.31 0.75 0.405],'FontSize',all_font_sz,'FontWeight','bold');
annotation('textbox',[0.8 0.7 0.1 0.1],'String', sprintf('n=%0.f',n_clust),'LineStyle','none','Color','k','FontSize',all_font_sz,'FontWeight','bold');
set(gca,'FontSize',all_font_sz,'FontWeight','Bold');
set(gca, 'XScale', 'log');
set(gcf,'color','w');
save_name = fullfile(save_dir_full,['Excitatory tau vs BF individual rooms  for ',get_neurons,' neurons NPSP<',num2str(NPSP_th),' ',model,' kernel.png']);
export_fig(save_name);
close all;