function plot_fra_pixels(fra_psth,clusters,save_dir)

%% Define plotting parameters
row = 5;
col = 4;
per = 0.005;
edgel = 0.035; edger = per; edgeh = per; edgeb = 0.05; space_h = per; space_v = 0.01;
[pos]=subplot_pos(row,col,edgel,edger,edgeh,edgeb,space_h,space_v);
num_clust = numel(clusters);
freqs = fra_psth(1).params.freqs;
num_freqs = numel(freqs);
for jj = 1:numel(freqs)
    x_labels{jj} = num2str(freqs(jj),'%.1f');
end
dB_lvls = fra_psth(1).params.dB_levels;
num_dB_lvls = numel(dB_lvls);
%% Plot the correct FRAs
clust_spacing = row*col - 1;
num_groups = ceil(num_clust/(clust_spacing+1));
last_spacing = num_clust - (num_groups-1)*(clust_spacing+1) - 1;
first_clust_ix = [1:clust_spacing+1:(num_groups)*(clust_spacing+1)];


for group = 1:num_groups
    
    first_clust = first_clust_ix(group);
    last_clust = first_clust + clust_spacing;
    
    if group == num_groups
        last_clust = first_clust + last_spacing;
    end
    num_subplots = numel(first_clust:last_clust);
    plot_ix = [first_clust:last_clust];
    
    figure('units','normalized','outerposition',[0 0 1 1]);
    for ii = 1:num_subplots
        subplot('position',pos{ii});
        ix_clust = plot_ix(ii);
        clust = clusters(ix_clust);
        ix = find(fra_psth(1).cluster_id == clust);
        imagesc(fra_psth(ix).X_dbft);
        colormap('jet');
        if ii == (row - 1)*col + 1
            xticks([1:num_freqs]);
            xticklabels(x_labels);
            yticks(1:num_dB_lvls);
            y_labels = string(dB_lvls);
            yticklabels(y_labels);
            set(gca,'FontName','Arial','FontSize',8,'FontWeight','Bold');
            xlabel('Frequency [kHz]','FontSize',10,'FontWeight','bold');
            ylabel('Sound Level [dB]','FontSize',16,'FontWeight','bold');
        elseif ii > (row - 1)*col + 1 && ii <= row*col
            xticks([1:num_freqs]);
            xticklabels(x_labels);
            set(gca,'FontName','Arial','FontSize',8,'FontWeight','Bold');
            xlabel('Frequency [kHz]','FontSize',10,'FontWeight','bold');
            set(gca,'ytick',[]);
        else
            set(gca,'ytick',[]);
            set(gca,'xtick',[]);
        end
    end
    if exist('save_dir','var')
        file_name = ['FRAs_Clusters_',num2str(clusters(plot_ix(1))),'_to_',num2str(clusters(plot_ix(end)))];
        save_name = [save_dir,file_name,'.png'];
        export_fig(save_name);
        close;
    end
end