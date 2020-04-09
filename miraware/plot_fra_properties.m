function plot_fra_properties(fra_psth,type,save_dir)

%% Define plotting parameters
switch type
    case 'exc'
        fras = {fra_psth.fra_prop_exc}';
        fras = {fras{~cellfun(@isempty,fras)}}';
    case 'inh'
        fras = {fra_psth.fra_prop_inh}';
        fras = {fras{~cellfun(@isempty,fras)}}';
end

row = 5;
col = 4;
per = 0.005;
edgel = 0.035; edger = per; edgeh = per; edgeb = 0.05; space_h = per; space_v = 0.01;
[pos]=subplot_pos(row,col,edgel,edger,edgeh,edgeb,space_h,space_v);
num_clust = numel(fras);
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

cl = 0;
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
        cl = cl+1;
        subplot('position',pos{ii});
        imagesc(fras{cl}.fra_smooth);
        colormap('inferno');
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
        %Plot the threshold, Q10, Q30 and CF lines
        hold on;
        line(fras{cl}.f_plot_ix,fras{cl}.lvl_plot_ix,'LineWidth',2,'Color', 'w');
        if ~isempty(fras{cl}.q10)
            line(fras{cl}.q10_f_ix,fras{cl}.q10_lvl_ix,'LineWidth',2,'Color', 'k');
        end
        if ~isempty(fras{cl}.q30)
            line(fras{cl}.q30_f_ix,fras{cl}.q30_lvl_ix,'LineWidth',2,'Color', 'k');
        end
        
        if ~isempty(fras{cl}.cf)
            xline(fras{cl}.cf_ix,'-',{[num2str(fras{cl}.cf,'%.1f'),'kHz']},'LineWidth',2,'Color', 'w');
        end
    end
    if exist('save_dir','var')
        file_name = ['FRAs_Clusters_batch',num2str(group),];
        save_name = fullfile(save_dir,[file_name,'.png']);
        export_fig(save_name);
        close;
    end
end