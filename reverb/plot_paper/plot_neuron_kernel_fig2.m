clear all; close all;
%% Params
model = 'ridge';
animal_names{1} = 'Ronnie'; animal_names{2} = 'Ronnie'; animal_names{3} = 'Ronnie';  
pen_names{1} = 'P06'; pen_names{2} = 'P06'; pen_names{3} = 'P13'; 
clusters{1} = '213'; clusters{2} = '221'; clusters{3} = '528'; 
kernel_dir = '/mnt/40086D4C086D41D0/Reverb_neuronal_data/Kernel_fits/kfolds/perfreq/ridge/10ms/200ms';
n_ker = length(animal_names);
rooms{1} = 'big';
rooms{2} = 'small';
save_dir = '/mnt/40086D4C086D41D0/Reverb_paper/fig_2';
%% Plot
font_type = 'Liberation Sans';
sz = 25;
y_font_sz = sz;
x_font_sz = sz;
all_font_sz = sz;

h_lim_ms = 210;
row = 2;
col = 3;
lw = 3;
per = 0.01;
edgel = 0.07; edger = per; edgeh = 0.03; edgeb = 0.1; space_h = 0.025; space_v = 0.03;
[pos]=subplot_pos(row,col,edgel,edger,edgeh,edgeb,space_h,space_v);

kernels = load(fullfile(kernel_dir,strjoin({animal_names{1},pen_names{1},clusters{1}},'_')));
freqs = fliplr(kernels.kernel.freqs);
n_f = length(freqs);
n_h = kernels.kernel.n_h;
dt = round(kernels.kernel.dt_ms);
h_max_ms = kernels.kernel.h_max_ms;
h_steps = [0:dt:h_max_ms-dt];
h_lim = round(h_lim_ms/dt);
skip_f = 5;
skip_h = 5;

count = 0;
for f = 1:skip_f:n_f
    count = count+1;
    f_labels{count} = num2str(freqs(f)./1000,'%.1f');
end

count = 0;
for t = 1:skip_h:n_h
    count = count+1;
    h_labels{count} = num2str(h_steps(t),'%.0f');
end



figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);
for k = 1:n_ker
    kernels = load(fullfile(kernel_dir,strjoin({animal_names{k},pen_names{k},clusters{k}},'_')));
    for r = 1:2
        room = rooms{r};
        switch model
            case {'sep','sep_kh'}
                k_fh = kernels.kernel.(room).k_f*fliplr(kernels.kernel.(room).k_h');
                k_h.(room) = fliplr(kernels.kernel.(room).k_h');
                k_h.(room) = k_h.(room)/max(abs(k_h.(room)));
            case {'lasso','ridge'}
                k_fh = fliplr(kernels.kernel.(room).main{end}.k_fh);
        end
        subplot('position',pos{(r-1)*col + k});
        k_fh = k_fh./max(abs(k_fh(:)));
        imagesc(k_fh);
        caxis([-1 1]);
        colormap('redblue');
        set(gca,'XTick',[], 'YTick', []);
        if k == 1 || k == 4
            yticks([1:skip_f:n_f]);
            yticklabels(f_labels);
            ylabel('Freqeuncy [kHz]','FontSize',y_font_sz,'FontWeight','bold');
        end
        
        if r == 2
            xticks([1:skip_h:h_max_ms+dt]);
            xticklabels(h_labels);
            xlabel('History [ms]','FontSize',x_font_sz,'FontWeight','bold');
        end
        set(gcf,'color','w');
        set(gca,'FontSize',all_font_sz,'FontWeight','Bold');
        xlim([0 h_lim]);
        
    end
end

save_name = fullfile(save_dir,['Neuronal_kernels_STRF_examples.png']);
export_fig(save_name);
close all;

figure('units','normalized','outerposition',[0 0 1 1],'DefaultTextFontName', font_type, 'DefaultAxesFontName',font_type);
for k = 1:n_ker
    kernels = load(fullfile(kernel_dir,strjoin({animal_names{k},pen_names{k},clusters{k}},'_')));
    for r = 1:2
        room = rooms{r};
        switch model
            case {'sep','sep_kh'}
                k_fh = kernels.kernel.(room).k_f*fliplr(kernels.kernel.(room).k_h');
                k_h.(room) = fliplr(kernels.kernel.(room).k_h');
                k_h.(room) = k_h.(room)/max(abs(k_h.(room)));
            case {'lasso','ridge'}
                k_fh = fliplr(kernels.kernel.(room).main{end}.k_fh);
                k_fh_neg = abs(min(k_fh,0));
                k_fh_pos = abs(max(k_fh,0));
                k_h_neg.(room) = mean(k_fh_neg);
                k_h_pos.(room) = mean(k_fh_pos);
                k_h_neg_f.(room) = k_h_neg.(room)./sum(k_h_neg.(room)(:)); %Scale the values to sum to 1 for inhibition
                k_h_pos_f.(room) = k_h_pos.(room)./sum(k_h_pos.(room)(:)); %Scale the values to sum to 1 for inhibition
        end
        tau_neg.(room) = (k_h_neg_f.(room)*h_steps'); %Compute a weighted sum of all values
        tau_pos.(room) = (k_h_pos_f.(room)*h_steps'); %Compute a weighted sum of all values
    end
    %Inhibition
    subplot('position',pos{k});
    set(gca,'XTick',[]);
    hold on;
    max_val = max([k_h_neg.big(:);k_h_neg.small(:)]);
    plot([0 -k_h_neg.big./max_val],'Color','b','LineWidth',lw);
    plot([0 -k_h_neg.small./max_val],'Color',[0.0588 1 1],'LineWidth',lw);
    hold off;
    if k == 1
        ylabel('Weight [AU]','FontSize',y_font_sz,'FontWeight','bold');
    else
        set(gca,'YTick', []);
    end
    set(gca,'FontSize',all_font_sz,'FontWeight','Bold');
    legend(['\tau_{large}=',num2str(tau_neg.big,'%0.f'),'ms'],['\tau_{small}=',num2str(tau_neg.small,'%0.f'),'ms'],'Location','southeast')
    set(gcf,'color','w');
    xlim([0 h_lim]);
    %Excitation
    subplot('position',pos{col + k});
    set(gca,'XTick',[]);
    hold on;
    max_val = max([k_h_pos.big(:);k_h_pos.small(:)]);
    plot([0 k_h_pos.big./max_val],'Color','r','LineWidth',lw);
    plot([0 k_h_pos.small./max_val],'Color',[1 0.7 0.7],'LineWidth',lw);
    hold off;
    if k == 1
        ylabel('Weigth [AU]');
    else
        set(gca,'YTick', []);
    end
    
    if r == 2
        xticks([1:skip_h:h_max_ms]);
        xticklabels(h_labels);
        xlabel('History [ms]','FontSize',x_font_sz,'FontWeight','bold');
    end
    
    set(gca,'FontSize',all_font_sz,'FontWeight','Bold');
    legend(['\tau_{large}=',num2str(tau_pos.big,'%0.f'),'ms'],['\tau_{small}=',num2str(tau_pos.small,'%0.f'),'ms'])
    set(gcf,'color','w');
    xlim([0 h_lim]);
end

save_name = fullfile(save_dir,['Neuronal_kernels_EI_examples.png']);
export_fig(save_name);
close all;