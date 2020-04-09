file_name{1} = '/home/alex/Desktop/Kernels/Ronnie/P13/data/sep_kernel_adjbeta';
list_n{1} = [2 3 7];
file_name{2} = '/home/alex/Desktop/Kernels/Ronnie/P06/data/sep_kernel_adjbeta';
list_n{2} = [2 3 4];
file_name{3} = '/home/alex/Desktop/Kernels/Derry/P08/data/sep_kernel_adjbeta';
list_n{3} = [3 13];
file_name{4} = '/home/alex/Desktop/Kernels/Kilkenny/P05/data/sep_kernel_adjbeta';
list_n{4} = 1;
save_dir = '/home/alex/Dropbox/SfN/';
n_h = 30;
time_bin_ms = 10;
time_spacing = 5;

count = 0;
for jj = 1:numel(file_name)
    load(file_name{jj});
    for ii = 1:numel(list_n{jj})
        count = count + 1;
        neuron = list_n{jj}(ii);
        k_h(count).anech = sep_kernel(neuron).anech_k_h.k_h;
        k_h(count).reverb1 = sep_kernel(neuron).reverb1_k_h.k_h;
        k_h(count).reverb2 = sep_kernel(neuron).reverb2_k_h.k_h;
        mao(count) = max(abs([k_h(count).anech(:);k_h(count).reverb1(:);k_h(count).reverb2(:)]));
        k_h(count).anech = k_h(count).anech./mao(count);
        k_h(count).reverb1 = k_h(count).reverb1./mao(count);
        k_h(count).reverb2 = k_h(count).reverb2./mao(count);
    end
end


%Plot the k_hs
time_ms = [-(n_h*time_bin_ms - time_bin_ms):time_spacing*time_bin_ms:0];
for jj = 1:numel(time_ms)
    x_labels{jj} = num2str(time_ms(jj),'%.0f');
end
num_columns = 3;
num_rows = 3;
per = 0.02;
edgel = 0.05; edger = per; edgeh = per; edgeb = 0.07; space_h = 0.01; space_v =0.03;
[pos]=subplot_pos(num_rows,num_columns,edgel,edger,edgeh,edgeb,space_h,space_v);
figure('units','normalized','outerposition',[0 0 1 1]);
font_y = 16;
font_x = 16;
for cluster = 1:count
    subplot('position',pos{cluster});
    hold on;
    plot(k_h(cluster).anech)
    plot(k_h(cluster).reverb1,'k');
    plot(k_h(cluster).reverb2,'r');
    ylim([-1 1]);
    hold off;
    set(findall(gca, 'Type', 'Line'),'LineWidth',3);
    hline = refline(0,0);
    hline.Color = [0.5 0.5 0.5];
    set(hline,'LineWidth',2);
    set(gca,'xtick',[]);
    if cluster == 1 || cluster == 4 || cluster == 7
        ylabel('Weight value [AU]','FontSize',font_y,'FontWeight','bold');
        set(gca,'FontName','Arial','FontSize',font_y,'FontWeight','Bold');
    else
        set(gca,'ytick',[]);
    end
    if cluster >= 7
        xticks([2:time_spacing:n_h]);
        xticklabels(x_labels);
        set(gca,'FontName','Arial','FontSize',font_x,'FontWeight','Bold');
        xlabel('History [ms]','FontSize',font_x,'FontWeight','bold');
%         legend('Anechoic Room','Small Room','Big Room','Location','northwest');
    end
    if cluster == 7
        legend('Anechoic Room','Small Room','Big Room','Location','northwest');
    end
end
set(gcf,'color','w');
save_name = [save_dir,'k_hs.jpg'];
export_fig(save_name);