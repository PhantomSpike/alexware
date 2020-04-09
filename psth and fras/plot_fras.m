function plot_fras(fra)

cluster_list = [1:100];
per = 0.02;
edgel = 0.03; edger = per; edgeh = per; edgeb = 0.05; space_h = 0.01; space_v =0.01;
[pos]=subplot_pos(10,10,edgel,edger,edgeh,edgeb,space_h,space_v);
% norm_fra = fras(1:5,:,:,11:21) - fras(1:5,:,:,31:41);
% norm_fra = norm_fra - mean(fras(1:5,:,:,1:10),4);
for cluster = 1:100
    subplot('position',pos{cluster});
    %     norm_fra = sum(fras(1:5,:,cluser_list(cluster),1:10) - fras(1:5,:,cluser_list(cluster),20:29));
    %     imagesc(flipud(norm_fra),4);
    cluster_no = cluster_list(cluster);
    imagesc(fra(:,:,cluster_no));
    %     imagesc(flipud(sum(fras(1:5,:,cluser_list(cluster),2:10),4)));
    %     max_resp = max(fras(1:5,:,cluser_list(cluster),1:10));
    %     imagesc(max_resp(1:5,:,cluser_list(cluster),1),4);
    colormap('jet');
    if cluster == 91
        xlabel('Frequency');
        ylabel('dB Level');
        set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold');
    end
end