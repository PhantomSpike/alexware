close all;
no_clusters = numel(strf);

for cluster_no = 1:no_clusters
    sum_anech(cluster_no) = sum(min(strf(cluster_no).anech.kernel.k_fh(:),0).^2);
    mat_anech(:,:,cluster_no) = strf(cluster_no).anech.kernel.k_fh;
    sum_reverb1(cluster_no) = sum(min(strf(cluster_no).reverb1.kernel.k_fh(:),0).^2);
    mat_reverb1(:,:,cluster_no) = strf(cluster_no).reverb1.kernel.k_fh;
    sum_reverb2(cluster_no) = sum(min(strf(cluster_no).reverb2.kernel.k_fh(:),0).^2);
    mat_reverb2(:,:,cluster_no) = strf(cluster_no).reverb2.kernel.k_fh;
end

avg_anech = mean(min(mat_anech,0),3);
avg_reverb1 = mean(min(mat_reverb1,0),3);
avg_reverb2 = mean(min(mat_reverb2,0),3);
% avg_anech = mean(mat_anech,3);
% avg_reverb1 = mean(mat_reverb1,3);
% avg_reverb2 = mean(mat_reverb2,3);

max_val_strf = max([avg_anech(:);avg_reverb1(:);avg_reverb2(:)]);
min_val_strf = min([avg_anech(:);avg_reverb1(:);avg_reverb2(:)]);
ie_ratio_anech = sum(abs(avg_anech(avg_anech(:)<0)))/sum(avg_anech(avg_anech(:)>0));
ie_ratio_reverb1 = sum(abs(avg_reverb1(avg_reverb1(:)<0)))/sum(avg_reverb1(avg_reverb1(:)>0));
ie_ratio_reverb2 = sum(abs(avg_reverb2(avg_reverb2(:)<0)))/sum(avg_reverb2(avg_reverb2(:)>0));

avg_f_anech = mean(avg_anech);
avg_f_reverb1 = mean(avg_reverb1);
avg_f_reverb2 = mean(avg_reverb2);
max_val = max([avg_f_anech(:);avg_f_reverb1(:);avg_f_reverb2(:)]);
min_val = min([avg_f_anech(:);avg_f_reverb1(:);avg_f_reverb2(:)]);

wsize = 2;
stdev = 1;
figure;
subplot(2,3,1);
title('Anechoic Average STRF');
output = smoothts(avg_anech,'g',wsize,stdev);
mao = max(abs(output(:)));
imagesc(output,[-mao mao]);

subplot(2,3,2);
title('Little Reverb Average STRF');
output = smoothts(avg_reverb1,'g',wsize,stdev);
mao = max(abs(output(:)));
imagesc(output,[-mao mao]);

subplot(2,3,3);
title('Loads Reverb Average STRF');
output = smoothts(avg_reverb2,'g',wsize,stdev);
mao = max(abs(output(:)));
imagesc(output,[-mao mao]);

colormap('redblue');
% caxis([min_val_strf max_val_strf]);

subplot(2,3,4);
plot(avg_f_anech);
ylim([min_val max_val]);
subplot(2,3,5);
plot(avg_f_reverb1);
ylim([min_val max_val]);
subplot(2,3,6);
plot(avg_f_reverb2)
ylim([min_val max_val]);

% figure;
% hold on;
% plot(avg_f_anech);
% ylim([min_val max_val]);
% plot(avg_f_reverb1);
% ylim([min_val max_val]);
% plot(avg_f_reverb2)
% ylim([min_val max_val]);
% 
figure;
hold on;
plot((avg_f_anech-min(avg_f_anech))/(max(avg_f_anech)-min(avg_f_anech)));
plot((avg_f_reverb1-min(avg_f_reverb1))/(max(avg_f_reverb1)-min(avg_f_reverb1)));
plot((avg_f_reverb2-min(avg_f_reverb2))/(max(avg_f_reverb2)-min(avg_f_reverb2)));


% edges = [10^-6:10^-6:0.3*10^-3];
% cen = 0:0.0001:0.0025;
% figure;
% title('Anechoic');
% hist(sum_anech,cen);
% figure;
% title('Reverb1');
% hist(sum_reverb1,cen);
% figure;
% title('Reverb2');
% hist(sum_reverb2,cen);
