function plot_betas_model(pen)
dir_name = '/home/alex/Desktop/Kernels/';
full_name = [dir_name,pen,'/data/model'];
load(full_name);
num_clusters = numel(model_kernel);
for clust = 1:num_clusters
    beta_small(clust) = model_kernel(clust).dexp_reverb1.beta;
    beta_big(clust) = model_kernel(clust).dexp_reverb2.beta;
    beta_anech(clust) = model_kernel(clust).dexp_anech.beta;
    NPSP(clust) = model_kernel(clust).NPSP;
end
sz = 30;
ratio_big_small = beta_big./beta_small;
figure;
plot(NPSP,beta_anech,'b.','MarkerSize',sz);
hold on;
plot(NPSP,beta_small,'k.','MarkerSize',sz);
plot(NPSP,beta_big,'r.','MarkerSize',sz);
ylim([0 200]);
xlim([0 40]);
figure;
plot(NPSP,ratio_big_small,'.','MarkerSize',sz);
hold on;
plot(NPSP,ones(numel(NPSP),1),'k--');
ylim([0 5]);
xlim([0 40]);