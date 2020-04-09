cluster_id = 455+1;
t_bin_s = 0.02;
edges = [0:t_bin_s:36];
for stim = 1:6
    for rep = 1:10
        resp(stim,rep,:) = histc(data(cluster_id).stim(stim).repeat(rep).spiketimes,edges);
    end
end

for j = 1:6
    subplot(3,2,j);
    imagesc(squeeze(resp(j,:,:)));
end