%% Define and load vars

% LNP kernels
cluster_dir = '/mnt/40086D4C086D41D0/Reverb_neuronal_data/Kernel_fits/kfolds/LNP_model/perfreq/ridge/10ms/200ms'; %Absolute path to the kernel used for the prediction
%Get the cochleagrams
coch_file = '/mnt/40086D4C086D41D0/Reverb_neuronal_data/Cochleagrams/Cochleagrams/New_th/Derry_Kilkenny_Cork/specpower/10ms/coch_all_conditions_specpower.mat';
load(coch_file,'coch');

save_dir = '/mnt/40086D4C086D41D0/Reverb_neuronal_data/Kernel_fits/kfolds/LNP_model/perfreq/ridge/10ms/Predicted_PSTHs';

NPSP_th = 40;
n_h = 20; %Number of history steps
n_trials = 10; %How many trials to simulate
%% Select the clusters to fit
load(fullfile(cluster_dir,'info.mat'),'info'); %Load the info file with the meta data
ix = info.NPSP<NPSP_th; %Find only the neurons that are below the NPSP th
%Get the necessary fields
cluster_ids = info.cluster_id(ix);
pen_names = info.pen_name(ix);
animal_names = info.animal_name(ix);
qualities = info.quality(ix);
NPSPs = info.NPSP(ix);
n_clusters = sum(ix);


%% Normalize the cochleagrams
small1_X_ft = coch(3).X_ft;
small2_X_ft = coch(4).X_ft;
big1_X_ft = coch(5).X_ft;
big2_X_ft = coch(6).X_ft;

mean_reverb = mean([small1_X_ft, small2_X_ft, big1_X_ft, big2_X_ft],2); std_reverb = std([small1_X_ft, small2_X_ft, big1_X_ft, big2_X_ft],[],2);
        
small1_X_ft_norm = (small1_X_ft - mean_reverb)./std_reverb;
small2_X_ft_norm = (small2_X_ft - mean_reverb)./std_reverb;
big1_X_ft_norm = (big1_X_ft - mean_reverb)./std_reverb;
big2_X_ft_norm = (big2_X_ft - mean_reverb)./std_reverb;


%% Tensorize
fprintf('== Tensorizing the cochleagrams ==\n');tic;
small1_X_fht = tensorize(small1_X_ft_norm,n_h);
small2_X_fht = tensorize(small2_X_ft_norm,n_h);
big1_X_fht = tensorize(big1_X_ft_norm,n_h);
big2_X_fht = tensorize(big2_X_ft_norm,n_h);
fprintf('== Done! This took %0.fs ==\n',toc);

%% Concatenate together
small_X_fht = cat(3,small1_X_fht,small2_X_fht);
big_X_fht = cat(3,big1_X_fht,big2_X_fht);
n_t = size(small_X_fht,3);

%% Make simulated firing rate 
tic;
parfor c = 1:n_clusters
    fprintf('== Processing cluster %0.f/%0.f ==\n',c,n_clusters);
    %Initialise vars
    y_t_small_pois = [];
    y_t_big_pois = [];
    temp = [];

    
    %Get the kernel for every cluster
    temp = load(fullfile(cluster_dir,strjoin({animal_names{c},pen_names{c},num2str(cluster_ids(c))},'_')));
    ker = temp.kernel.reverb.main{end};
    
    %Make the linear prediction
    y_t_small = kernelconv(small_X_fht, ker);
    y_t_big = kernelconv(big_X_fht, ker);
    
    %Make all -ve values 0
    y_t_small = max(y_t_small,0);
    y_t_big = max(y_t_big,0);
    
    %Add Poisson noise
    for tr = 1:n_trials
        y_t_small_pois(tr,:) = poissrnd(y_t_small);
        y_t_big_pois(tr,:) = poissrnd(y_t_big);
    end
    
    %Average the virtual trials
    temp_psth{c}.y_t_small_mean = mean(y_t_small_pois);
    temp_psth{c}.y_t_big_mean = mean(y_t_big_pois);

end


%Save the results
fprintf('== Saving the results ==\n');
for cl = 1:n_clusters
    psth.y_t_small_mean = temp_psth{cl}.y_t_small_mean;
    psth.y_t_big_mean = temp_psth{cl}.y_t_small_mean;
    save(fullfile(save_dir,strjoin({animal_names{cl},pen_names{cl},num2str(cluster_ids(cl))},'_')),'psth');
end

clear info;
info.cluster_id = cluster_ids;
info.NPSP = NPSPs;
info.animal_name = animal_names;
info.quality = qualities;
info.pen_name = pen_names;
save(fullfile(save_dir,'info'),'info');

fprintf('== Done! This took %0.fs ==\n',toc);
