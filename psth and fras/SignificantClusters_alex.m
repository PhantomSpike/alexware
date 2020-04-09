clear all; close all;
%Load the spike times and cluster ids + the file containing information
%about the stimulus which was used
data_path = '/home/phant0msp1ke/Desktop/Data/raw_ephys/SC_Mouse_30_06_17/Orientation/Orientation_car_filt_100_1900_sorted';
stim_file = '/home/phant0msp1ke/Desktop/Code/ForAlex/StimData/2017_06_30/2017_06_30_11_58_35_classic_ori_trig_gray_gaps_NewTrigger_SmartProbes.mat';
fs = 30000;
start_time_s = 0;
end_time_s = 1800;
chunk_s = 100;
synch_ch = load_synch(data_path,start_time_s,end_time_s,chunk_s);

%Find the difference between every n+1 sample - n sample. This tells us the
%the beginning/end of each stim presentation
diff_sig = diff(synch_ch); 

%Find the time in seconds for the beginning of each stim presentation
start_ix = find(diff_sig==1);
stim_start_s = start_ix./fs;

spikes = get_spike_times(data_path,'nonoise'); %Get only the spike times which are not classified as noise

spike_times_s = spikes(:,1); %Spike times in seconds
clusters = spikes(:,2);
cluster_id = unique(clusters); %Sort the clusters which are good in ascending order
total_no_clusters = length(cluster_id); %The total number of unqiue clusters

% load stimulus info
s = load(stim_file);
stiminfo = s.rand_log;

oris = [];

for rr = 1:length(stiminfo.repetition)
    oris = [oris stiminfo.repetition(rr).eye.ori];
end

unique_oris = unique(oris);
t = -2:.05:3;
%%
Sig(1:total_no_clusters)=NaN;
OSI_clusters(1:total_no_clusters)=NaN;
DSI_clusters(1:total_no_clusters)=NaN;
response_strength_reps=[];
response_strength_reps2=[];

% spike times for cluster whose number is cluster_id, in samples
%cluster_no is the counting variable and cluster_id is the actual cluster
%id
for cluster_no = 1:total_no_clusters
    cluster_times_sec = spike_times_s(clusters==cluster_id(cluster_no));
    n_spikes = length(cluster_times_sec);
    fprintf('Total number of spikes = %d\n', n_spikes);
    
    psthes = zeros(length(t), length(unique_oris));
    response_strength = zeros(1, length(unique_oris));
    response_strength_reps=[];
    response_strength_reps2=[];
    
    for oo = 1:length(oris)
        ori = oris(oo);
        start_time = stim_start_s(oo);
        
        ori_idx = find(unique_oris==ori);
        bins = start_time + t;
        %fprintf('%d %0.2f %0.2f\n', ori, min(bins), max(bins));
        psth = histc(cluster_times_sec, bins);
        psthes(:, ori_idx) = psthes(:, ori_idx) + psth;
        psthes_reps(:, oo) =  psth;
        
        
        
        response_strength(ori_idx) = response_strength(ori_idx) + sum(cluster_times_sec>start_time-1 & cluster_times_sec<(start_time+1));
        response_strength_reps(oo) =  sum(cluster_times_sec>(start_time-1 ) & cluster_times_sec<(start_time+1));
        
    end
    [o oo]=sort(oris);
    psthes_reps=permute(reshape(psthes_reps(:,oo),[101],s.n_repetitions,8),[1,3,2]); % [bins, ori, repetitions]
    response_strength_reps2=reshape(response_strength_reps(oo),s.n_repetitions,8)'; % ori x repetitions
    % Calculate significance
    Sig(cluster_no)=anova1(response_strength_reps2',[],'off');
    
    [OSI, DSI, PrefDir, OppDir, OrthoDir1, OrthoDir2]=OsiDsi(mean(response_strength_reps2,2),unique_oris);
    OSI_clusters(cluster_no)=OSI;
    DSI_clusters(cluster_no)=DSI;
    response_strength_reps3(cluster_no,:,:)=response_strength_reps2;
    
    psthes2(cluster_no,:,:)=psthes;
    
    %Shows you different clusters and pauses, so you can select them
%         figure(1);
%         for ii = 1:length(unique_oris)
%             subplot(4, 2, ii);
%             bar(t, psthes(:, ii));
%             axis tight;
%             title(unique_oris(ii));
%             ylim([0, 40]);
%         end
%         figure(2);
%         th = deg2rad(unique_oris);
%         %r = sum(psthes, 1);
%         r = response_strength;
%         polarplot([th th(1)], [r r(1)]);
%         pause;
end
display([ 'Number of clusters: ' num2str(sum(~isnan(Sig)))])
display([ 'Number of significant clusters: ' num2str(sum(Sig(~isnan(Sig))<0.05))])
display([ 'Percentage of significant clusters: ' num2str(sum(Sig(~isnan(Sig))<0.05)./sum(~isnan(Sig)).*100)])
% distribution of OSI and DSI
% Repeat for all recordings, remember to change the number of repetitions and carefull when the recording failed
% check timings!!!
figure, hist(OSI_clusters(Sig<0.05))
sig_clusters_ix = find(Sig<0.05);
sig_clusters_id = cluster_id(sig_clusters_ix);
%% Plots
%Tell which clusters you want for plotting orientation tunning curves
cluster_ids_examples= [10,50,160];
Example_tuning_curves(response_strength_reps3,psthes2,unique_oris,cluster_ids_examples);