%Give the correct directory with Kilosorted files
data_root = '/media/alex/4TB SSD/DATA/Ronnie_23_01_2018/bin_files/Preprocessed/CRA_done/P06-reverb_with_noise_same_cra/P06-reverb_with_noise_same_cra_sorted';
%Give the correct Benware Grid folder
grid_root = '/media/alex/4TB SSD/DATA/Ronnie_23_01_2018/benware_gridinfo/P06-reverb_with_noise_same'; 
sound_files_root = '/media/alex/4TB SSD/DATA/Ronnie_23_01_2018/reverb_with_noise_stimuli/wav_files/reverb_with_noise_same';
sahani_path = '/media/alex/4TB SSD/DATA/Ronnie_23_01_2018/metadata/P06-reverb_with_noise_same_psth_sahani';
synch_path = '/media/alex/4TB SSD/DATA/Ronnie_23_01_2018/metadata/synch_ch_reverbwithnoisesame_P06';

kernel_type = 'insep'; %Choose either 'sep' for separable or 'insep' for inseparable
NPSP_val = 40; %The NPSP val that is used fir thresholding the 'good clusters'
time_bin_ms = 10; %The size of the time bin for the cochleagram and psth in ms
n_h = 30; %The number of history steps for the STRF


min_trig_length_s = 20; %The minmum trigger length in seconds. This will depend on the stimulus. It is best practice to check the synch channel always
fs = 30000; %Sampling rate in seconds
qualityOpt = 'nonoise'; %Which untis you want to analyse
%% Generate the cochleagram
[coch,t_s] = cochleagram_alex(sound_files_root,time_bin_ms,n_h);
%% Get the synch channel
load(synch_path);
%% Get the triggers
[start_time_ms] = get_triggers2(synch_ch,min_trig_length_s,fs);
%% Get the spike times
Y = get_spike_times(data_root,qualityOpt);
%% Load the psth_sahani and sort clusters
load(sahani_path);
count = 0;
for cluster_no = 1:numel(psth_sahani)
    if psth_sahani(cluster_no).NPSP < NPSP_val
        count = count + 1;
        good_clusters(count,1) = psth_sahani(cluster_no).cluster_id;
        good_clusters(count,2) = psth_sahani(cluster_no).NPSP;
    end
end
good_clusters = sortrows(good_clusters,2);
good_cluster_id = good_clusters(:,1);
no_good_clusters = numel(good_cluster_id);
%% Generate the PSTHs
n_h = coch(1).params.n_h;
psth_sturf = get_psth_strf(Y,start_time_ms,t_s,grid_root,good_cluster_id,n_h);
%% Generate the STRFs
strf(no_good_clusters).anech.kernel = [];
strf(no_good_clusters).reverb1.kernel = [];
strf(no_good_clusters).reverb2.kernel = [];
strf(no_good_clusters).cluster_id = [];


switch kernel_type
    case 'sep'
        fprintf('== Generating STRFs ==\n');tic;
        parfor cluster = 1:no_good_clusters
            strf(cluster).anech.kernel = sepkerneltensor2(coch(1).anech, psth_sturf(cluster).anech);
            strf(cluster).reverb1.kernel = sepkerneltensor2(coch(1).reverb1, psth_sturf(cluster).reverb1);
            strf(cluster).reverb2.kernel = sepkerneltensor2(coch(1).reverb2, psth_sturf(cluster).reverb2);
            strf(cluster).cluster_id = good_cluster_id(cluster);
        end
        fprintf('== Done! Processing took %.0fs ==\n',toc);
    case 'insep'
        fprintf('== Generating STRFs ==\n');tic;
        parfor cluster = 1:no_good_clusters
            strf(cluster).anech.kernel = elnet_fht(coch(1).anech, psth_sturf(cluster).anech,'lasso');
            strf(cluster).reverb1.kernel = elnet_fht(coch(1).reverb1, psth_sturf(cluster).reverb1,'lasso');
            strf(cluster).reverb2.kernel = elnet_fht(coch(1).reverb2, psth_sturf(cluster).reverb2,'lasso');
            strf(cluster).cluster_id = good_cluster_id(cluster);
        end
        fprintf('== Done! Processing took %.0fs ==\n',toc);
        
end

%% Plot the STRFs
plot_ix = repmat([1:no_good_clusters],1,3);
num_columns = no_good_clusters;
num_rows = 3;
per = 0.02;
edgel = 0.03; edger = per; edgeh = per; edgeb = 0.05; space_h = 0.01; space_v =0.01;
[pos]=subplot_pos(num_rows,num_columns,edgel,edger,edgeh,edgeb,space_h,space_v);

switch kernel_type
    case 'sep'
        for ii = 1:num_rows*num_columns
            subplot('position',pos{ii});
            ix = plot_ix(ii);
            if ii <= no_good_clusters
                data = strf(ix).anech.kernel.k_f*strf(ix).anech.kernel.k_h';
                imagesc(data);
                colormap(jet);
%                 plot(strf(ix).anech.kernel.k_h);
                
            elseif ii > no_good_clusters && ii <= 2*no_good_clusters
                data = strf(ix).reverb1.kernel.k_f*strf(ix).reverb1.kernel.k_h';
                imagesc(data);
                colormap(jet);
%                 plot(strf(ix).reverb1.kernel.k_h);
                
            elseif ii > 2*no_good_clusters && ii <= 3*no_good_clusters
                data = strf(ix).reverb2.kernel.k_f*strf(ix).reverb2.kernel.k_h';
                imagesc(data);
                colormap(jet);
%                 plot(strf(ix).reverb2.kernel.k_h);
            end
        end
    case 'insep'
        for ii = 1:num_rows*num_columns
            subplot('position',pos{ii});
            ix = plot_ix(ii);
            if ii <= no_good_clusters
                data = strf(ix).anech.kernel.k_fh;
                imagesc(data);
                colormap(jet);
                
            elseif ii > no_good_clusters && ii <= 2*no_good_clusters
                data = strf(ix).reverb1.kernel.k_fh;
                imagesc(data);
                colormap(jet);
       
            elseif ii > 2*no_good_clusters && ii <= 3*no_good_clusters
                data = strf(ix).reverb2.kernel.k_fh;
                imagesc(data);
                colormap(jet);
            end
        end
end