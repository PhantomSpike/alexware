function [strf,psth_sturf,coch] = get_strf2(kernel_type,plot,store)

%Give the correct directory with Kilosorted files
data_root = ['/mnt/40086D4C086D41D0/Reverb_data/Derry/P08_Derry/P08-alex_reverb_wtih_noise_same.3/CRA/'];
%Give the correct Benware Grid folder
grid_root = ['/mnt/40086D4C086D41D0/Reverb_data/Derry/benware_gridinfo/P08-alex_reverb_wtih_noise_same.3/'];
sound_files_root = '/mnt/40086D4C086D41D0/Reverb_data/Derry_Kilkenny_Cork_stims/reverb_with_noise_stimuli/wav_files';
sahani_path = ['/mnt/40086D4C086D41D0/Reverb_data/Derry/metadata/P08-reverb_with_noise_same_psth_sahani'];
synch_path = ['/mnt/40086D4C086D41D0/Reverb_data/Derry/metadata/synch_ch_reverbwithnoisesame_P08'];


NPSP_val = 40; %The NPSP val that is used fir thresholding the 'good clusters'
time_bin_ms = 10; %The size of the time bin for the cochleagram and psth in ms
n_h = 30; %The number of history steps for the STRF

% folder_name = [pen,'/',kernel_type,'/'];


min_trig_length_s = 40; %The minmum trigger length in seconds. This will depend on the stimulus. It is best practice to check the synch channel always
min_inter_trig_length_s = 0.1;
fs = 30000; %Sampling rate in seconds
qualityOpt = 'nonoise'; %Which untis you want to analyse
%% Generate the cochleagram
[coch,t_s] = cochleagram_alex(sound_files_root,time_bin_ms,n_h);
%% Get the synch channel
load(synch_path);
%% Get the triggers
[start_time_ms] = get_triggers_new(synch_ch,min_trig_length_s,min_inter_trig_length_s,fs);
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
if plot
    gen_name = '/home/alex/Desktop/Kernels/';
    clust_spacing = 9;
    num_groups = ceil(no_good_clusters/(clust_spacing+1));
    last_spacing = no_good_clusters - (num_groups-1)*(clust_spacing+1) - 1;
    first_clust_ix = [1:clust_spacing+1:(num_groups)*(clust_spacing+1)];
    for group = 1:num_groups
        
        first_clust = first_clust_ix(group);
        last_clust = first_clust + clust_spacing;
        if group == num_groups
            last_clust = first_clust + last_spacing;
        end
        num_clust = numel(first_clust:last_clust);
        plot_ix = repmat([first_clust:last_clust],1,3);
        num_columns = num_clust;
        num_rows = 3;
        per = 0.02;
        edgel = 0.04; edger = per; edgeh = per; edgeb = 0.05; space_h = 0.01; space_v =0.01;
        [pos]=subplot_pos(num_rows,num_columns,edgel,edger,edgeh,edgeb,space_h,space_v);
        
        spacing_freq = 2;
        spacing_time = 6;
        freqs = coch(1).params.f;
        num_freqs = numel(freqs);
        freqs = ceil(freqs)/1000; %Convert the frequencies into kHz
        freqs = freqs([1:spacing_freq:num_freqs]);
        time_ms = [-(n_h*time_bin_ms - time_bin_ms):spacing_time*time_bin_ms:0];
        
        for jj = 1:numel(time_ms)
            x_labels{jj} = num2str(time_ms(jj),'%.0f');
        end
        
        for ii = 1:numel(freqs)
            y_labels{ii} = num2str(freqs(ii),'%.1f');
        end
        
        
        switch kernel_type
            case 'sep'
                if store
                    figure('units','normalized','outerposition',[0 0 1 1]);
                else
                    figure;
                end
                for ii = 1:num_rows*num_columns
                    subplot('position',pos{ii});
                    ix = plot_ix(ii);
                    if ii <= num_clust
                        data = strf(ix).anech.kernel.k_f*strf(ix).anech.kernel.k_h';
                        mao = max(abs(data(:)));
                        imagesc(data,[-mao mao]);
                        colormap('redblue');
                        set(gca,'xtick',[]);
                        set(gca,'ytick',[]);
                        
                        if ii == 1
                            yticks([1:spacing_freq:num_freqs]);
                            yticklabels(y_labels);
                            set(gca,'FontName','Arial','FontSize',10,'FontWeight','Bold');
                            ylabel('Frequency [kHz]','FontSize',16,'FontWeight','bold');
                        end
                        
                    elseif ii > num_clust && ii <= 2*num_clust
                        data = strf(ix).reverb1.kernel.k_f*strf(ix).reverb1.kernel.k_h';
                        mao = max(abs(data(:)));
                        imagesc(data,[-mao mao]);
                        colormap('redblue');
                        set(gca,'xtick',[]);
                        set(gca,'ytick',[]);
                        
                        if ii == num_clust + 1
                            yticks([1:spacing_freq:num_freqs]);
                            yticklabels(y_labels);
                            set(gca,'FontName','Arial','FontSize',10,'FontWeight','Bold');
                            ylabel('Frequency [kHz]','FontSize',16,'FontWeight','bold');
                        end
                        
                    elseif ii > 2*num_clust && ii <= 3*num_clust
                        data = strf(ix).reverb2.kernel.k_f*strf(ix).reverb2.kernel.k_h';
                        mao = max(abs(data(:)));
                        imagesc(data,[-mao mao]);
                        colormap('redblue');
                        colormap('redblue');
                        xticks([1:spacing_time:n_h]);
                        xticklabels(x_labels);
                        set(gca,'ytick',[]);
                        set(gca,'FontName','Arial','FontSize',10,'FontWeight','Bold');
                        xlabel('History [ms]','FontSize',14,'FontWeight','bold');
                        
                        if ii == 2*num_clust + 1
                            yticks([1:spacing_freq:num_freqs]);
                            yticklabels(y_labels);
                            ylabel('Frequency [kHz]','FontSize',16,'FontWeight','bold');
                        end
                        
                    end
                end
                
            case 'insep'
                if store
                    figure('units','normalized','outerposition',[0 0 1 1]);
                else
                    figure;
                end
                for ii = 1:num_rows*num_columns
                    subplot('position',pos{ii});
                    ix = plot_ix(ii);
                    if ii <= num_clust
                        data = strf(ix).anech.kernel.k_fh;
                        mao = max(abs(data(:)));
                        imagesc(data,[-mao mao]);
                        colormap('redblue');
                        set(gca,'xtick',[]);
                        set(gca,'ytick',[]);
                        
                        if ii == 1
                            yticks([1:spacing_freq:num_freqs]);
                            yticklabels(y_labels);
                            set(gca,'FontName','Arial','FontSize',10,'FontWeight','Bold');
                            ylabel('Frequency [kHz]','FontSize',16,'FontWeight','bold');
                        end
                        
                    elseif ii > num_clust && ii <= 2*num_clust
                        data = strf(ix).reverb1.kernel.k_fh;
                        mao = max(abs(data(:)));
                        imagesc(data,[-mao mao]);
                        colormap('redblue');
                        set(gca,'xtick',[]);
                        set(gca,'ytick',[]);
                        
                        if ii == num_clust + 1
                            yticks([1:spacing_freq:num_freqs]);
                            yticklabels(y_labels);
                            set(gca,'FontName','Arial','FontSize',10,'FontWeight','Bold');
                            ylabel('Frequency [kHz]','FontSize',16,'FontWeight','bold');
                        end
                        
                    elseif ii > 2*num_clust && ii <= 3*num_clust
                        data = strf(ix).reverb2.kernel.k_fh;
                        mao = max(abs(data(:)));
                        if mao == 0
                            mao = 1; %Hack to fix empty STRFs
                        end
                        imagesc(data,[-mao mao]);
                        colormap('redblue');
                        xticks([1:spacing_time:n_h]);
                        xticklabels(x_labels);
                        set(gca,'ytick',[]);
                        set(gca,'FontName','Arial','FontSize',10,'FontWeight','Bold');
                        xlabel('History [ms]','FontSize',14,'FontWeight','bold');
                        
                        if ii == 2*num_clust + 1
                            yticks([1:spacing_freq:num_freqs]);
                            yticklabels(y_labels);
                            ylabel('Frequency [kHz]','FontSize',16,'FontWeight','bold');
                        end
                        
                    end
                end
        end
        file_name = [gen_name,'Clusters_',num2str(first_clust),'_',num2str(last_clust),'.png'];
        if store
            export_fig(file_name);
            close;
        end
    end
end
end