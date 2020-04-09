function model_kernel = model_neuron(animal,pen,model,method,plot_dexp,plot_com,plot_res,store)

%% Get the seperable kernel using all the data
dir_name = ['/home/alex/Desktop/Kernels/',animal,'/'];
data_dir = ['/mnt/40086D4C086D41D0/Reverb_data/',animal,'/'];
meta_dir = ['/mnt/40086D4C086D41D0/Reverb_data/',animal,'/metadata/'];
grid_dir = ['/mnt/40086D4C086D41D0/Reverb_data/',animal,'/benware/'];
%Give the correct directory with Kilosorted files
data_root = [data_dir,pen,'_',animal,'/',pen,'-alex_reverb_wtih_noise_same/CRA/'];
%Give the correct Benware Grid folder
grid_path = [grid_dir,pen,'-alex_reverb_with_noise_same'];
sound_files_root = '/mnt/40086D4C086D41D0/Reverb_data/Derry_Kilkenny_Cork_stims/reverb_with_noise_stimuli/wav_files';
sahani_path = [meta_dir,pen,'-reverb_with_noise_same_psth_sahani'];
synch_path = [meta_dir,'synch_ch_reverbwithnoisesame_',pen];


NPSP_val = 40; %The NPSP val that is used fir thresholding the 'good clusters'
time_bin_ms = 10; %The size of the time bin for the cochleagram and psth in ms
n_h = 30; %The number of history steps for the STRF
time_spacing = 2;
freq_spacing = 2;
betainit_coeff = 6;
alphainit_coeff = 1;
extra_per_A = 0;
extra_per_B = 0;
lim_val = 140;


min_trig_length_s = 40; %The minmum trigger length in seconds. This will depend on the stimulus. It is best practice to check the synch channel always
min_inter_trig_length_s = 0.1;
fs = 30000; %Sampling rate in seconds
qualityOpt = 'nonoise'; %Which untis you want to analyse
%% Generate the cochleagram
[coch,t_s] = cochleagram_alex(sound_files_root,time_bin_ms,n_h);
%% Get the synch channel
load(synch_path);
%% Get the triggers
if strcmp(animal,'Ronnie') && strcmp(pen,'P06')
    [start_time_ms] = get_triggers2(synch_ch,min_trig_length_s,fs);
else
    [start_time_ms] = get_triggers_new(synch_ch,min_trig_length_s,min_inter_trig_length_s,fs);
end
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
no_clusters = numel(good_cluster_id);
%% Generate the PSTHs
n_h = coch(1).params.n_h;
psth_sturf = get_psth_strf(Y,start_time_ms,t_s,grid_path,good_cluster_id,n_h);

%% Generate the model responses
coch_total = cat(3,coch(1).anech,coch(1).reverb1,coch(1).reverb2); %Concatenate all the cochleagrams together

model_kernel(no_clusters).full_kernel_real = [];
model_kernel(no_clusters).full_kernel_model = [];
model_kernel(no_clusters).k_h_model_anech = [];
model_kernel(no_clusters).k_h_model_reverb1 = [];
model_kernel(no_clusters).k_h_model_reverb2 = [];
model_kernel(no_clusters).anech_com = [];
model_kernel(no_clusters).reverb1_com = [];
model_kernel(no_clusters).reverb2_com = [];
model_kernel(no_clusters).cluster_id = [];
model_kernel(no_clusters).NPSP = [];
model_kernel(no_clusters).dexp_anech = [];
model_kernel(no_clusters).dexp_reverb1 = [];
model_kernel(no_clusters).dexp_reverb2 = [];
model_kernel(1).params.n_h = n_h;
model_kernel(1).params.time_bin_ms = time_bin_ms;
model_kernel(1).params.NPSP_threshold = NPSP_val;

switch model
    case 'lin'
        model_kernel(1).model_name = 'Linear';
    case 'ln'
        model_kernel(1).model_name = 'LN';
    case 'poisson'
        model_kernel(1).model_name = 'LN-Poisson';
end

n_h = coch(1).params.n_h;
time_bin_ms = coch(1).params.time_bin_ms;
time_ms = [-(n_h*time_bin_ms - time_bin_ms):time_bin_ms:0];

parfor clust = 1:no_clusters
    fprintf('== Generating real and model kernel %0.f/%0.f ==\n',clust,no_clusters);tic;
    psth_length = numel(psth_sturf(clust).anech);
    current_psth = [psth_sturf(clust).anech psth_sturf(clust).reverb1 psth_sturf(clust).reverb2];
    model_kernel(clust).full_kernel_real = sepkerneltensor2(coch_total,current_psth);
    
    
    %Do the convolution to ge the predicted response
    y_hat_t_all = sepconv(coch_total, model_kernel(clust).full_kernel_real);
    
    switch model
        case 'lin'
            y_hat_model_all = y_hat_t_all;
            y_hat_model_anech = y_hat_model_all(1:psth_length);
            y_hat_model_reverb1 = y_hat_model_all(psth_length+1:2*psth_length);
            y_hat_model_reverb2 = y_hat_model_all(2*psth_length+1:3*psth_length);
        case 'ln'
            %Get the ln model using all the real data and predicted model neuron
            model_full_ln = getlnmodel2(y_hat_t_all,current_psth);
            
            %Get the ln model neuron response for all the data together i.e. all
            %conditions will have a single non-linearity
            y_hat_model_all = lnmodelresp(model_full_ln.params, y_hat_t_all);
            
            %Split it into the three reverb conditions
            y_hat_model_anech = y_hat_model_all(1:psth_length);
            y_hat_model_reverb1 = y_hat_model_all(psth_length+1:2*psth_length);
            y_hat_model_reverb2 = y_hat_model_all(2*psth_length+1:3*psth_length);
            
        case 'poisson'
            %Get the ln model using all the real data and predicted model neuron
            model_full_ln = getlnmodel2(y_hat_t_all,current_psth);
            
            %Get the ln model neuron response for all the data together i.e. all
            %conditions will have a single non-linearity
            y_hat_model_ln = lnmodelresp(model_full_ln.params, y_hat_t_all);
            
            %Make a Poisson process neuron
            rand_mat = rand(10,numel(y_hat_model_ln));
            rep_mat_y_hat = repmat(y_hat_model_ln,10,1);
            poisson_mat = rep_mat_y_hat > rand_mat;
            y_hat_model_all = mean(poisson_mat);
            
            %Split it into the three reverb conditions
            y_hat_model_anech = y_hat_model_all(1:psth_length);
            y_hat_model_reverb1 = y_hat_model_all(psth_length+1:2*psth_length);
            y_hat_model_reverb2 = y_hat_model_all(2*psth_length+1:3*psth_length);
    end
    
    %Fit a single sep kernel using all of this new data
    model_kernel(clust).full_kernel_model = sepkerneltensor2(coch_total,y_hat_model_all);
    
    %Fit separate k_hs for the three different conditions
    k_h_model_anech_str = get_k_h(coch(1).anech,y_hat_model_anech, model_kernel(clust).full_kernel_model.k_f);
    k_h_model_anech = k_h_model_anech_str.k_h;
    model_kernel(clust).k_h_model_anech = k_h_model_anech_str;
    
    k_h_model_reverb1_str = get_k_h(coch(1).reverb1,y_hat_model_reverb1, model_kernel(clust).full_kernel_model.k_f);
    k_h_model_reverb1 = k_h_model_reverb1_str.k_h;
    model_kernel(clust).k_h_model_reverb1 = k_h_model_reverb1_str;
    
    k_h_model_reverb2_str = get_k_h(coch(1).reverb2,y_hat_model_reverb2, model_kernel(clust).full_kernel_model.k_f);
    k_h_model_reverb2 = k_h_model_reverb2_str.k_h;
    model_kernel(clust).k_h_model_reverb2 = k_h_model_reverb2_str;
    
    %Compute the double exponential fits
    k_h_model_anech_str.k_h = fliplr(k_h_model_anech_str.k_h');
    [max_val_anech, max_ix_anech] = max(k_h_model_anech_str.k_h);
    [min_val_anech, min_ix_anech] = min(k_h_model_anech_str.k_h);
    Ainit_coeff = max_val_anech + max_val_anech*extra_per_A;
    Binit_coeff = abs(min_val_anech) + abs(min_val_anech)*extra_per_B;
    dexp_anech = fitdexpalex(k_h_model_anech_str.k_h(max_ix_anech:end),n_h,time_bin_ms,Ainit_coeff,Binit_coeff,alphainit_coeff,betainit_coeff,method);
    dexp_anech.beta = dexp_anech.beta + (max_ix_anech-1)*time_bin_ms; % Adjust the time constant to take into account the shifted start
    model_kernel(clust).dexp_anech = dexp_anech;
    
    k_h_model_reverb1_str.k_h = fliplr(k_h_model_reverb1_str.k_h');
    [max_val_reverb1, max_ix_reverb1] = max(k_h_model_reverb1_str.k_h);
    [min_val_reverb1, min_ix_reverb1] = min(k_h_model_reverb1_str.k_h);
    Ainit_coeff = max_val_reverb1 + max_val_reverb1*extra_per_A;
    Binit_coeff = abs(min_val_reverb1) + abs(min_val_reverb1)*extra_per_B;
    dexp_reverb1 = fitdexpalex(k_h_model_reverb1_str.k_h(max_ix_reverb1:end),n_h,time_bin_ms,Ainit_coeff,Binit_coeff,alphainit_coeff,betainit_coeff,method);
    dexp_reverb1.beta = dexp_reverb1.beta + (max_ix_reverb1-1)*time_bin_ms; % Adjust the time constant to take into account the shifted start
    model_kernel(clust).dexp_reverb1 = dexp_reverb1;
    
     k_h_model_reverb2_str.k_h = fliplr(k_h_model_reverb2_str.k_h');
    [max_val_reverb2, max_ix_reverb2] = max(k_h_model_reverb2_str.k_h);
    [min_val_reverb2, min_ix_reverb2] = min(k_h_model_reverb2_str.k_h);
    Ainit_coeff = max_val_reverb2 + max_val_reverb2*extra_per_A;
    Binit_coeff = abs(min_val_reverb2) + abs(min_val_reverb2)*extra_per_B;
    dexp_reverb2 = fitdexpalex(k_h_model_reverb2_str.k_h(max_ix_reverb2:end),n_h,time_bin_ms,Ainit_coeff,Binit_coeff,alphainit_coeff,betainit_coeff,method);
    dexp_reverb2.beta = dexp_reverb2.beta + (max_ix_reverb2 -1)*time_bin_ms; % Adjust the time constant to take into account the shifted start
    model_kernel(clust).dexp_reverb2 = dexp_reverb2;
    
    %Compute the center of mass
    k_h_model_anech = min(k_h_model_anech,0);
    model_kernel(clust).anech_com = time_ms*(k_h_model_anech./sum(k_h_model_anech));
    
    k_h_model_reverb1 = min(k_h_model_reverb1,0);
    model_kernel(clust).reverb1_com = time_ms*(k_h_model_reverb1./sum(k_h_model_reverb1));
    
    k_h_model_reverb2 = min(k_h_model_reverb2,0);
    model_kernel(clust).reverb2_com = time_ms*(k_h_model_reverb2./sum(k_h_model_reverb2));
    
    model_kernel(clust).cluster_id = good_cluster_id(clust);
    model_kernel(clust).NPSP = good_clusters(clust,2);
end
fprintf('== Done! Processing took %.0fs ==\n',toc);
file_name2  = [dir_name,pen,'/data/model_adjbeta'];
save(file_name2,'model_kernel');
%% Plot scatter plots for the double exponential for the three different conditions
if plot_dexp
    
    for cluster = 1:numel(model_kernel)
        anech_dexp_plot(cluster) = model_kernel(cluster).dexp_anech.beta;
        reverb1_dexp_plot(cluster) = model_kernel(cluster).dexp_reverb1.beta;
        reverb2_dexp_plot(cluster) = model_kernel(cluster).dexp_reverb2.beta;
        NPSP(cluster) = model_kernel(cluster).NPSP;
    end
    
    c = NPSP;
    sz = 70;
    
    row = 1;
    col = 3;
    per = 0.03;
    edgel = 0.04; edger = per; edgeh = per; edgeb = 0.05; space_h = per; space_v = per;
    [pos]=subplot_pos(row,col,edgel,edger,edgeh,edgeb,space_h,space_v);
    if store
        figure('units','normalized','outerposition',[0 0 1 1]);
    else
        figure;
    end
    subplot('position',pos{1});
    scatter(reverb1_dexp_plot,reverb2_dexp_plot,sz,c,'filled');
    set(gca,'FontName','Arial','FontSize',10,'FontWeight','Bold');
    xlabel('Double exp Reverb Small [ms]','FontSize',16,'FontWeight','bold');
    ylabel('Double exp Reverb Big [ms]','FontSize',16,'FontWeight','bold');
    axis equal;
    %     xl = xlim;
    %     yl = ylim;
    %     l = [min(xl(1),yl(1)) max(xl(2),yl(2))];
    l = [0 lim_val];
    xlim(l);
    ylim(l);
    
    hline = refline(1,0);
    hline.Color = 'k';
    subplot('position',pos{2});
    scatter(anech_dexp_plot,reverb2_dexp_plot,sz,c,'filled');
    set(gca,'FontName','Arial','FontSize',10,'FontWeight','Bold');
    xlabel('Double exp Anechoic [ms]','FontSize',16,'FontWeight','bold');
    ylabel('Double exp Reverb Big [ms]','FontSize',16,'FontWeight','bold');
    axis equal;
    %     xl = xlim;
    %     yl = ylim;
    %     l = [min(xl(1),yl(1)) max(xl(2),yl(2))];
    xlim(l);
    ylim(l);
    hline = refline(1,0);
    hline.Color = 'k';
    
    subplot('position',pos{3});
    scatter(anech_dexp_plot,reverb1_dexp_plot,sz,c,'filled');
    set(gca,'FontName','Arial','FontSize',10,'FontWeight','Bold');
    xlabel('Double exp Anechoic [ms]','FontSize',16,'FontWeight','bold');
    ylabel('Double exp Reverb Small [ms]','FontSize',16,'FontWeight','bold');
    axis equal;
    %     xl = xlim;
    %     yl = ylim;
    %     l = [min(xl(1),yl(1)) max(xl(2),yl(2))];
    xlim(l);
    ylim(l);
    hline = refline(1,0);
    hline.Color = 'k';
    colormap('jet');
    colorbar;
    
    if store
        file_name  = [dir_name,pen,'/Dexp/',model_kernel(1).model_name,'_Model_neurons_Dexp_',method,'_init_',pen];
        export_fig(file_name);
        close;
    end
    %% Plot scatter plots for the three different conditions for the double exponential fits
    clust_spacing = 9;
    num_groups = ceil(no_clusters/(clust_spacing));
    last_spacing = no_clusters - (num_groups-1)*(clust_spacing) - 1;
    first_clust_ix = [1:clust_spacing:(num_groups)*(clust_spacing)];
    for group = 1:num_groups
        if store
            figure('units','normalized','outerposition',[0 0 1 1]);
        else
            figure;
        end
        first_clust = first_clust_ix(group);
        last_clust = first_clust + clust_spacing-1;
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
        for ii = 1:num_rows*num_columns
            subplot('position',pos{ii});
            ix = plot_ix(ii);
            if ii <= num_clust
                time_ms = [0:time_bin_ms:(numel(model_kernel(ix).k_h_model_anech.k_h) -1)*time_bin_ms];
                plot(time_ms,fliplr(model_kernel(ix).k_h_model_anech.k_h'),'k.','MarkerSize',10);
                hold on;
                time_ms = time_ms(1+numel(model_kernel(ix).k_h_model_anech.k_h) - numel(model_kernel(ix).dexp_anech.F):end);
                plot(time_ms,model_kernel(ix).dexp_anech.F);
                hold off;
                axis tight;
                set(findall(gca, 'Type', 'Line'),'LineWidth',2);
                
                if ii == 1
                    set(gca,'FontName','Arial','FontSize',10,'FontWeight','Bold');
                    ylabel('Weight value','FontSize',16,'FontWeight','bold');
                end
                
            elseif ii > num_clust && ii <= 2*num_clust
                time_ms = [0:time_bin_ms:(numel(model_kernel(ix).k_h_model_reverb1.k_h) - 1)*time_bin_ms];
                plot(time_ms,fliplr(model_kernel(ix).k_h_model_reverb1.k_h'),'k.','MarkerSize',10);
                hold on;
                time_ms = time_ms(1+numel(model_kernel(ix).k_h_model_reverb1.k_h) - numel(model_kernel(ix).dexp_reverb1.F):end);
                plot(time_ms,model_kernel(ix).dexp_reverb1.F,'k');
                hold off;
                axis tight;
                set(findall(gca, 'Type', 'Line'),'LineWidth',2);
                
                if ii == num_clust + 1
                    set(gca,'FontName','Arial','FontSize',10,'FontWeight','Bold');
                    ylabel('Weight value','FontSize',16,'FontWeight','bold');
                end
                
            elseif ii > 2*num_clust && ii <= 3*num_clust
                
                time_ms = [0:time_bin_ms:(numel(model_kernel(ix).k_h_model_reverb1.k_h) - 1)*time_bin_ms];
                plot(time_ms,fliplr(model_kernel(ix).k_h_model_reverb1.k_h'),'k.','MarkerSize',10);
                hold on;
                time_ms = time_ms(1+numel(model_kernel(ix).k_h_model_reverb1.k_h) - numel(model_kernel(ix).dexp_reverb2.F):end);
                plot(time_ms,model_kernel(ix).dexp_reverb2.F,'r');
                hold off;
                axis tight;
                set(findall(gca, 'Type', 'Line'),'LineWidth',2);
                set(gca,'FontName','Arial','FontSize',10,'FontWeight','Bold');
                xlabel('History [ms]','FontSize',14,'FontWeight','bold');
                if ii == 2*num_clust + 1
                    ylabel('Weight value','FontSize',16,'FontWeight','bold');
                end
                
            end
        end
        if store
            file_name  = [dir_name,pen,'/Dexp/',model_kernel(1).model_name,'_Model_neurons_Dexp_fits_',method,'_init_Clust_',num2str(first_clust),'_',num2str(last_clust),'_',pen];
            export_fig(file_name);
            close;
        end
    end
end
%% Plot scatter plot for the COM
if plot_com
    for cluster = 1:numel(model_kernel)
        anech_com_plot_lnmodel(cluster) = model_kernel(cluster).anech_com;
        reverb1_com_plot_lnmodel(cluster) = model_kernel(cluster).reverb1_com;
        reverb2_com_plot_lnmodel(cluster) = model_kernel(cluster).reverb2_com;
        NPSP(cluster) = model_kernel(cluster).NPSP;
    end
    
    c = NPSP;
    sz = 70;
    
    row = 1;
    col = 3;
    per = 0.03;
    edgel = 0.04; edger = per; edgeh = per; edgeb = 0.05; space_h = per; space_v = per;
    [pos]=subplot_pos(row,col,edgel,edger,edgeh,edgeb,space_h,space_v);
    figure('units','normalized','outerposition',[0 0 1 1]);
    subplot('position',pos{1});
    scatter(reverb1_com_plot_lnmodel,reverb2_com_plot_lnmodel,sz,c,'filled');
    set(gca,'FontName','Arial','FontSize',10,'FontWeight','Bold');
    xlabel('COM Reverb Small Model [ms]','FontSize',16,'FontWeight','bold');
    ylabel('COM Reverb Big Model [ms]','FontSize',16,'FontWeight','bold');
    axis equal;
    xl = xlim;
    yl = ylim;
    l = [min(xl(1),yl(1)) max(xl(2),yl(2))];
    xlim(l);
    ylim(l);
    
    hline = refline(1,0);
    hline.Color = 'k';
    subplot('position',pos{2});
    scatter(anech_com_plot_lnmodel,reverb2_com_plot_lnmodel,sz,c,'filled');
    set(gca,'FontName','Arial','FontSize',10,'FontWeight','Bold');
    xlabel('COM Anechoic Model [ms]','FontSize',16,'FontWeight','bold');
    ylabel('COM Reverb Big Model [ms]','FontSize',16,'FontWeight','bold');
    axis equal;
    xl = xlim;
    yl = ylim;
    l = [min(xl(1),yl(1)) max(xl(2),yl(2))];
    xlim(l);
    ylim(l);
    hline = refline(1,0);
    hline.Color = 'k';
    
    subplot('position',pos{3});
    scatter(anech_com_plot_lnmodel,reverb1_com_plot_lnmodel,sz,c,'filled');
    set(gca,'FontName','Arial','FontSize',10,'FontWeight','Bold');
    xlabel('COM Anechoic Model [ms]','FontSize',16,'FontWeight','bold');
    ylabel('COM Reverb Small Model [ms]','FontSize',16,'FontWeight','bold');
    axis equal;
    xl = xlim;
    yl = ylim;
    l = [min(xl(1),yl(1)) max(xl(2),yl(2))];
    xlim(l);
    ylim(l);
    hline = refline(1);
    hline.Color = 'k';
    colormap('jet');
    colorbar;
    if store
        file_name  = [dir_name,pen,'/COM/',model_kernel(1).model_name,'_Model_neurons_COM_',pen];
        export_fig(file_name);
        close;
    end
end
%% Plotting k_hs and full kernels
if plot_res
    clust_spacing = 9;
    time_spacing = 3;
    freq_spacing = 2;
    num_groups = ceil(no_clusters/(clust_spacing));
    last_spacing = no_clusters - (num_groups-1)*(clust_spacing) - 1;
    first_clust_ix = [1:clust_spacing:(num_groups)*(clust_spacing)];
    num_rows = 3;
    num_columns = 3;
    
    %% Plot the k_h component
    for group = 1:num_groups
        
        first_clust = first_clust_ix(group);
        last_clust = first_clust + clust_spacing-1;
        if group == num_groups
            last_clust = first_clust + last_spacing;
        end
        plot_vec = [first_clust:last_clust];
        num_clust = numel(first_clust:last_clust);
        per = 0.02;
        edgel = 0.04; edger = per; edgeh = per; edgeb = 0.05; space_h = 0.01; space_v =0.01;
        [pos]=subplot_pos(num_rows,num_columns,edgel,edger,edgeh,edgeb,space_h,space_v);
        if store
            figure('units','normalized','outerposition',[0 0 1 1]);
        else
            figure;
        end
        for cluster = 1:num_clust
            subplot('position',pos{cluster});
            hold on;
            clust_val = plot_vec(cluster);
            plot(model_kernel(clust_val).k_h_model_anech.k_h);
            plot(model_kernel(clust_val).k_h_model_reverb1.k_h,'k');
            plot(model_kernel(clust_val).k_h_model_reverb2.k_h,'r');
            hold off;
            set(findall(gca, 'Type', 'Line'),'LineWidth',2);
            hline = refline(0,0);
            hline.Color = [0.5 0.5 0.5];
            
            xticks([2:time_spacing:n_h]);
            time_ms = [-(n_h*time_bin_ms - time_bin_ms):time_spacing*time_bin_ms:0];
            for jj = 1:numel(time_ms)
                x_labels{jj} = num2str(time_ms(jj),'%.0f');
            end
            xticklabels(x_labels);
            
            if cluster == num_columns*(num_rows - 1) + 1
                axis on;
                set(gca,'FontName','Arial','FontSize',10,'FontWeight','Bold');
                xlabel('History [ms]','FontSize',14,'FontWeight','bold');
                ylabel('Weight value','FontSize',16,'FontWeight','bold');
                legend('Anech','Small','Big','Location','northwest');
            end
            xlim([1 n_h]);
        end
        if store
            file_name  = [dir_name,pen,'/k_h/',model_kernel(1).model_name,'_Model_Neurons_Clusters_',num2str(first_clust),'_',num2str(last_clust),'_',pen];
            export_fig(file_name);
            close;
        end
        %% Plot the full kernel
        if store
            figure('units','normalized','outerposition',[0 0 1 1]);
        else
            figure;
        end
        for cluster = 1:num_clust
            subplot('position',pos{cluster});
            clust_val = plot_vec(cluster);
            output = model_kernel(clust_val).full_kernel_model.k_f*model_kernel(clust_val).full_kernel_model.k_h';
            mao = max(abs(output(:)));
            imagesc(output,[-mao mao]);
            colormap('redblue');
            
            freqs = coch(1).params.f;
            num_freqs = numel(freqs);
            xticks([1:time_spacing:n_h]);
            freqs = ceil(freqs)/1000; %Convert the frequencies into kHz
            freqs = freqs([1:freq_spacing:num_freqs]);
            time_ms = [-(n_h*time_bin_ms - time_bin_ms):time_spacing*time_bin_ms:0];
            
            for jj = 1:numel(time_ms)
                x_labels{jj} = num2str(time_ms(jj),'%.0f');
            end
            
            for ii = 1:numel(freqs)
                y_labels{ii} = num2str(freqs(ii),'%.1f');
            end
            
            xticklabels(x_labels);
            yticks([1:freq_spacing:num_freqs]);
            yticklabels(y_labels);
            
            if cluster == num_columns*(num_rows - 1) + 1
                axis on;
                set(gca,'FontName','Arial','FontSize',10,'FontWeight','Bold');
                xlabel('History [ms]','FontSize',14,'FontWeight','bold');
                ylabel('Frequency [kHz]','FontSize',16,'FontWeight','bold');
            end
        end
        if store
            file_name  = [dir_name,pen,'/sep/',model_kernel(1).model_name,'_Model_Neurons_Clusters_',num2str(first_clust),'_',num2str(last_clust),'_',pen];
            export_fig(file_name);
            close;
        end
    end
end
end