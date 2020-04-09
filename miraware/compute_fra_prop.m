function [fra_prop] = compute_fra_prop(fra,base_rate,params)
%[clust_info] = compute_fra_stats(fra_psth)
%This function computes characteristic frequency (CF), best frequency (BF), Q10
%and Q30. We first smooth the FRA with a smoothing window. Then we
%threshold the resulting FRA with the spontaneous rate of each neuron + 20%
%of the maximum activity from the smoothed FRA. We then define FRA 'bounds'
%as the lowest level at each freqeuncy that is significantaly higher above this threshold 
%>>INPUT>>
%fra_psth
%<<OUTPUT<<
freqs = params.freqs;
db_levels = params.dB_levels;
%% First perform smoothing of the FRA
s_win=[0.25 0.5 0.25;0.5 1 0.5; 0.25 0.5 0.25]; %Consider to pad this properly !!!!!!
s_win = s_win/sum(s_win(:)); % Normalise so it sums to 1
fra_smooth = conv2(fra,s_win,'same'); % Smooth the FRA

%% Define base rate and compute bounds
%Define the threshold firing rate. This will be the spontaneous (base) rate
% + 20% of the max smoothed FRA rate
th_rate = base_rate + 1/5*max(fra_smooth(:));
fra_th = fra_smooth - th_rate; %Subtract the threshold from all FRA entries
ix_sig = find(fra_th>0); %Find the index of the freq-level combinations that are significantly higher than the threshold
[lvl_ix,f_ix] = ind2sub(size(fra),ix_sig); %Convert the linear ix to rows (levels) and cols (frequencies) indices
sig_f_ix = unique(f_ix); %Find all the freqeuncies that have at least one significant bin in the fra

%For every significant freqeuncy find the lowest dB level that is
%significant
if ~isempty(sig_f_ix)
    for f = 1:numel(sig_f_ix)
        min_lvl_ixs(f) = max(lvl_ix(f_ix==sig_f_ix(f))); %Note because lvls go high->low but indexing is 1->n we use the max
    end
else
    min_lvl_ixs = [];
end

%% Find the BF, CF, Q10 and Q30 (if possible)
%BF
[~,bf_ix] = max(mean(fra));
bf = freqs(bf_ix);

%CF
[cf_lvl_ix,~] = max(min_lvl_ixs); %Find the maximum value for the index of the dB level at the lowest threshold
if ~isempty(cf_lvl_ix)
    cf_ix = sig_f_ix(min_lvl_ixs==cf_lvl_ix);
    cf = freqs(cf_ix); %Get the freqeuncy/ies that correspond to this lowest level/s
else
    cf = [];
    cf_ix = [];
end


if length(cf)>1
    cf = 2^mean(log2(cf)); %If there is more than one freq at the lowest threshold find the log-weighted mean
    cf_ix = 2^mean(log2(sig_f_ix(min_lvl_ixs==cf_lvl_ix))); %Get the index of cf in the same way for plotting purposes
end

%Q10
q10_lvl_ix = cf_lvl_ix - 1;
if q10_lvl_ix>0
    q10_f_ix = unique(f_ix(lvl_ix==q10_lvl_ix));
    q10_f = freqs(q10_f_ix);
    bw_10 = max(q10_f) - min(q10_f); %First find the bandwidth at 10dB louder than the lowest threshold level
    q10 = cf/bw_10;
else
    q10 = [];
    bw_10 = [];
    q10_f_ix = [];
end

%Q30
q30_lvl_ix = cf_lvl_ix - 3;
if q30_lvl_ix>0
    q30_f_ix = unique(f_ix(lvl_ix==q30_lvl_ix));
    q30_f = freqs(q30_f_ix);
    bw_30 = max(q30_f) - min(q30_f); %Then find the bandwidth at 30dB louder than the lowest threshold level
    q30 = cf/bw_30;
else
    q30 = [];
    bw_30 = [];
    q30_f_ix = [];
end

%% Save the results
lvl_plot_ix = zeros(length(freqs),1);
lvl_plot_ix(sig_f_ix) = min_lvl_ixs;

fra_prop.lvl_plot_ix = lvl_plot_ix;
fra_prop.f_plot_ix = [1:length(freqs)]';
fra_prop.bf = bf; 
fra_prop.cf = cf; 
fra_prop.bf_ix = bf_ix; 
fra_prop.cf_ix = cf_ix; 
fra_prop.q10 = q10; 
fra_prop.q30 = q30; 
fra_prop.bw_10 = bw_10; 
fra_prop.bw_30 = bw_30;
fra_prop.q10_f_ix = q10_f_ix;
fra_prop.q10_lvl_ix = q10_lvl_ix*ones(length(q10_f_ix));
fra_prop.q30_f_ix = q30_f_ix;
fra_prop.q30_lvl_ix = q30_lvl_ix*ones(length(q30_f_ix));
fra_prop.db_th = db_levels(cf_lvl_ix);
fra_prop.fra_smooth = fra_smooth;
