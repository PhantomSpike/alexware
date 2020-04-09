function [results] = fast_reg_old(wfs)
num_clust = size(wfs,1);
fs = 30000;
for jj = 1:num_clust
    wf = [];
    wf = wfs(jj,:);
    [~,ix_abs_max] = max(abs(wf)); %Find ix of the maximum abs value
    if wf(ix_abs_max)<0
        wf = -wf; %If the max is negative, invert the spike
    end
    [peak_val,ix_peak] = max(wf); %Find the peak
    [trough_val,ix_trough] = min(wf); %Find the trough
    spike_width_samp = abs(ix_peak - ix_trough); %Find spike width in samples
    results.spike_width_ms(jj,1) = (spike_width_samp/fs)*1000; %Convert to ms
    results.peak_trough_ratio(jj,1) = peak_val./abs(trough_val); %Find the peak_trough_ratio
end