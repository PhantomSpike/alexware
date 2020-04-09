function [clust_info] = feature_extract(clust_info)
%[clust_info] = feature_extract(clust_info)
%This function extracts features from the raw waveforms for later
%clustering
%>>INPUT>>
%clust_info - The clust_info var containing the raw wfs and a lot of
%different metadata
%<<OUTPUT<<
%clust_info - Same var but with the results from the function appended in
%the results field

delay_ms = 0.01; %The delay from the trough before calculating the slope
slope_period_ms = 0.2; %The time window over whcih to calculate the slope
num_clust = size(clust_info.good_units,1);
delta_t_ms = clust_info.int_time_ms(2) - clust_info.int_time_ms(1); %Compute length in ms of one sample
total_dur_samples = length(clust_info.int_time_ms);
count=0; %Counter to count only the 'normal' wfs for indexing purposes
for clust = 1:num_clust
    wf = [];
    wf = clust_info.good_int_wf(clust,:);
    
    [~,ix_abs_max] = max(abs(wf)); %Find ix of the maximum abs value
    
    if wf(ix_abs_max)<0 && ix_abs_max/total_dur_samples<0.6
        count = count + 1;
        clust_info.results.normal_int_wf(count,:) = wf; %Save the 'normal' waveform for each cluster
        wf = -wf; %If the max is -ve then the spike is 'normal'. Invert it for computing purposes. Also make sure the peak is not very late (not later than 60% from the beginning of the waveform)
        clust_info.results.normal_cluster_id(count,1) = clust_info.good_units(clust,1); %Save the id of the waveform if it has a normal (non-inverted) shape
        clust_info.good_units(clust,3) = 1; %Also save logical in the third column of the good_units field
    else
        clust_info.good_units(clust,3) = 0; %logical in the third column of the good_units field
        continue %If the max is +ve i.e. inverted spike, then don't process it
    end
    
    [peak_val,ix_peak] = max(wf); %Find the peak 
    [trough_val,ix_trough] = min(wf(ix_peak:end)); %Find the trough

    %Find spike width
    spike_width_samp = ix_trough; %Find spike width in samples, because the ix is relative to the maximum so it will be already the spike width
    clust_info.results.spike_width_ms(count,1) = spike_width_samp*delta_t_ms; %Convert to ms
    
    %Find the peak-to-trough ratio
    clust_info.results.peak_trough_ratio(count,1) = peak_val/abs(trough_val); 
    
    %Find the absolute spike magnitude
    clust_info.results.abs_spike_mag(count,1) = peak_val + abs(trough_val); 
    
    %Find the end slope
    ix_trough = ix_trough + ix_peak - 1; %Find the ix of the trough relative to the beginning of the waveform
    delay_samples = round(delay_ms/delta_t_ms);
    slope_period_samples = round(slope_period_ms/delta_t_ms); %Find the number of samples for the end-slope calculation 
    ix_endslope_start = ix_trough+delay_samples; %The index in samples for when to begin calcuating the slope
    ix_endslope_end = ix_endslope_start+slope_period_samples; %The index in samples for when to begin calcuating the slope
    if ix_endslope_end>total_dur_samples %Check if the last point of the slope will be further than the end of the wf. If so, make it the end of the wf
        ix_endslope_end = total_dur_samples;
    end
    P = polyfit(clust_info.int_time_ms(ix_endslope_start:ix_endslope_end),wf(ix_endslope_start:ix_endslope_end),1); %Calculate a first order polynomial (y=mx+c) fit to the data 

    clust_info.results.end_slope(count,1) = -P(1); %Save the slope = end slope. Chnage the sign because of the inversion
    
    %Find the width at 1/2 height for the peak (depolarization part of AP)
    half_peak = peak_val/2;
    [~,left_bound_ix_p] = min(abs(wf(1:ix_peak)-half_peak)); %Find the left bound ix i.e. the time point at which the wf is closest to the half peak value before the peak (to the left)
    [~,right_bound_ix_p] = min(abs(wf(ix_peak:end)-half_peak)); %Find the right bound ix i.e. the time point at which the wf is closest to the half peak value after the peak (to the right)
    right_bound_ix_p = right_bound_ix_p + ix_peak -1; %Find the ix of the right bound relative to the begining of the waveform
    width_hh_peak_samples = right_bound_ix_p - left_bound_ix_p; %Find the width at 1/2 height for the peak in samples
    width_hh_peak_ms = width_hh_peak_samples*delta_t_ms; %Convert to ms
    clust_info.results.width_hh_peak_ms(count,1) = width_hh_peak_ms; %Save the width at 1/2 height for the peak
    
    %Find the width at 1/2 height for the trough (hyper-polarization part of AP)
    half_trough = trough_val/2;
    [~,left_bound_ix_t] = min(abs(wf(ix_peak:ix_trough)-half_trough)); %Find the left bound ix i.e. the time point at which the wf is closest to the half peak value before the peak (to the left)
    left_bound_ix_t = left_bound_ix_t + ix_peak - 1; %Get the index relative to the start of  the waveform
    [~,right_bound_ix_t] = min(abs(wf(ix_trough:end)-half_trough)); %Find the right bound ix i.e. the time point at which the wf is closest to the half peak value after the peak (to the right)
    right_bound_ix_t = right_bound_ix_t + ix_trough -1; %Find the ix of the right bound relative to the begining of the waveform
    width_hh_trough_samples = right_bound_ix_t - left_bound_ix_t; %Find the width at 1/2 height for the peak in samples
    width_hh_trough_ms = width_hh_trough_samples*delta_t_ms; %Convert to ms
    clust_info.results.width_hh_trough_ms(count,1) = width_hh_trough_ms; %Save the width at 1/2 height for the trough
    
    %Save other useful results for testing of the algorithm later
    clust_info.results.ix_peak(count,1) = ix_peak;
    clust_info.results.ix_trough(count,1) = ix_trough;
    clust_info.results.y_int(count,1) = -P(2); %Save the y-intercept. Change the sign because of the inversion
    clust_info.results.ix_endslope_start(count,1) = ix_endslope_start;
    clust_info.results.ix_endslope_end(count,1) = ix_endslope_end;
    clust_info.results.left_bound_ix_p(count,1) = left_bound_ix_p;
    clust_info.results.right_bound_ix_p(count,1) = right_bound_ix_p;
    clust_info.results.left_bound_ix_t(count,1) = left_bound_ix_t;
    clust_info.results.right_bound_ix_t(count,1) = right_bound_ix_t;
end