function CC = get_cc_cor(stn1_s,stn2_s,t_bin_ms,tau)

t_shift = [-tau:tau];
t_end_s = max([stn1_s(:);stn2_s(:)]) + 1;
t_bin_s = t_bin_ms/1000;
t_s = [0:t_bin_s:t_end_s];

psth_n1 = histcounts(stn1_s,t_s);
psth_n2 = histcounts(stn2_s,t_s);
n_bins = length(psth_n1); %Find the length of the histogram

for shift = 1:length(t_shift)
    tau = t_shift(shift);
    psth_n2_temp = circshift(psth_n2,tau);
    mean_psth_n2_temp = mean(psth_n2_temp);
    mean_psth_n1 = mean(psth_n1);
    CC(shift) = (1/n_bins)*(psth_n2_temp-mean_psth_n2_temp)*(psth_n1-mean_psth_n1)';
end

zero_bin = ceil(length(t_shift)/2);
CC = CC./CC(zero_bin);
CC(zero_bin) = 0;
