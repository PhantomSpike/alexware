function CC = get_cc(stn1_s,stn2_s,t_bin_ms,t_w_ms)
%This was the way that Kerry proposed to run the correlations. I don't
%think this is correct

num_spikes_n1 = length(stn1_s);
t_ms = [-t_w_ms:t_bin_ms:t_w_ms];
stn1_ms = stn1_s.*1000;
stn2_ms = stn2_s.*1000;
cc_psth = zeros(num_spikes_n1,numel(t_ms));

for spike = 1:num_spikes_n1
    cc_psth(spike,:) = histc(stn2_ms,stn1_ms(spike) + t_ms);
end

CC = mean(cc_psth);
zero_bin = ceil(length(t_ms)/2);
CC(zero_bin) = 0;
