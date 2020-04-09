function kurt = ir_kurtosis(data,Fs,t_win_ms)

t_win_samp = (t_win_ms/1000)*Fs;
sz = length(data);
kurt = zeros(sz-t_win_samp,1);

for tau = 1:(sz-t_win_samp)
    kurt(tau) = kurtosis(data(tau:tau + t_win_samp - 1),0);
end



