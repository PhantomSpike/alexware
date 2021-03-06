count = 0;
%% Model 1
save_dir = '/mnt/40086D4C086D41D0/Reverb_neuronal_data/Kernel_fits/kfolds';
model = 'ridge';
normz = 'perfreq';
kfold = 1;
dt_ms = 10;
h_max_ms = 200;
n_cores = 20;

% try
make_simlnp_model(model,h_max_ms,dt_ms,normz,save_dir,kfold,n_cores);
% catch
%     count = count + 1;
%     neuro_fail{count} = [num2str(model),'_',num2str(h_max_ms),'ms_',fname];
% end

