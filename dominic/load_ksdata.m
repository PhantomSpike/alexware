sorted_dir = '/media/basile/09F6A5BC59F7E39C/Dominic_Data/Ronnie/P06/P06-quning';
chan_map = '/home/basile/Desktop/Code/Preprocessing_pipeline_Neuropixel/Kilosort2/configFiles/neuropixPhase3A_kilosortChanMap.mat';
setenv('NEUROPIXEL_MAP_FILE',chan_map);
ks = Neuropixel.KiloSortDataset(sorted_dir);
ks.load();
stats = ks.computeBasicStats();
ks.printBasicStats();
metrics = ks.computeMetrics();