%% Main script for running KS2
pen_root = '/mnt/40086D4C086D41D0/Kilosort_tute/PLP_copy'; % folder containing all stimuli to be sorted, aranged in 1 folder per penetration, with 1 config folder for this peentration and cthe stimuli from that penetration in separate folders
plot_figs = 1; %Whether to plot figures while running KS. If you plot figs the computer becomes unusable 
addpath(genpath('/home/alex/Desktop/Code/Kilosort2')); % path to kilosort folder
addpath('/home/alex/Desktop/Code/npy-matlab');

penList = dir(pen_root);
penList = penList(~ismember({penList.name},{'.','..'}));

for pen = 1:numel(penList)
    penpath = fullfile(penList(pen).folder,penList(pen).name);
    stimList = dir(penpath);
    stimList = stimList(~ismember({stimList.name},{'.','..','config_dir'}));
    meta_dir = dir(fullfile(fullfile(stimList(1).folder,stimList(1).name), '*ap*.meta')); % Get the ap.meta filename
    meta_info = readSpikeGLXmeta(fullfile(meta_dir.folder, meta_dir.name)); % Load the meta data
    option = meta_info.imProbeOpt;
    pathToYourConfigFile = [penpath '/config_dir']; % take from Github folder and put it somewhere else (together with the master_file)
    run(fullfile(pathToYourConfigFile, 'StandardConfig.m'));
    ops.NchanTOT = meta_info.nSavedChans;
    
    if plot_figs
        ops.fig = 1;
    else
        ops.fig = 0;
    end
    
    for stim = 1:numel(stimList)
        % the binary file is in this folder
        fpath = fullfile(stimList(stim).folder,stimList(stim).name);
        ops.fproc       = fullfile(fpath, 'temp_wh.dat'); % proc file on a fast SSD
        
        ops.trange = [0 Inf]; % time range to sort
        
        %% this block runs all the steps of the algorithm
        fprintf('Looking for data inside %s \n', fpath)
        
        % find the binary file
        fs          = [dir(fullfile(fpath, '*.ap.bin')) dir(fullfile(fpath, '*.dat'))];
        ops.fbinary = fullfile(fpath, fs(1).name);
        
        % preprocess data to create temp_wh.dat
        rez = preprocessDataSub(ops);
        
        % time-reordering as a function of drift
        rez = clusterSingleBatches(rez);
        save(fullfile(fpath, 'rez.mat'), 'rez', '-v7.3');
        
        % main tracking and template matching algorithm
        rez = learnAndSolve8b(rez);
        
        % final merges
        rez = find_merges(rez, 1);
        
        % final splits by SVD
        rez = splitAllClusters(rez, 1);
        
        % final splits by amplitudes
        rez = splitAllClusters(rez, 0);
        
        % decide on cutoff
        rez = set_cutoff(rez);
        
        fprintf('found %d good units \n', sum(rez.good>0))
        
        % write to Phy
        fprintf('Saving results to Phy  \n')
        rezToPhy(rez, fpath);
        
        %% if you want to save the results to a Matlab file...
        
        % discard features in final rez file (too slow to save)
        rez.cProj = [];
        rez.cProjPC = [];
        
        % save final results as rez2
        fprintf('Saving final results in rez2  \n')
        fname = fullfile(fpath, 'rez2.mat');
        save(fname, 'rez', '-v7.3');
    end
end
