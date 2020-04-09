% This script does batch pre-processing and spike sorting of raw data from
% the Neuropixeles collected with SpikeGLX. 
cra_on = 1; %Logical statement that determines whether we perform CRA on the data
automerge = 0; %Logical statement that determines whether we use automering of the tempaltes

pen_root = '/mnt/40086D4C086D41D0/Kilosort_optimize/Msc_Students/Data/'; % folder containing all stimuli to be sorted, aranged in 1 folder per penetration, with 1 config folder for this peentration and cthe stimuli from that penetration in separate folders

addpath(genpath('/home/alex/Desktop/Code/KiloSort')) % path to kilosort folder
addpath(genpath('/home/alex/Desktop/Code/KiloSort/npy-matlab')) % path to npy-matlab scripts

penList = dir(pen_root);
penList = penList(~ismember({penList.name},{'.','..'}));

for pen = 1:numel(penList)
    penpath = fullfile(penList(pen).folder,penList(pen).name);
    
    stimList = dir(penpath);
    stimList = stimList(~ismember({stimList.name},{'.','..','config_dir'}));
    
    pathToYourConfigFile = [penpath '/config_dir']; % take from Github folder and put it somewhere else (together with the master_file)
    run(fullfile(pathToYourConfigFile, 'StandardConfig.m'));
    
    switch option
        case 'option1'
            ChanTotal = 385;
        case 'option2'
            ChanTotal = 385;
        case 'option3'
            ChanTotal = 385;
        case 'option4'
            ChanTotal = 277;
    end
    
    for stim = 1:numel(stimList)
        fpath = fullfile(stimList(stim).folder,stimList(stim).name);
        if cra_on
            % common reference averaging
            process_raw_spikedata(fpath,ChanTotal);
        end
        
        craFilePath = dir(fullfile(fpath,'CRA','*.ap.bin'));
        craFilePath = craFilePath(~ismember({craFilePath.name},{'.','..'}));
        
        % change StandardConfig for each different penetrations
        ops.fbinary             = fullfile(craFilePath.folder,craFilePath.name); % will be created for 'openEphys'
        ops.fproc               = fullfile(craFilePath.folder,'/temp_wh.dat'); % residual from RAM of preprocessed data
        ops.root                = fpath; % 'openEphys' only: where raw files are
        
        % run KiloSort for each penetration
        tic; % start timer
        %
        if ops.GPU
            gpuDevice(1); % initialize GPU (will erase any existing GPU arrays)
        end
        
        if strcmp(ops.datatype , 'openEphys')
            ops = convertOpenEphysToRawBInary(ops);  % convert data, only for OpenEphys
        end
        %
        [rez, DATA, uproj] = preprocessData(ops); % preprocess data and extract spikes for initialization
        rez                = fitTemplates(rez, DATA, uproj);  % fit templates iteratively
        rez                = fullMPMU(rez, DATA);% extract final spike times (overlapping extraction)
        
        % AutoMerge. rez2Phy will use for clusters the new 5th column of st3 if you run this)
        if automerge
            rez = merge_posthoc2(rez); %AI might want to change this
        end
        
        % save matlab results file
        save(fullfile(ops.root,'CRA','rez.mat'), 'rez', '-v7.3');
        
        % Saving the spike times in mat files
        spike_times = rez.st3(:,1);
        save(fullfile(ops.root,'CRA','spike_times.mat'),'spike_times');
        
        % save python results file for Phy
        rezToPhy(rez, [ops.root,'/CRA']);
        
        % remove temporary file
        delete(ops.fproc);
    end
    
end