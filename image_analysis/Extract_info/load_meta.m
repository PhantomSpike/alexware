function [meta] = load_meta(raw_dir,save_dir)
%This function extracts useful meta info from many recordings so this can
%be easily read and fed to suite2p for processing later
%>>INPUT>>
%raw_dir - Folder containing the data from a given animal. 
%The assumed folder tree structure is Animal--->Areaxx--->Stimx
%save_dir - Where you want to save the excel and Matlab meta files
%<<OUTPUT<<
%meta - struct with fields corresponding to the different info of interest
%Also saves a .mat and .xls file with the same info

if ~exist('save_dir','var') %If user doesn't specify save_dir make it the one with the data
    save_dir = raw_dir;
end

areaList = dir(raw_dir);
areaList = areaList(~ismember({areaList.name},{'.','..'}));

for area = 1:numel(areaList)
    areapath = fullfile(areaList(area).folder,areaList(area).name);
    stimList = dir(areapath);
    stimList = stimList(~ismember({stimList.name},{'.','..'}));
    grabtemp = dir(fullfile(stimList(1).folder,stimList(1).name,'*GRABinfo.mat')); %Take the GRABInfo from the first folder in this area. THere might be super small diff between them but it should be roughly constant accross stimuli
    grabpath = fullfile(grabtemp.folder,grabtemp.name);
    load(grabpath);
    meta(area).name = areaList(area).name; %Get the name of the area
    
    %Get the number of planes. Sometimes if z-stack was collected before a
    %given area the number of slices will be wrongly written as the number
    %of z-stack planes. Here we correct for this
    if GRABinfo.stackNumSlices > 5
        meta(area).num_planes = 1;
    else
        meta(area).num_planes = GRABinfo.stackNumSlices;
    end
    
    meta(area).fs = round(GRABinfo.scanFrameRate,2); %Get the frame rate
    meta(area).fs_per_plane = round(GRABinfo.scanFrameRate/meta(area).num_planes,2); %Get the frame rate with 2 sig figs
    meta(area).zoom = GRABinfo.scanZoomFactor; %Get the zoom factor
    if meta(area).zoom == 2
        meta(area).diameter = 12; %Get the recommended diameter for this zoom factor
        meta(area).min_neuropil = 350; %Get the recommended diameter for this zoom factor
    elseif meta(area).zoom == 3
        meta(area).diameter = 17; %Get the recommended diameter for this zoom factor
        meta(area).min_neuropil = 686; %Get the recommended diameter for this zoom factor
    elseif meta(area).zoom == 4
        meta(area).diameter = 20; %Get the recommended diameter for this zoom factor
        meta(area).min_neuropil = 896; %Get the recommended diameter for this zoom factor
    end
end
meta(1).tau_6f = 0.7; %Save the tau values for all indicators
meta(1).tau_6m = 1;
meta(1).tau_6s = 1.25;
%Save as excel file
[~,animal_name] = fileparts(raw_dir);
filename = [animal_name,'_metainfo','.xlsx'];
save_name = fullfile(save_dir,filename);
writetable(struct2table(meta),save_name);
%Save as .mat file
filename = [animal_name,'_metainfo'];
save_name = fullfile(save_dir,filename);
save(save_name,'meta');

        
            
    
    
    