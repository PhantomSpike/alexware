%% Open file
fpath = '/mnt/40086D4C086D41D0/Suite2p_python_data/Raw/Bender/area06/(20181009_15_13_02)-_area06_tonesA/';
high_pass_f = 0.05; %High-pass cut-off frequency
% save_dir = '/media/alex/5FC39EAD5A6AA312/Micron_imaging/Last_Troubleshoot/Test_stain_last/Leela_U/zstack_processed/';
[fname,pname,tag]=uigetfile({'*.tif'},'Select File',fpath);

grabfname=[fname(1:end-8),'_GRABinfo.mat'];

load([pname grabfname])
fs = GRABinfo.scanFrameRate;
% Stimulus variables

nChannels=GRABinfo.channelsSave;

imWidth=GRABinfo.scanPixelsPerLine;
imHeight=GRABinfo.scanLinesPerFrame;
stackNumSlices = GRABinfo.stackNumSlices;
z_step_um = GRABinfo.stackZStepSize;
frames_per_z = GRABinfo.acqNumFrames;
total_frames = nChannels*stackNumSlices*frames_per_z;
frames_per_file = 900;


imList = dir([pname fname(1:end-8) '*.tif']);
nIm = length(imList);

frames_last_file = total_frames - frames_per_file*(nIm - 1);

% im=zeros(imHeight,imWidth,total_frames,'uint16');
resp = zeros(total_frames,1);
%% Open raw images
fprintf('Extracting images...\n');tic;
curr_frame = 0;
count=0;
for i = 1:nIm
    fprintf('== Processing image %0.f/%0.f ==\n',i,nIm);
    args(i).filename=[pname imList(i).name];
    args(i).info=imfinfo([pname imList(i).name]);
    args(i).pixelregion=[];
    
    if i < nIm
        nFrames = frames_per_file;
    else
        nFrames = frames_last_file;
    end
    
    for ii=1:nFrames
        count=count+1;
        args(i).index=ii;
        args(i).offset=args(i).info(ii).Offset;
        [curr_frame,~,~]=rtifc(args(i));
        resp(count) = mean(curr_frame(:));
    end
end
fprintf('== Done! Extraction took %0.fs ==\n',toc);
resp_filt = highpass(resp,fs,high_pass_f,4);