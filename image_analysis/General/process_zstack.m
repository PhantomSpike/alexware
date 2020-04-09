%% Open file
fpath = '/media/alex/5FC39EAD5A6AA312/Micron_imaging/Last_Troubleshoot/Test_stain_last/Leela_U/(20190418_12_04_35)-_Leela_U_zstack/';
save_dir = '/media/alex/5FC39EAD5A6AA312/Micron_imaging/Last_Troubleshoot/Test_stain_last/Leela_U/zstack_processed/';
[fname,pname,tag]=uigetfile({'*.tif'},'Select File',fpath);

grabfname=[fname(1:end-8),'_GRABinfo.mat'];

load([pname grabfname])
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

im=zeros(imHeight,imWidth,total_frames,'uint16');
%% Open raw images
fprintf('Extracting images...\n');tic;

count=0;
for i = 1:nIm
    fprintf('== Opening image %0.f/%0.f ==\n',i,nIm);
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
        [im(:,:,count),~,~]=rtifc(args(i));
    end
end
fprintf('== Done! Extraction took %0.fs ==\n',toc);

%% Store in cell aray according to z-plane
ix_z = [1:frames_per_z:total_frames];
img_mat = cell(stackNumSlices,1);

for slice = 1:stackNumSlices
    img_mat{slice} = im(:,:,ix_z(slice):ix_z(slice)+(frames_per_z-1));
end
%% Register
fprintf('Registering images...\n');tic;
parfor slice = 1:stackNumSlices
    fprintf('== Processing slice %0.f/%0.f ==\n',slice,stackNumSlices);
    img_mat{slice} = register_img(img_mat{slice});
    img_mat{slice} = mean(img_mat{slice},3);
end
fprintf('== Done! Registration took %0.fs ==\n',toc);

%% Save the results as a .tif
fprintf('Writing Tiff file...\n');tic;
outputFileName = [save_dir,'reg_zstack.tif'];
for slice=1:stackNumSlices
   %Convert to 16 bit image
   current_slice = uint16(double(intmax('uint16'))*mat2gray(img_mat{slice}));
   %Rotate by 90 degrees counter-clockwise
   current_slice = rot90(current_slice);
   %Write to tif file
   imwrite(current_slice, outputFileName, 'WriteMode', 'append',  'Compression','none');
end
fprintf('== Done! Writing took %0.fs ==\n',toc);