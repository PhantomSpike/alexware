function [img_xyz,info] = tiff_loader(fname)

info = imfinfo(fname);
num_frames = numel(info);
for z = 1:num_frames
    fprintf('== Loading frame %0.f/%0.f ==\n',z,num_frames);
    img_xyz(:,:,z) = imread(fname, z, 'Info', info);
    % ... Do something with the image ...
end