function [img_final_xyz] = reconstruct_coronal2(img_path)

%% Initialize things
theta = -65; %The angle for rotating the z-axis
ref_slice = 1; %Which slice to use for the ref image
extra_pix = 2; %How many extra pixels to add to the big matrix of zeros/NaN
save_folder = [img_path '/Processed'];
meta_folder = [img_path '/Meta'];

if exist(save_folder,'dir')
    rmdir(save_folder, 's');
end

% if exist(meta_folder,'dir')
%     rmdir(meta_folder, 's');
% end


%Get the file location and names
dir_info = dir(fullfile(img_path,'**/*.tif'));
num_slices = length(dir_info);

img_xyz = cell(1);

for slice = 1:num_slices
    fname = fullfile(dir_info(slice).folder,dir_info(slice).name);
    info1 = imfinfo(fname);
    y_dim(slice) = info1(1).Height; 
    x_dim(slice) = info1(1).Width; 
end

y_max = max(y_dim);
y_half = round(y_max/2);
x_max = max(x_dim);
x_half = round(x_max/2);
img_xyz = cell(num_slices,1);
img_mean_xy = cell(num_slices,1);

fprintf('== Running the Reconstruction ==\n');tic;

%% Load the raw files

for slice = 1:num_slices
    fprintf('== Processing slice %.0f/%.0f ==\n',slice,num_slices);
    fpath = fullfile(dir_info(slice).folder,dir_info(slice).name);
    [temp_xyz,~] = tiff_loader(fpath);
    sz = size(temp_xyz);
    temp_y = sz(1);
    temp_x = sz(2);
    temp_y_half = floor(temp_y/2);
    temp_x_half = floor(temp_x/2);
    img_xyz{slice} = zeros(y_max+extra_pix,x_max+extra_pix,sz(3));
    y_start = y_half-(temp_y_half-1);
    x_start = x_half-(temp_x_half-1);
    img_xyz{slice}(y_start:y_start+(temp_y-1),x_start:x_start+(temp_x-1),:) = temp_xyz;
%     img_mean_xy{slice} = mean(img_xyz{slice},3);
end

fprintf('== Done! Loading took %0.fs ==\n',toc);

%% Register the images
show_img = 0;
ref_img = mean(img_xyz{ref_slice},3);
slice_list = [1:num_slices];
slice_list(ref_slice) = [];
tform = cell(num_slices,1);
img_xyz_t = cell(num_slices,1);

parfor slice = slice_list
    fprintf('== Registering slice %.0f/%.0f ==\n',slice,num_slices);tic;
    [~,tform{slice}] = reg_img_alex(ref_img,mean(img_xyz{slice},3),show_img);
    fprintf('== Done! Registration took %.0fs ==\n',toc);
    num_zstacks = size(img_xyz{slice},3);
    fprintf('== Transforming slice %.0f/%.0f ==\n',slice,num_slices);tic;
    for z = 1:num_zstacks
        img_xyz_t{slice}(:,:,z) = imwarp(img_xyz{slice}(:,:,z),tform{slice},'OutputView',imref2d(size(ref_img)));
    end
    fprintf('== Done! Transformation took %.0fs ==\n',toc);
end



filename_met = [meta_folder '/tform.mat'];
filename_met2 = [meta_folder '/img_xyz_trans.mat'];
save(filename_met,'tform');
save(filename_met2,'img_xyz');

%% Rotate the image
fprintf('== Rotating the image ==\n');tic;
count = 0;
downsampling_factor = 1;
for slice = 1:num_slices
    num_zstacks = size(img_xyz{slice},3);
    for z = 1:num_zstacks
        count = count + 1;
        temp_img = imrotate(img_xyz{slice}(1:downsampling_factor:end,1:downsampling_factor:end,z),theta);
        temp_img = uint16(double(intmax('uint16'))*mat2gray(temp_img));
        img_final_xyz(:,:,count) = temp_img'; %Transpose the final image so it is loaded correctly in ImageJ
    end
end
fprintf('== Done! Rotation took %.0fs ==\n',toc);
clear img_xyz;

addpath(genpath('/home/alex/Fiji.app/scripts/'));
ImageJ;
IJM.show('img_final_xyz');
keyboard;
% mkdir(save_folder);
% filename = [save_folder '/Naja_processed.tif'];
% num_z = size(img_final_xyz,3);
% 
% for z = 1:num_z
%      imwrite(img_final_xyz(:,:,z), filename, 'WriteMode', 'append',  'Compression','none');
% end
% 
% fprintf('== Done! Writing took %.0fs ==\n',toc);

ij.IJ.run("Quit",""); %Close ImageJ

% theta_y = theta_y_deg*(pi/180); %Convert the angle to radians
% theta_z = theta_z_deg*(pi/180); %Convert the angle to radians
% R_x = [1 0 0 0; ...
%        0 cos(theta_x) -sin(theta_x) 0; ...
%        0 sin(theta_x) cos(theta_x) 0; ...
%        0 0 0 1];
   
% R_y = [cos(theta_y) 0 sin(theta_y) 0; ... %Matrix for rotation in 3D around the y-axis
%        0 1 0 0; ...
%        -sin(theta_y) 0 cos(theta_y) 0; ...
%        0 0 0 1];
%    
% R_z = [cos(theta_z) -sin(theta_z) 0 0; ... %Matrix for rotation in 3D around the z-axis
%        sin(theta_z) cos(theta_z) 0 0; ...
%        0 0 1 0; ...
%        0 0 0 1];


% mkdir(save_folder);