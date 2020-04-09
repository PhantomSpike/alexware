function [final_image,params] = image_stitcher2(tilepath)

save_folder = [tilepath '/Stitched'];
fov_um = 1000;

if exist(save_folder,'dir')
    rmdir(save_folder, 's');
end

%Get the file location and names
stimList = dir(tilepath);
stimList = stimList(~ismember({stimList.name},{'.','..'}));
num_tiles = numel(stimList);

meanimg_xy = cell(num_tiles,1);

fprintf('== Running the stitching on tiles ==\n');tic;
%Main loop
for tile = 1:num_tiles
    fprintf('== Processing tile no %.0f/%.0f ==\n',tile,num_tiles);
    %% Load the raw files
    fpath = fullfile(stimList(tile).folder,stimList(tile).name,[stimList(tile).name '.tif']);
    [rawimg_xyz,info] = tiff_loader(fpath);
    grabapth = fullfile(stimList(tile).folder,stimList(tile).name,[stimList(tile).name '_GRABinfo.mat']);
    load(grabapth);
    pos_xy_um(tile,:) = GRABinfo.xyzPosition(1,1:2);
    %% Take the mean
    meanimg_xy{tile} = mean(rawimg_xyz,3);
end
%% Find the size of the image in x and y
x_size_pixel = GRABinfo.scanPixelsPerLine;
y_size_pixel = GRABinfo.scanPixelsPerLine;
%Find the um per pixel and convert xy posiiton from um to pixels
um_per_pixel = fov_um/x_size_pixel;
pos_xy_pixel = round(pos_xy_um./um_per_pixel);

%Find the range of values
max_x = max(pos_xy_pixel(:,1)) + x_size_pixel;
min_x = min(pos_xy_pixel(:,1)) - x_size_pixel;
max_y = max(pos_xy_pixel(:,2)) + y_size_pixel;
min_y = min(pos_xy_pixel(:,2)) - y_size_pixel;

%Initialize big matrix of NaNs to store all the tiles for averaging
store_mat_yxt = nan(max_y - min_y,max_x - min_x,num_tiles);

%Fill the planes of the matrix t with the tiles according to their xy
%position
for tile = 1:num_tiles
    y_old = pos_xy_pixel(tile,2);
    x_old = pos_xy_pixel(tile,1);
    %Find the position in the new coordinate system
    y_new = abs(min_y) + y_old;
    x_new = max_x - x_old;
    %Insert the tile in the correct place in the big matrix
    store_mat_yxt(y_new:y_new+y_size_pixel-1,x_new:x_new+x_size_pixel-1,tile) = meanimg_xy{tile};
end

%Take the average of the big matrix
stitched_image = nanmax(store_mat_yxt,[],3);
%Smooth with Gaussian kernel
% stitched_smooth = imgaussfilt(stitched_image,1);
%Rotate by 90 degrees counter-clockwise
final_image = rot90(stitched_image);
%Replace NaN entries with 0s
final_image(isnan(final_image))=0;
%Save the image as a .tif file
mkdir(save_folder);
filename = [save_folder '/Stitched.tif'];
final_image = uint8(255*mat2gray(final_image));
imwrite(final_image, filename);
figure;
imagesc(final_image);
colormap('gray');
%Store all the parameters
params.max_x = max_x;
params.min_x = min_x;
params.max_y = max_y;
params.min_y = min_y;
save([save_folder '/params'],'params');
save([save_folder '/stitched_image'],'final_image');
fprintf('== Done! Stitching took %.0f s ==\n',toc);
