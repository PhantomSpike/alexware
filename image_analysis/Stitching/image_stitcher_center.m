function [stitched_image_final,params] = image_stitcher_center(tilepath,register,smooth)

save_folder = [tilepath '/Stitched'];
fov_zoom1_um = 1000;
im_class = 'uint16';

if exist(save_folder,'dir')
    rmdir(save_folder, 's');
end

%Get the file location and names
stimList = dir(tilepath);
stimList = stimList(~ismember({stimList.name},{'.','..'}));
num_tiles = numel(stimList);

img_xy = cell(num_tiles,1);

fprintf('== Running the stitching on tiles ==\n');tic;
%Main loop
for tile = 1:num_tiles
    fprintf('== Processing tile no %.0f/%.0f ==\n',tile,num_tiles);
    %% Load the raw files
    ftemp = dir(fullfile(stimList(tile).folder,stimList(tile).name,'*.tif'));
    fpath = fullfile(ftemp.folder,ftemp.name);
    [img_xy{tile},info] = tiff_loader(fpath);
    grabtemp = dir(fullfile(stimList(tile).folder,stimList(tile).name,'*GRABinfo.mat'));
    grabpath = fullfile(grabtemp.folder,grabtemp.name);
    load(grabpath);
    pos_xy_um(tile,:) = GRABinfo.xyzPosition(1,1:2);
    zoom_factor_tile(tile) = GRABinfo.scanZoomFactor; %Find the zoom factor for every tile
    %% Take the mean
%     meanimg_xy{tile} = mean(rawimg_xyz,3);
end
%% Register
if register
    parfor tile = 1:num_tiles
        img_xy{tile} = register_img(img_xy{tile});
        img_xy{tile} = mean(img_xy{tile},3);
    end
else
    for tile = 1:num_tiles
        img_xy{tile} = mean(img_xy{tile},3);
    end
end
%% Find the size of the image in x and y
size_pixels = GRABinfo.scanPixelsPerLine;
half_size_pixels = size_pixels/2;
%Find the um per pixel and convert xy posiiton from um to pixels
zoom_largest_FOV = min(zoom_factor_tile); %Find the zoom that corresponds to the largest FOV
fov_um = round(fov_zoom1_um/zoom_largest_FOV); %Find the field of view that corresponds to the largest zoom    
um_per_pixel = fov_um/size_pixels;
pos_xy_pixel = round(pos_xy_um./um_per_pixel);

%Find the range of values
max_x = max(pos_xy_pixel(:,1)) + half_size_pixels;
min_x = min(pos_xy_pixel(:,1)) - half_size_pixels;
max_y = max(pos_xy_pixel(:,2)) + half_size_pixels;
min_y = min(pos_xy_pixel(:,2)) - half_size_pixels;

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
    %Insert the tile in the correct place in the big matrix:
    %If the zoom of the tile is higher than the smallest zoom, resample the
    %image to be of the the right size to be inserted in the matrix
    if zoom_factor_tile(tile) > zoom_largest_FOV
        zoom_ratio = zoom_factor_tile(tile)/zoom_largest_FOV;
        size_pixels_scaled = round(size_pixels/zoom_ratio);
        half_size_pixels_scaled = round(size_pixels_scaled/2);
        F = griddedInterpolant(img_xy{tile});
        [sx,sy] = size(img_xy{tile});
        xq = (0:zoom_ratio:sx)';
        yq = (0:zoom_ratio:sy)';
        img_scaled = F({xq,yq});
        if length(y_new-(half_size_pixels_scaled-1):y_new+half_size_pixels_scaled) > length(img_scaled)
            store_mat_yxt(y_new-(half_size_pixels_scaled-1):y_new+(half_size_pixels_scaled-1),x_new-(half_size_pixels_scaled-1):x_new+(half_size_pixels_scaled-1),tile) = img_scaled;
        else
            store_mat_yxt(y_new-(half_size_pixels_scaled-1):y_new+half_size_pixels_scaled,x_new-(half_size_pixels_scaled-1):x_new+half_size_pixels_scaled,tile) = img_scaled;
        end
    else
        store_mat_yxt(y_new-(half_size_pixels-1):y_new+half_size_pixels,x_new-(half_size_pixels-1):x_new+half_size_pixels,tile) = img_xy{tile};
    end
end

%Take the average of the big matrix
stitched_image = nanmax(store_mat_yxt,[],3);

if smooth
    %Smooth with Gaussian kernel
    stitched_image = imgaussfilt(stitched_image,1);
end

%Rotate by 90 degrees counter-clockwise
stitched_image_final = rot90(stitched_image);
%Replace NaN entries with 0s
stitched_image_final(isnan(stitched_image_final)) = 0;
%Save the image as a .tif file
mkdir(save_folder);
filename = [save_folder '/Stitched.tif'];

switch im_class
    case 'uint8'
        stitched_image_final = uint8(double(intmax(im_class))*mat2gray(stitched_image_final));
    case 'uint16'
        stitched_image_final = uint16(double(intmax(im_class))*mat2gray(stitched_image_final));
end

imwrite(stitched_image_final, filename);
figure;
imagesc(stitched_image_final);
colormap('gray');
%Store all the parameters
params.max_x = max_x;
params.min_x = min_x;
params.max_y = max_y;
params.min_y = min_y;
params.um_per_pixel = um_per_pixel;
params.size_pixels = size_pixels;
params.zoom_largest_FOV = zoom_largest_FOV;
save([save_folder '/params'],'params');
save([save_folder '/stitched_image_final'],'stitched_image_final');
fprintf('== Done! Stitching took %.0f s ==\n',toc);