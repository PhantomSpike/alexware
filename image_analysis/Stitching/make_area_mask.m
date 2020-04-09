function [area_mask,area_names] = make_area_mask(fpath,params,save_dir)


%Params
max_x = params.max_x;
min_x = params.min_x;
max_y = params.max_y;
min_y = params.min_y; 
um_per_pixel = params.um_per_pixel;
size_pixels = params.size_pixels;
zoom_largest_FOV = params.zoom_largest_FOV;
%Get the file location and names
areaList = dir(fpath);
areaList = areaList(~ismember({areaList.name},{'.','..','Tilling','Zstack','misc','processed','proc_log.txt'}));
num_areas = numel(areaList);
area_names = {};


%Main loop
for area = 1:num_areas
    area_names{area} = areaList(area).name(end-1:end); %Get the exact name of the area from the golfer name
    
    if ~strcmp(area_names{area},'10')
        area_names{area} = strrep(area_names{area},'0','');
    end
    
    insideList = dir(fullfile(areaList(area).folder,areaList(area).name));
    insideList = insideList(~ismember({insideList.name},{'.','..'}));
    grabapth = fullfile(insideList(1).folder,insideList(1).name);
    load(grabapth);
    pos_xy_um(area,:) = GRABinfo.xyzPosition(1,1:2);
    zoom_factor(area) = GRABinfo.scanZoomFactor;
end

%Convert from um to pixel space
pos_xy_pixel = round(pos_xy_um./um_per_pixel);

%Initialize big matrix of NaNs to store all the tiles for averaging
store_mat_yxt = zeros(max_y - min_y,max_x - min_x,num_areas);

for area = 1:num_areas
    size_pixels_ind = (size_pixels/zoom_factor(area))*zoom_largest_FOV; %Divide the size of area if zoom was smaller than the image matrix
    half_size_pixels = round(size_pixels_ind/2);
    y_old = pos_xy_pixel(area,2);
    x_old = pos_xy_pixel(area,1);
    %Find the position in the new coordinate system
    y_new = abs(min_y) + y_old;
    x_new = max_x - x_old;
    %Insert the tile in the correct place in the big matrix
    store_mat_yxt(y_new-(half_size_pixels-1):y_new+half_size_pixels,x_new-(half_size_pixels-1):x_new+half_size_pixels,area) = area.*ones(2*half_size_pixels,2*half_size_pixels);
end

%Rotate by 90 degrees counter-clockwise
area_mask = rot90(store_mat_yxt);
save([save_dir '/area_mask.mat'],'area_mask');