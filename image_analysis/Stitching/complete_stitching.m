tilepath = '/home/alex/Desktop/Data/Tile_Anna';
fpath =  '/home/alex/Desktop/Data/Anna_data/';
register = 1; %Whether to register the tiles or not
save_dir = [tilepath '/Stitched'];
divfactor = 0.3;

[stitched_image,params] = image_stitcher_center(tilepath,register);
close all;
area_mask = make_area_mask(fpath,params,save_dir);

displayimgnmask2_alex(stitched_image,area_mask,save_dir,divfactor);