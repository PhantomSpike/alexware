function stitch_fun_areas(tilepath,areapath)
%Enter the folder name where the tiles are saved and the folder name where
%the GRABinfo files for the different areas are saved

%Don't change these
save_dir = [tilepath '/Stitched'];
divfactor = 0.3;
load([save_dir '/params.mat']);
load([save_dir '/stitched_image_final.mat']);
[area_mask,area_names] = make_area_mask(areapath,params,save_dir);
displayimgnmask2_alex(stitched_image_final,area_mask,save_dir,area_names,divfactor);