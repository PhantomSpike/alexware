img_name = '/media/alex/5FC39EAD5A6AA312/Micron_imaging/Reconstruction_ferrets/Naja/Converted/Meta/img_final_xyz.mat';
downsmaple_f = 10;
ij.IJ.run("Quit",""); %Close ImageJ
load(img_name);
temp_img = double(img_final_xyz);
temp_img = (temp_img./max(temp_img(:)))*255;
F = griddedInterpolant(temp_img);
[sx,sy,sz] = size(temp_img);
xq = (0:downsmaple_f:sx)';
yq = (0:downsmaple_f:sy)';
zq = (1:sz)';
vq = uint8(F({xq,yq,zq}));
addpath(genpath('/home/alex/Fiji.app/scripts/'));
ImageJ;
IJM.show('vq');
