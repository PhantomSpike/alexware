function reg_img = get_reg_img(folder_name)

file_name = fullfile(folder_name,'Fall.mat');
load(file_name);
reg_img = ops.meanImg;