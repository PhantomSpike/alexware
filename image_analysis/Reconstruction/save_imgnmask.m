folder_name = '/home/alex/Desktop/Data/Zoidberg_data/area01/suite2p/plane0/';
save_folder = '/home/alex/Desktop/Data/Zoidberg_data/img_mask/';
save_name = 'Zoidberg_area01';

ROIMask = get_roi_mask(folder_name);
reg_img = get_reg_img(folder_name);

filename = [save_folder,save_name,'.tif'];
final_image = uint8(255*mat2gray(reg_img));
imwrite(final_image, filename);
filename2 = [save_folder,'ROIMask.mat'];
save(filename2,'ROIMask');



