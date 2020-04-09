function fix_scanlines(file_path,shift_pix)

[img_xyz,info] = tiff_loader(file_path);
[save_dir,filename] = fileparts(file_path);

im_new = img_xyz;
im_new(:,1:2:end) = circshift(im_new(:,1:2:end),shift_pix);
figure;
imagesc(img_xyz);
title('Original image');
figure;
imagesc(im_new);
title(['Corrected image with ',num2str(shift_pix),' pixels shift']);

% Save the corrected image as a .tif
outputFileName = fullfile(save_dir,[filename,'_cor_',num2str(shift_pix),'pixshift.tif']);
%Convert to 16 bit image
final_img = uint16(double(intmax('uint16'))*mat2gray(im_new));
%Write to tif file
imwrite(final_img, outputFileName, 'WriteMode', 'append',  'Compression','none');


