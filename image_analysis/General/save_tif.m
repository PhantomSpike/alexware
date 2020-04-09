function save_tif(img,filename,save_dir)

% Save the corrected image as a .tif
outputFileName = [save_dir,filename,'.tif'];
%Convert to 16 bit image
final_img = uint16(double(intmax('uint16'))*mat2gray(img));
%Write to tif file
imwrite(final_img, outputFileName, 'WriteMode', 'append',  'Compression','none');