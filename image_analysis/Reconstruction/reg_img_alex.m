function [reg_moving,tform] = reg_img_alex(fixed,moving,show_img)
%This function registers one image to another and outputs the transformed
%image and the transformation
fprintf('== Registering image ==\n'); tic;
option_optim = 'monomodal';
% 'monomodal' - Monomodal images have similar brightness and contrast. The images are captured on the same type of scanner or sensor
% 'multimodal' - Multimodal images have different brightness and contrast. The images can come from two different types of devices, 
%such as two camera models or two types of medical imaging modalities (like CT and MRI). The images can also come from a single 
%device, such as a camera using different exposure settings, or an MRI scanner using different imaging sequences.
transform_type = 'rigid';
%'translation' - (x,y) translation in 2-D, or (x,y,z) translation in 3-D.
%'rigid' - Rigid transformation consisting of translation and rotation.
%'similarity' - Nonreflective similarity transformation consisting of translation, rotation, and scale.
%'affine' - Affine transformation consisting of translation, rotation, scale, and shear.

[optimizer, metric] = imregconfig(option_optim); %This function gets the optimization object and metric for the registration

tform = imregtform(moving, fixed, transform_type, optimizer, metric); %This function provides the best geometric transformation of the moving image to the ref fixed image based on the parameters provided

fprintf('== Done! Registration took %0.fs ==\n',toc); 

fprintf('== Warping image ==\n'); tic;

reg_moving = imwarp(moving,tform,'OutputView',imref2d(size(fixed)));

fprintf('== Done! Warping took %0.fs ==\n',toc); 

if show_img
    figure;
    imshowpair(fixed, reg_moving,'Scaling','joint');
end