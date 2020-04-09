function displayimgnmask2_alex(Image,area_mask,save_dir,area_names,divfactor)

if nargin<5
    divfactor=0.3;
end

num_areas = size(area_mask,3);
im_class = class(Image);

figure('units','normalized','outerposition',[0 0 1 1]);
Image = imresize(Image,divfactor);

area_mask = imresize(area_mask,divfactor,'nearest');

for area = 1:num_areas
    perim_mask_log(:,:,area) = bwperim(area_mask(:,:,area));
end


perim_mask_log = max(perim_mask_log,[],3);

switch im_class
    case 'uint8'
        perim_mask = uint8(perim_mask_log).*intmax(im_class);
    case 'uint16'
        perim_mask = uint16(perim_mask_log).*intmax(im_class);
end

Image_areas = max(Image,perim_mask);

%     imshow((I_in.*~bwperim(lbl_mask>0))',[min(I_in(:)) max(I_in(:))./divfactor],...
%         'InitialMagnification','fit');
imshow((Image_areas),...
    'InitialMagnification','fit'); %Alex 31/05/2018


for i=1:num_areas
    current_mask = squeeze(area_mask(:,:,i));
    centd_mask=zeros(num_areas,2);
    [r,c]=find(current_mask==i);
    rc=[c,r];
    centd_mask(i,:)=mean(rc,1);
    text(centd_mask(i,1)-2,centd_mask(i,2),['\fontsize{9}\color{red}' num2str(area_names{i})]);
end

set(gca,'box','off');
export_fig([save_dir '/Stitched_mask.tif']);
end
