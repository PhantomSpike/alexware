function fix_scanlines_test(file_path,shift_pix,save_dir)

[img_xyz,info] = tiff_loader(file_path);

shifts = [-shift_pix:1:shift_pix];
num_shifts = length(shifts);
im_mat = cell(num_shifts);

for j = 1:num_shifts
    curr_shift = shifts(j);
    im_mat{j} = img_xyz;
    im_mat{j}(:,1:2:end) = circshift(im_mat{j}(:,1:2:end),curr_shift);
%     im_mat{j}(1:2:end,:) = circshift(im_mat{j}(1:2:end,:),curr_shift);
%     %Shift in rows instead
    figure('units','normalized','outerposition',[0 0 1 1]);
    imagesc(im_mat{j});
    file_name = ['Shift of',num2str(curr_shift),' pixels'];
    save_name = fullfile(save_dir,[file_name,'.jpg']);
    export_fig(save_name);
    close all;
end








