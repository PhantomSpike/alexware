function ROIMask = get_roi_mask(folder_name)

file_name = fullfile(folder_name,'Fall.mat');
load(file_name);
file_name = fullfile(folder_name,'iscell.npy');
iscell = readNPY(file_name);
ROIMask = zeros([512 512]);

for i = 1:length(stat)
    if iscell(i,1)
        ipix = stat{i}.ypix + (stat{i}.xpix-1) .* 512;
        ROIMask(ipix) = i;
    end
end