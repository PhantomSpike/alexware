main_dir = '/mnt/40086D4C086D41D0/Reverb_data/Done/Process';

animal_list = dir(main_dir);
animal_list = animal_list(~ismember({animal_list.name},{'.','..'}));

for animal = 1:numel(animal_list)
    fprintf('== Processing animal %s ==\n',animal_list(animal).name);
    pen_dir = dir([fullfile(animal_list(animal).folder,animal_list(animal).name),'/**/*ap.bin']);
    for pen = 1:length(pen_dir)
        sorted_dir = pen_dir(pen).folder;
        fprintf('== File: %s ==\n',pen_dir(pen).name);
        plot_cl_info(sorted_dir);
    end
end