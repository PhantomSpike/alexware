tmp_folder = '/home/alex/Desktop/Data/ImageJ_test/temp/'; %Folder where the macro is stored
fname = '/media/alex/5FC39EAD5A6AA312/GAD_labelling/Nala/'; %The name of the folder with all the images
process = 'lr'; %What conversion to do to the data: lr - flip horizontally; ud - flip vertically

javaaddpath '/home/alex/java/jar/mij.jar'; %Add this so we can run Miji
try
    MIJ.exit;
catch
    fprintf('IMageJ already closed\n');
end

%% Get the names of all the .oir files
root_dir_info = dir(fname);
root_dir_info = root_dir_info(~ismember({root_dir_info.name},{'.','..'}));
num_folders = length(root_dir_info);
fid = fopen([tmp_folder,'/tmp.ijm'],'wt'); %Open a pointer to the .ijm file that will be saved

for f = 1:num_folders
    imgf_info = dir(fullfile(root_dir_info(f).folder,root_dir_info(f).name,'*.oir'));
    num_images = length(imgf_info);
    %% This part makes the ImageJ Script. Maybe I can run some sort of for loop in here for each file
    
    for img = 1:num_images
        fname = fullfile(imgf_info(img).folder,imgf_info(img).name);
        fprintf(fid,['run("Bio-Formats Importer", "open=',fname,' autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");\n']);
        save_name = fullfile(imgf_info(img).folder,[root_dir_info(f).name,'_',num2str(img),'.tif']);
        switch process
            case 'lr'
                fprintf(fid,'run("Flip Horizontally", "stack");\n'); %Flip the slide horizontally (fliplr) so un-do the horizontal flipping doen by the microscope
            case 'ud'
                fprintf(fid,'run("Flip Vertically", "stack");\n'); %Flip the slide vertically (flipud) so un-do the vertical flipping doen by the microscope
        end
        fprintf(fid,['saveAs("Tiff","',save_name,'");\n']); %Write command to open the file
        fprintf(fid,'close();\n'); %Open the Image in Matlab
    end
    
end

fclose(fid); %Close the pointer

Miji; %Launch Miji
MIJ.run('Install...', ['install=',tmp_folder,'/tmp.ijm']); %Install the Macro
fprintf('== Processing .oir files ==\n');tic;
MIJ.run('tmp'); %Run the Macro
MIJ.exit; %Close Miji
fprintf('== Done! Writing took %.1f sec ==\n',toc); 