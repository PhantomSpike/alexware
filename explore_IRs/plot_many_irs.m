function plot_many_irs(ir_path,save_dir)

ir_list = dir(ir_path);
ir_list = ir_list(~ismember({ir_list.name},{'.','..'}));

num_ir = length(ir_list);

for jj = 1:num_ir
    file_path = fullfile(ir_list(jj).folder,ir_list(jj).name);
    plot_IR(file_path,save_dir);
end