function write_spikedata(data,dir_name,processed_name)

new_name = fullfile(dir_name,processed_name);
new_name = [new_name,'.ap.bin'];
fid2 = fopen(new_name,'w');
fprintf('== Writing file %s ==\n',new_name);tic;
fwrite(fid2,data,'int16');
fprintf('== Done! Writing took %f sec ==\n',toc);
fclose(fid2);

end