function coch_sound(sound_file,ear,place)
[data,fs] = audioread(sound_file);
% place = 'Original';
switch ear
    case 'left'
        data = data(:,1);
    case 'right'
        data = data(:,2);
end
example_coch(data,fs);
export_fig(['/home/alex/Desktop/Cricket/',place,'/',ear,'_coch.jpg'])