function fft_sound(sound_file,ear,place)

% place = 'Original';

[data,fs] = audioread(sound_file);
switch ear
    case 'left'
        data = data(:,1);
    case 'right'
        data = data(:,2);
end
[P1,f] = fft_run(data,fs);
export_fig(['/home/alex/Desktop/Cricket/',place,'/',ear,'_ft.jpg']);