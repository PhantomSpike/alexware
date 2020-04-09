sound_file = '/home/alex/Desktop/Cricket/crickets.wav';
ear = 'left';
Fstop = 1000;
Fpass = 1100;
Astop = 60;
pb_ripple = 0.5;
place = 'Filtered';
[data,fs] = audioread(sound_file);

switch ear
    case 'left'
        data = data(:,1);
    case 'right'
        data = data(:,2);
end

filt_data = highpass_iir(data,Fpass,Fstop,Astop,pb_ripple,fs);

close;
example_coch(filt_data,fs);
export_fig(['/home/alex/Desktop/Cricket/',place,'/',ear,'_coch_',num2str(Fstop),'.jpg']);
close;
[P1,f] = fft_run(filt_data,fs);
export_fig(['/home/alex/Desktop/Cricket/',place,'/',ear,'_ft_',num2str(Fstop),'.jpg']);
close;