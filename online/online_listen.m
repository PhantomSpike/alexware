function online_listen(start_time,end_time,ch_no,data_root,probe)
%This function loads a selected channel from a given bin data set and plays
%it from the speakers 
%% Load the data
fs = 30000;
[cropped_data] = load_spikedata(data_root,start_time,end_time,probe);
cropped_data = double(cropped_data);
cropped_data = cropped_data - median(cropped_data,2);
%% Plot the channel
t = [start_time:1/fs:end_time];
figure;
plot(t,cropped_data);
xlabel('Time [s]');
ylabel('Amplitude');
title(['Channel #',]);
%% Play the sound
ch = num2str(ch_no);
fprintf('== Playing sound from channel %s ==\n', ch);
sound(cropped_data(ch_no,:),fs);
end