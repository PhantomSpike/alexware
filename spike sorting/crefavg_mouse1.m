function processed_data = crefavg_mouse1(data,correction)
% Process the data to do Common Reference Averaging (CRA) i.e. subtract the
% median from each channel across all time from itself and then subtract the median at each time point across all channels from every channel.
%>>Input>>
%data - The data that you want to process. The format assumed is
%NchannelsxNtimepoints
%correction - A flag that determines whether to remove the common median
%for every time point across channels or not
%<<Output<<
%processed_data - The processed data


fprintf('== CRA processing ==\n');tic;

if nargin < 2
    correction = false;
end

%Extract the synch channel and remove it from the data
synch_ch = data(end,:);
data(end,:) = [];
%Find the median across all time for each channel and subtract it from that channel to remove offset
data = data - median(data,2);

if correction
    %Find the median across all channels for each time point to remove common
    %noise across channels and subtract it from the data
    data = data - median(data,1);
end


data(end+1,:) = synch_ch;
processed_data = flipud(data);
fprintf('== Done! Processing took %f sec ==\n',toc);