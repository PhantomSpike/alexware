function processed_data = crefavg_hack(data,ChanTotal,Chan_border)
% processed_data = crefavg(data,nChan,groupSpan)
% Process the data to do Common Reference Averaging (CRA) i.e. subtract the median from each channel across all time from itself and then subtract the median at each time point across all channels from every channel
% set default values of the arguments if don't exist
% usually 12 is a good value for groupSpan

fprintf('== CRA processing ==\n');tic;
%Extract the synch channel and remove it from the data
synch_ch = data(ChanTotal,:);
data(ChanTotal,:) = [];
%Extract the channels that are outside of the brain
out_ch = data(1:Chan_border,:);
data(1:Chan_border,:) = [];
%Find the median across all time for each channel and subtract it from that channel to remove offset
%median_data = median(data,2);
data = data - median(data,2);
%Find the median across all channels for each time point to remove common
%noise across channels and subtract it from the data
%median_allchs_time = median(data,1);
data = data - median(data,1);

%Flip the channels up-down
data = flipud(data);
out_ch = flipud(out_ch);

%Now we add the normalized data
data = [data;out_ch];
data(ChanTotal,:) = synch_ch;
processed_data = data;
fprintf('== Done! Processing took %f sec ==\n',toc);