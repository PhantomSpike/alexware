function processed_data = crefavg2(data,ChanTotal,nChan,groupSize)
% processed_data = crefavg(data,nChan,groupSpan)
% Process the data to do Common Reference Averaging (CRA) i.e. subtract the median from each channel across all time from itself and then subtract the median at each time point across all channels from every channel
% set default values of the arguments if don't exist
% usually 12 is a good value for groupSpan
if ~exist('nChan')
    nChan = size(data,1)-1;
end
if ~exist('groupSize')
    grouping = 'off';
else
    grouping = 'on';
end
fprintf('== CRA processing ==\n');tic;
%Extract the synch channel and remove it from the data
synch_ch = data(ChanTotal,:);
data(ChanTotal,:) = [];
%Extract the channels that are outside of the brain
out_ch = data(nChan+1:ChanTotal-1,:);
data(nChan+1:ChanTotal-1,:) = [];
%Find the median across all time for each channel and subtract it from that channel to remove offset
%median_data = median(data,2);
data = data - median(data,2);
%Find the median across all channels for each time point to remove common
%noise across channels and subtract it from the data
%median_allchs_time = median(data,1); 
if strcmp(grouping,'off')
    data = data - median(data,1);
elseif strcmp(grouping,'on')
    numGroup = floor(nChan/groupSize);
    for group = 1:numGroup-1
        group_ind = [1:groupSize]+(group-1)*groupSize;
        data(group_ind,:) = data(group_ind,:) - median(data(group_ind,:),1);
    end
    group_ind = [(group_ind(end)+1):nChan];
    data(group_ind,:) = data(group_ind,:) - median(data(group_ind,:),1);    
elseif strcmp(grouping,'sliding')
    new_data = zeros(nChan,size(data,2));
    groupSpan = floor(groupSize/2);
    for chan = 1:nChan
        if chan <= groupSpan
            group_ind = [-chan+1:groupSpan]+chan;
            new_data(chan,:) = data(chan,:)-median(data(group_ind,:),1);
        elseif chan >= nChan - groupSpan
            group_ind = [-groupSpan:nChan-chan]+chan;
            new_data(chan,:) = data(chan,:)-median(data(group_ind,:),1);
        else
            group_ind = [-groupSpan:groupSpan]+chan;
            new_data(chan,:) = data(chan,:)-median(data(group_ind,:),1);        
        end
    end
end
%Now we add the normalized data
data(nChan+1:ChanTotal-1,:) = out_ch;
data(ChanTotal,:) = synch_ch;
processed_data = data;
fprintf('== Done! Processing took %f sec ==\n',toc);