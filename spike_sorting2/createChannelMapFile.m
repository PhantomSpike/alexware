function createChannelMapFile(Nchan,option)


switch option
    case 'option1'
        Nchannels = 385; % 385th is the synchronisation signal
        xcoords   = [repmat([32 11 39 18], [1 (Nchannels-1)/4]) 0]'; %These are the x-coords for option 1
    case 'option2'
        Nchannels = 385; % 385th is the synchronisation signal
        xcoords   = [repmat([43 11 59 27], [1 (Nchannels-1)/4]) 0]';
    case 'option3'
        Nchannels = 385; % 385th is the synchronisation signal
        xcoords   = [repmat([43 11 59 27], [1 (Nchannels-1)/4]) 0]';
    case 'option4'
        Nchannels = 277; % 277th is the synchronisation signal
        xcoords   = [repmat([43 11 59 27], [1 (Nchannels-1)/4]) 0]';
    otherwise
        error('Invalid option entry');
end

connected = true(Nchannels, 1);
connected(Nchan+1:Nchannels-1) = false;
connected(end) = false;

if strcmp(option,'option4')
    connected([37 76 113 152 189 228 265]) = false;
else
    connected([37 76 113 152 189 228 265 304 341 380]) = false;
end

chanMap   = (1:Nchannels)';
chanMap0ind = chanMap - 1;

ycoords   = [floor((0:Nchannels-2)/2)*20 0]';
% kcoords   = ones(Nchannels,1); % grouping of channels (i.e. tetrode groups)

shankInd = ones(Nchannels,1); 
name = ['Neuropixels Phase 3a ',option];
fs = 30000; % sampling frequency
save('chanMap.mat', ...
    'chanMap','chanMap0ind','connected','name','shankInd','xcoords', 'ycoords','fs');


%%

% kcoords is used to forcefully restrict templates to channels in the same
% channel group. An option can be set in the master_file to allow a fraction 
% of all templates to span more channel groups, so that they can capture shared 
% noise across all channels. This option is

% ops.criterionNoiseChannels = 0.2; 

% if this number is less than 1, it will be treated as a fraction of the total number of clusters

% if this number is larger than 1, it will be treated as the "effective
% number" of channel groups at which to set the threshold. So if a template
% occupies more than this many channel groups, it will not be restricted to
% a single channel group. 