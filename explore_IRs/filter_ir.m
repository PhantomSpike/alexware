function data_bp = filter_ir(data,Fs,varargin)
%Function filters IRs 
%>>INPUT>>
%3 optional input arguments
%lowfreq
%highfreq
%filter_order


% only want 3 optional inputs at most
numvarargs = length(varargin);

if numvarargs > 3
    error('filter_ir:TooManyInputs', ...
        'requires at most 3 optional inputs');
end

% set defaults for optional inputs
optargs = {200 20000 8};

% now put these defaults into the valuesToUse cell array, 
% and overwrite the ones specified in varargin.
optargs(1:numvarargs) = varargin;

% Place optional args in memorable variable names
[lowfreq, highfreq, filt_order] = optargs{:};


%Bandpass both ears
for channel = 1:2
    data_bp(:,channel) = bandpass(data(:,channel),Fs,lowfreq,highfreq,filt_order);
end

