function [thisDat] = getLFPdata(lfpFilename,nChansInFile)
% function [lfpByChannel, allPowerEst, F] = lfpBandPower(lfpFilename, lfpFs, nChansInFile, freqBand)
% Computes the power in particular bands, and across all frequencies, across the recording
% samples 10 segments of 10 sec each to compute these things. 


% load nClips one-sec samples
d = dir(lfpFilename); 
nSamps = d.bytes/2/nChansInFile;


mmf = memmapfile(lfpFilename, 'Format', {'int16', [nChansInFile nSamps], 'x'});



thisDat = double(mmf.Data.x(:, (1:end)));
    
thisDat = bsxfun(@minus, thisDat, mean(thisDat,2)); %median subtraction
    
    


    