tilepath = '/home/alex/Desktop/Data/Tiles';
register = 0;
%Get the file location and names
stimList = dir(tilepath);
stimList = stimList(~ismember({stimList.name},{'.','..','config_dir'}));
num_tiles = numel(stimList);

meanimg_xy = cell(num_tiles,1);
unit8_max = double(intmax('uint8'));
%Main loop
parfor tile = 1:num_tiles
    fprintf('== Tile no %.0f/%.0f ==\n',tile,num_tiles);
    %% Load the raw files
    fpath = fullfile(stimList(tile).folder,stimList(tile).name,[stimList(tile).name '.tif']);
    [rawimg_xyz,info] = tiff_loader(fpath);
    
    %% Register
%     if register
%         sz=size(rawimg_xyz);
%         frtoavg = 20;
%         Reg_stack = zeros(sz(1),sz(2),1,'uint16');
%         reg = zeros(sz(1),sz(2),frtoavg,'uint16');
%         xyshift = zeros(frtoavg, 2);
%         avg=rawimg_xyz(:,:,1);
%         for fr=1:sz(3);
%             [output, Greg ] = dftregistration(fft2(avg(50:end-50,50:end-50)),fft2(rawimg_xyz(50:end-50,50:end-50,fr)),1);%IG cropping edges
%             xyshift(fr, :)  = [output(3) output(4)];
%             reg_error(fr)=output(1);
%             transvec=[output(3),output(4)];%MFI
%             se=translate(strel(1),transvec);%MFI
%             Reg_stack(:,:,1) = uint16(imdilate(rawimg_xyz(:,:,fr),se));%MFI
%             reg(:,:,fr) = Reg_stack(:,:,1);
%             if fr<= frtoavg
%                 avg=mean(reg,3);
%             end
%         end
%         figure, imagesc(avg), colormap('gray');
%     end

    %% Take the mean
    mean_rawimg_xy = mean(rawimg_xyz,3);
    corr_factor = max(mean_rawimg_xy(:))/unit8_max;
    meanimg_xy{tile} = uint8(mean_rawimg_xy./corr_factor);
end

%%  Register Image Pairs
points = detectSURFFeatures(meanimg_xy{1});
[features, points] = extractFeatures(meanimg_xy{1}, points);

% Initialize all the transforms to the identity matrix. Note that the
% projective transform is used here because the building images are fairly
% close to the camera. Had the scene been captured from a further distance,
% an affine transform would suffice.
numImages = num_tiles;
tforms(numImages) = projective2d(eye(3));

% Initialize variable to hold image sizes.
imageSize = zeros(numImages,2);

% Iterate over remaining image pairs
for n = 2:numImages

    % Store points and features for I(n-1).
    pointsPrevious = points;
    featuresPrevious = features;

    % Save image size.
    imageSize(n,:) = size(meanimg_xy{n});

    % Detect and extract SURF features for I(n).
    points = detectSURFFeatures(meanimg_xy{n});
    [features, points] = extractFeatures(meanimg_xy{n}, points);

    % Find correspondences between I(n) and I(n-1).
    indexPairs = matchFeatures(features, featuresPrevious, 'Unique', true);

    matchedPoints = points(indexPairs(:,1), :);
    matchedPointsPrev = pointsPrevious(indexPairs(:,2), :);

    % Estimate the transformation between I(n) and I(n-1).
    tforms(n) = estimateGeometricTransform(matchedPoints, matchedPointsPrev,...
        'projective', 'Confidence', 99.9, 'MaxNumTrials', 2000);

    % Compute T(n) * T(n-1) * ... * T(1)
    tforms(n).T = tforms(n).T * tforms(n-1).T;
end

% Compute the output limits  for each transform
for i = 1:numel(tforms)
    [xlim(i,:), ylim(i,:)] = outputLimits(tforms(i), [1 imageSize(i,2)], [1 imageSize(i,1)]);
end

% Next, compute the average X limits for each transforms and find the image that is in the center. 
% Only the X limits are used here because the scene is known to be horizontal. 
% If another set of images are used, both the X and Y limits may need to be used to find the center image.

avgXLim = mean(xlim, 2);
avgYLim = mean(ylim, 2);

[~, idx] = sort(avgXLim);
[~, idy] = sort(avgYLim);

centerIdx = floor((numel(tforms)+1)/2);

centerImageIdx = idx(centerIdx);

