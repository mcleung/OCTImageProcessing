% outputs rProfile, a mxn matrix where m = # of radius steps and n = # of
% slices analyzed

function rProfile = extractRadialProfile(cell3D, cellMask, minPNslice, ...
    maxPNslice, rList, xCenter, yCenter)



% extract intensity profile as a function of distance from center
% Use 10 bins
rIntervals = 0:.1:1;
rProfile = zeros(length(rIntervals)-1, (maxPNslice - minPNslice - 1));

% calculate current distance from center
[iC jC] = meshgrid(1:size(cell3D,2), 1:size(cell3D,1));
currDistCenter = sqrt((iC - xCenter).^2 + (jC - yCenter).^2);


for i = minPNslice+1:maxPNslice-1
    
    currImBlur = imfilter(cell3D(:,:,i), fspecial('gaussian', 20, 10));
    currImMedFilt = medfilt2(currImBlur, [15 15]);
    currImMedFilt = medfilt2(currImMedFilt, [15 15]);
    currImMedFilt = imadjust(currImMedFilt, [.2, .6], [0 1]);
    
    % also calculate radial profile outward
    currR = rList(i);
    rStepList = rIntervals .* currR;
    
    for j = 1:length(rIntervals)-1
        
        % curr radius range is from rStepList(j):rStepList(j+1)
        currSliceMask = (currDistCenter > rStepList(j)) & (currDistCenter < rStepList(j+1));
        intensityMask = currImMedFilt .* currSliceMask .* cellMask(:,:,i);
        rProfile(j,i-minPNslice) = mean(intensityMask(intensityMask > 0));
        
    end
    
    
end

