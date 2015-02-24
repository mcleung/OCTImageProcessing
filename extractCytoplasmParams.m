% This function extracts cytoplasmic "clumping" parameters for a given
% embryo
%
% Usage: paramList = extractCytoplasmParams(cell3D, cellMask, minPNslice, maxPNslice)
%
% paramList contains:


function paramList = extractCytoplasmParams(cell3D, cellMask, ...
    minPNslice, maxPNslice, rList, xCenter, yCenter)

% make version of cellMask with NaNs instead of just logical
cellMaskN = double(cellMask);
cellMaskN(cellMask == 0) = NaN;

% for i = 1:size(cellMask,3)
%     currSlice = cellMaskN(:,:,i);
%     currSlice(currSlice == 0) = NaN;
%     cellMaskN(:,:,i) = currSlice;
% end

% extract stdev and entropy of all pixels in cytoplasm for each slice
paramList = struct();
stdList = zeros(1,(maxPNslice - minPNslice - 1));
entropyList = zeros(1,(maxPNslice - minPNslice - 1));
entropyValsList = zeros(1,(maxPNslice - minPNslice - 1));

% get indices 4,7,15,18 ... or all indices
offsets = [[zeros(10,1) , 3*[1:10]']; [3*[1:10]' , zeros(10,1)]];
% glcmHomogeneityList = zeros(4, (maxPNslice - minPNslice - 1));

% extract intensity profile as a function of distance from center
% Use 20 bins
rIntervals = 0:.05:1;
rProfile = zeros(length(rIntervals)-1, (maxPNslice - minPNslice - 1));

% calculate current distance from center
[iC jC] = meshgrid(1:size(cell3D,2), 1:size(cell3D,1));
currDistCenter = sqrt((iC - xCenter).^2 + (jC - yCenter).^2);

% fig = figure;
% slicePlotList = round(minPNslice+1:1:maxPNslice-1);

for i = minPNslice+1:maxPNslice-1
    
    % find entropy and stdev of pixels inside mask of blurred and median
    % filtered image
    currImBlur = imfilter(cell3D(:,:,i), fspecial('gaussian', 20, 10));
    currImMedFilt = medfilt2(currImBlur, [15 15]);
    currImMedFilt = medfilt2(currImMedFilt, [15 15]);
    currImMedFilt = imadjust(currImMedFilt, [.2, .6], [0 1]);
    
    se = strel('disk', 15);
    currEntropy = entropyfilt(currImMedFilt, se.getnhood);
    
    %     imageThresh = graythresh(currImMedFilt(~isnan(cellMaskN(:,:,i))));
    %     currEntropyThresh = entropyfilt((currImMedFilt > imageThresh), se.getnhood);
    
    G = graycomatrix(cell3D(:,:,i).*cellMaskN(:,:,i), ...
        'Offset', offsets, 'Symmetric', true);
    stats = graycoprops(G, 'Contrast Correlation Energy Homogeneity');
    
    glcmHomogeneityList(:,i-minPNslice) = stats.Homogeneity';%([4 7 15 18])';
    glcmEnergyList(:,i-minPNslice) = stats.Energy';
    
    % also calculate radial profile outward
    currR = rList(i);
    rStepList = rIntervals * currR;
    
    for j = 1:length(rIntervals)-1
        
        % curr radius range is from rStepList(j):rStepList(j+1)
        currSliceMask = (currDistCenter > rStepList(j)) & ...
            (currDistCenter < rStepList(j+1));
        intensityMask = currImMedFilt .* currSliceMask .* cellMask(:,:,i);
        rProfile(j,i-minPNslice) = mean(intensityMask(...
            ~isnan(intensityMask) & intensityMask > 0));
        
    end
    
    %     figure(fig);
    %     subplot(2,2,1);
    % %     imagesc(cell3D(:,:,i));
    % %     colorbar;
    %     plot(stats.Contrast); title('Contrast');
    %     subplot(2,2,2);
    % %     imagesc(currImMedFilt);
    % %     colorbar;
    %     plot(stats.Correlation); title('Correlation');
    %     subplot(2,2,3);
    % %     imagesc(currEntropy .* cellMask(:,:,i));
    % %     colorbar;
    %     plot(stats.Energy); title('Energy');
    %     subplot(2,2,4);
    %     plot(stats.Homogeneity); title('Homogeneity');
    
    %     imageThresh
    %     pause(.1);
    
    % plot some of the slices along with profiles
%     if (~isempty(intersect(i,slicePlotList)))
%         figure(fig);
%         subplot(2,1,1);
%         imshow(currImMedFilt .* cellMask(:,:,i));
%         subplot(2,1,2);
%         hold on;
%         plot(rProfile(:,i-minPNslice))
%         ylim([.2 .8]);
%     end
%     
    entropyValsList(i-minPNslice) = mean(mean(currEntropy(~isnan...
        (cellMaskN(:,:,i)))));
    entropyList(i-minPNslice) = entropy(currImMedFilt.*cellMaskN(:,:,i));
    
    stdPixelList = currImMedFilt.*cellMaskN(:,:,i);
    stdPixelList = stdPixelList(~isnan(stdPixelList));
    stdList(i-minPNslice) = std(stdPixelList);
    
end

paramList.stdev = stdList;
paramList.entropy = entropyList;
paramList.entropyfilt = entropyValsList;
paramList.homogeneity = glcmHomogeneityList;
paramList.energy = glcmEnergyList;
paramList.rProfile = rProfile;













