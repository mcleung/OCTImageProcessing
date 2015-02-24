function dataFinal = removeLineArtifacts(dataGauss, coords)

xmin = coords(1);
xmax = coords(2);
ymin = coords(3);
ymax = coords(4);

% find noise levels in rectangle with noise only
% this is the reference frame to scale all other noise to
backgroundLevel = squeeze(mean(mean(dataGauss(ymin:ymax, xmin:xmax, :),1),2));

% first replace all points with sum over z above a certain threshold (45)

dataSum = sum(dataGauss,3);
dataSum = varNorm(dataSum);

% figure, surf(dataSum)
% shading interp;

% record (x,y) with values above threshold (look for outliers in image
% hist) (find 2 largest peaks)

[counts, x] = imhist(dataSum);
[p h] = mspeaks(x, counts);

peakLocs = p(:,1);
peakHeights = p(:,2);

[pHeightSort, inds] = sort(peakHeights, 'descend');
pLocSort = peakLocs(inds);
h = h(inds,:);

pLocSort = pLocSort(1:2);
h = h(1:2,:);

if pLocSort(1) > pLocSort(2)
    cutoff = pLocSort(1) + 4*(h(1,2) - pLocSort(1));
    cutoffTaper = pLocSort(1);% + 2*(h(1,2) - pLocSort(1));
else
    cutoff = pLocSort(2) + 4*(h(2,2) - pLocSort(2));
    cutoffTaper = pLocSort(2);% + 2*(h(2,2) - pLocSort(2));
end

mask = dataSum > cutoff;
% figure, imshow(mask)
mask = imdilate(mask, strel('disk', 5));
iList = 1:size(dataGauss,1);
jList = 1:size(dataGauss,2);
[I,J] = meshgrid(iList,jList);
I = I';
J = J';

iList = I(mask);
jList = J(mask);
se = strel('disk', 15);

% for each z location, take those x,y, convolve with gaussian, and subtract
% from image, then take abs value

% that should remove all remaining artifacts

dataGaussR = reshape(dataGauss,size(dataGauss,1)*size(dataGauss,2), ...
    size(dataGauss,3));
dataFinalR = dataGaussR;
maskR = reshape(mask, size(dataGauss,1)*size(dataGauss,2), 1);

% OR, blur mask, scale each region by local height of original artifact,
% and then subtract from original image

% % 1. Find max value in original image of each region
% dataFinal = dataGauss;
% maskDilate = imdilate(mask, strel('disk', 10));
% L = bwlabel(maskDilate);
% LR = reshape(L, size(dataGauss,1)*size(dataGauss,2), 1);
% S = regionprops(maskDilate);
% maskBlur = imfilter(double(mask), fspecial('gaussian', 20, 5));
% maskBlurR = reshape(maskBlur, size(dataGauss,1)*size(dataGauss,2), 1);
% regionMaxOrig = zeros(length(S), size(dataGauss,3));
% regionMaxBlur = zeros(length(S), 1);
% 
% % find max of each region after blurring, and max of each region in slices
% % then multiply blurred image by original mask with each region scaled
% % appropriately
% for ii = 1:length(S)
%     regionMaxOrig(ii,:) = max(dataGaussR(LR==ii,:),[],1) - backgroundLevel';
%     regionMaxBlur(ii,:) = max(maskBlurR(LR==ii));
% end
% 
% for k = 50%1:size(dataGauss,3)
%    
%     k
%     maskBlurScaledR = maskBlurR;
%     
%     % scale the brightness of each blurred artifact by its brightness in
%     % that slice of dataGauss
%     for ii = 1:length(S)
%         maskBlurScaledR(LR==ii) = maskBlurScaledR(LR==ii)*...
%             regionMaxOrig(ii,k)/regionMaxBlur(ii);
%     end
%     
%     maskBlurScaled = reshape(maskBlurScaledR, size(dataGauss,1), size(dataGauss,2));
%     dataFinal(:,:,k) = abs(dataGauss(:,:,k) - maskBlurScaled);
%     
%     figure, imshow(dataFinal(:,:,k));
%     
% end
%

% do all slices at once; iterate over bright points in mask
for ii = 1:length(iList)
        
    i = iList(ii);
    j = jList(ii);
    
    % if point is artifact
    % replace it by median of points surrounding it in a circle
    % that are not part of artifact
    mask2 = zeros(size(dataSum));
    mask2(i,j) = 1;
    
    % to find where mask2(i,j) maps to when reshaped
    mask2R = reshape(mask2, size(dataGauss,1)*size(dataGauss,2), 1);
    mask2R = (mask2R == 1);
    
    % now dilate and reshape
    mask2D = imdilate(mask2, se);
    mask2DR = reshape(mask2D, size(dataGauss,1)*size(dataGauss,2), 1);
    mask2DR = (mask2DR == 1);
    
    % make sure no NaNs in output
    if sum(mask2DR & ~maskR) == 0 % & ~maskR
        % now dilate and reshape
        mask2D = imdilate(mask2, strel('disk', 20));
        mask2DR = reshape(mask2D, size(dataGauss,1)*size(dataGauss,2), 1);
        mask2DR = (mask2DR == 1);
    end
    
    % take median along first dim and replace point value in all slices at
    % once
    % actually, don't replace, but do weighted sum based on how high above
    % threshold original value is
    origVals = dataGaussR(mask2R,:);
    newVals = mean(dataGaussR(mask2DR & ~maskR,:),1); % & ~maskR
    
%     if i == 463 && j == 619
%         1
%     end
    
    a = sigmf(origVals, [20, cutoffTaper]);
%     a = 1 - exp(-20*(origVals - cutoffTaper));
%     a(a > 1) = 1; 
%     a(a < 0) = 0;
    dataFinalR(mask2R,:) = a.*newVals + (1-a).*origVals;
    
    if mod(ii,100)==0
        fprintf('\n %f percent complete' , ceil(ii/length(iList)*100));
    end
    
end

dataFinal = reshape(dataFinalR, size(dataGauss,1), size(dataGauss,2), ...
    size(dataGauss,3));

for k = 1:size(dataGauss,3)
    
    dataFinal(:,:,k) = medfilt2(dataFinal(:,:,k),[5 5]);
    
end

fprintf('\n\n');
















