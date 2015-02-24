% xpts and ypts are coordinates of best fit circle returned by Hough
% transform

% edges is the output of canny edge detection on the original image

function [xfit yfit] = findBestEllipse(x, y, im2, sliceEdge)

% want to find an ellipse that is close to xpts, ypts and best fits edges
edgeDilate = imdilate(sliceEdge, strel('disk', 10));

A = zeros(size(sliceEdge,1), size(sliceEdge,2));
for i = 1:length(x)
    A(y(i), x(i)) = 1;
end

% do some thresholding to only keep points of high intensity that are close
% to detected circle contour
ptsDilate = imdilate(A, strel('disk', 10));
ptsToFit = ptsDilate & edgeDilate;
ptsIntensity = ptsToFit .* im2;
ptsPos = ptsIntensity(ptsIntensity > 0);
ptsIntensity(ptsIntensity < graythresh(ptsPos)) = 0;
ptsIntensity(ptsIntensity > 0) = 1;
contourToFit = imerode(ptsIntensity, strel('disk', 3));

[xPts yPts] = find(contourToFit > 0);

% discard some points so we have fewer to fit
for i = 1:length(xPts)
   if mod(i,5) ~= 0
       xPts(i) = NaN;
       yPts(i) = NaN;
   end
end

xPts = xPts(~isnan(xPts));
yPts = yPts(~isnan(yPts));

% now find optimal ellipse with max amount of overlap with ptsToFit
[z, a, b, alpha] = fitellipse([xPts' ; yPts']);

% figure, imshow(slices);
hold on;
plotellipse(z,a,b,alpha);%, 'LineWidth', 2);