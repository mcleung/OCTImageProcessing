% calculate the density in the neighborhood of each pixel and filter out
% pixels in low density areas

function currW_NS2 = filterNoise(currWAdj, coords)

xmin = coords(1);
xmax = coords(2);
ymin = coords(3);
ymax = coords(4);

% make filter to average pixels in a 10x10 window around a given pixel
H = fspecial('average', [20 20]);
currWAvg = imfilter(currWAdj, H);
% figure, imshow(currWAvg);

noiseLevel = mean(mean(currWAvg(ymin:ymax, xmin:xmax))) + ...
    3*std(std(currWAdj(ymin:ymax, xmin:xmax)));

% figure, imhist(currWAvg(ymin:ymax, xmin:xmax));

% suppress all pixels that don't have nearby bright pixels
% set their values to the gaussian filtered values divided by some factor
currW_NS = currWAdj;
currW_NS(currWAvg < noiseLevel) = currWAvg(currWAvg < noiseLevel).^2;

% figure, imshow(currW_NS);
currWAvg2 = imfilter(currW_NS, fspecial('average', [40 40]));
% figure, imshow(currWAvg2);

% maybe instead of thresholding it taper it off more based on how far from
% the threshold value a given pixel is (if it's above threshold, raise it
% to a smaller power, if it's below threshold raise it to a larger power)
noiseLevel2 = mean(mean(currWAvg2(ymin:ymax, xmin:xmax))) + ...
    3*std(std(currWAdj(ymin:ymax, xmin:xmax)));

% see how far away the value at currWAvg2 is from the noise level
% pixesl 
currWAvg2_thresh = exp(1*(noiseLevel2 - currWAvg2));
currW_NS2 = currWAdj .^ (currWAvg2_thresh); % currWAvg2(currWAvg2 < noiseLevel2).^2;
% figure, imshow(currW_NS2);
