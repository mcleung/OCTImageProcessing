function dataAdj = AdjCBFilt(dataBackSub, coords)

xmin = coords(1);
xmax = coords(2);
ymin = coords(3);
ymax = coords(4);

dataAdj = dataBackSub;

% find noise levels in rectangle with noise only
% this is the reference frame to scale all other noise to
refFrame = dataBackSub(ymin:ymax, xmin:xmax, ceil((size(dataBackSub,3))/2));

% for each slice, adjust the contrast and filter it
for i =1:size(dataBackSub,3)    
    %i
    currSlice = dataBackSub(:,:,i);
%     figure, imshow(imadjust(currSlice));
    noiseFrame = currSlice(ymin:ymax, xmin:xmax);

    % scale image so its noise level is the same as that of the reference
    % frame
    currSliceAdj = imadjust(currSlice, [min(noiseFrame(:)) ...
        max(noiseFrame(:))], [min(refFrame(:)) max(refFrame(:))]);
    
    % wiener filter to reduce noise and adjust window again
    currW = wiener2(currSliceAdj, [5 5]);
    currWAdj = imadjust(currW, [0 min(max(currW(:)),1)], [0 1]);
%     figure, imshow(currWAdj);
    
    % filter out noise (pixels in low density areas)
    currW_NS = filterNoise(currWAdj, coords);
    dataAdj(:,:,i) = currW_NS;
    %if mod(i,100)==0 waitbar(i/size(dataBackSub,3));end;
end