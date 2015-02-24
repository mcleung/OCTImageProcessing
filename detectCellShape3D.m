function [rList jListOut xcPoints ycPoints zcPoints] = ...
    detectCellShape3D(cell3D, jList, jStep, numCells, xCenter, yCenter, ...
    radiiBest, initSeed)

% next, find cell body location at each step
close all;
rList = [];
xcPoints = [];
ycPoints = [];
zcPoints = [];
jListOut = [];

% first seed is very important
if (nargin < 8)
    rLast = 35;
else
    rLast = initSeed;
end
% or can just generally detect all circles
useSeed = 1;

for j = jList

    j
    % close all;
    slices = cell3D(:,:,j);
    slicesG = imfilter(slices, fspecial('gaussian', [10 10], 4));
    slicesG = imadjust(slicesG);
%     im = medfiltL(slicesG, 15);
%     im2 = imadjust(im);
%     figure, imshow(im2);
    sliceEdge = edge(slicesG, 'canny');%, [.005 .05]);
%     figure, imshow(sliceEdge);
    
    % dilate and multiply by local intensity
    edgesDilate = imdilate(sliceEdge, strel('disk', 5));
    edgesBright = edgesDilate .* slicesG;
    
    % keep only edges with moderate intensity values
    lowThresh = graythresh(slicesG(10:end-10,10:end-10));
    edgesBright(edgesBright < lowThresh) = 0;
    
    highThresh = graythresh(edgesBright(edgesBright > 0));
    edgesBright(edgesBright > highThresh) = 0;
%     figure, imshow(edgesBright);
    
    if size(cell3D,3) < 60
        multFactor = 8; %8
    else
        multFactor = 4;
    end

    % Carry out a Hough Transform to detect cell body 
    if useSeed && j ~= jList(1)
        if j < size(cell3D,3) / 2 - multFactor
            radii = rLast - multFactor/2 : min(rLast + jStep*multFactor,105);
        elseif j > size(cell3D,3) / 2 + multFactor
            radii = rLast - multFactor : min(rLast + multFactor,105);
        else
            radii = min(radiiBest - multFactor/2,103) : min(radiiBest + multFactor/2,105);
        end
    else
        radii = 30:105;
    end
    
    h = circle_hough(edgesBright, radii, 'same', 'normalise');
    
    % Find some peaks in the accumulator
    nHoodXY=51;
    nHoodR=51;
    nPeaks=numCells*2;
    peaks = circle_houghpeaks(h, radii, 'nhoodxy', nHoodXY, 'nhoodr',...
        nHoodR, 'npeaks', nPeaks, 'smoothxy', 4);
    
    % Look at the results
    clf; 
    imagesc(slices); colormap gray;
    hold on;
    overlap = [];
    
    for i = 1:size(peaks,2)
        
        peak = peaks(:,i);
        
        [x, y] = circlepoints(peak(3));
        x = x+peak(1);
        y = y+peak(2);
%         plot(x, y, 'b-');
%         hold on;
        
        outOfBounds = (x < 1) | (x > size(edgesBright,2)) | (y < 1) | ...
            (y > size(edgesBright,1));
        
        x(outOfBounds) = NaN;
        y(outOfBounds) = NaN;
        
        x = x(~isnan(x));
        y = y(~isnan(y));
        currOverlap = 0;
        
        for k = 1:length(x)     
            currOverlap = currOverlap + edgesBright(y(k),x(k));           
        end        
        %     currOverlap = currOverlap;% / peak(3);  
        
        peakLocationX = min(abs(peak(1) - xCenter));
        peakLocationY = min(abs(peak(2) - yCenter));
         
        if peakLocationX > 20 || peakLocationY > 20
            currOverlap = 0;
        end
        
        overlap = [overlap currOverlap];
        %     pause;
        
    end
    
    bestIO = find(overlap == max(overlap))
    
    % plot if circle found
    if overlap(bestIO) > 0
        [x, y] = circlepoints(peaks(3,bestIO));
        x = x + peaks(1,bestIO);
        y = y + peaks(2,bestIO);
        
        % find best fit ellipse based on best fit circle
        %     findBestEllipse(x, y, im2, sliceEdge);
        
        hold on;
        plot(x, y, 'g-', 'LineWidth', 3);
        rLast = peaks(3,bestIO)
        
        rList = [rList peaks(3,bestIO)];
        xcPoints = [xcPoints x];
        ycPoints = [ycPoints y];
        zcPoints = [zcPoints j*ones(1, length(x))];
        jListOut = [jListOut j];
        
    end
    
    pause(1);
    
end