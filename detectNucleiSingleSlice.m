% detect nuclei

function [areaList overlapList] = detectNucleiSingleSlice(cell3D, numCells, ...
    xCenter, yCenter, rList)

areaList = zeros(numCells,1);
overlapList = zeros(numCells,1);

% close all;
slices = mean(cell3D,3);
slicesG = imfilter(slices, fspecial('gaussian', [10 10], 4));
slicesG = imadjust(slicesG);

% dilated edges around graythresh boundaries
slicesBin = slicesG < graythresh(slicesG);
slicesBinD = imdilate(slicesBin, strel('disk', 7));
slicesBinE = imerode(slicesBin, strel('disk', 7));
slicesBin = slicesBinD - slicesBinE;

% multiplier of distance from xCenter, yCenter
[xMat, yMat] = meshgrid(1:size(cell3D,2), 1:size(cell3D,1));
rMat = sqrt((xMat - xCenter).^2 + (yMat - yCenter).^2);
rListScaled = 1 - (rList*.9 - min(min(rMat)))/(max(max(rMat)) - min(min(rMat)));
rMat = (rMat - min(min(rMat)))/(max(max(rMat)) - min(min(rMat)));
rMat = 1 - rMat;
rMat(rMat < rListScaled) = 0;

% dilate and multiply by local intensity
edgesBright = rMat .* slicesBin; 

for cellNum = 1:numCells
    
    % Carry out a Hough Transform to detect cell body
    % select radii for each nucleus
    radii = 18:25;
    h = circle_hough(edgesBright, radii, 'same');%, 'normalise');
    
    % Find some peaks in the accumulator
    nHoodXY = 31;
    nHoodR = 31;
    nPeaks = 10;
    peaks = circle_houghpeaks(h, radii, 'nhoodxy', nHoodXY, 'nhoodr',...
        nHoodR, 'npeaks', nPeaks, 'smoothxy', 4);
    
    overlap = [];
    [iC jC] = meshgrid(1:size(edgesBright,2), 1:size(edgesBright,1));
    areaInt = [];
    slicesM = medfilt2(slicesG, [10 10]);
    
    % go through peaks detected from hough transform
    for i = 1:size(peaks,2)
        
        peak = peaks(:,i);
        
        [x, y] = circlepoints(peak(3));
        x = x+peak(1);
        y = y+peak(2);
        outOfBounds = (x < 1) | (x > size(edgesBright,2)) | (y < 1) | ...
            (y > size(edgesBright,1));
        
        x(outOfBounds) = NaN;
        y(outOfBounds) = NaN;
        
        x = x(~isnan(x));
        y = y(~isnan(y));
        
        % find intensity inside circle
        currDist = sqrt((iC - peak(1)).^2 + (jC - peak(2)).^2);
        currArea = sum(sum(slicesM(currDist < peak(3))/(pi*peak(3)^2)));
        
        % don't want to detect circles too far from cell center in
        % xy
        distFromCenter = sqrt((xCenter(cellNum) - peak(1))^2 + ...
            (yCenter(cellNum) - peak(2))^2);
        
%         % if it's outside the cell, reject it
%         if distFromCenter > rList(cellNum)/2
%             currArea = 0;
%         end
        
        areaInt = [areaInt currArea*distFromCenter];
        
        % now find how well the circles overlap edges
        currOverlap = 0;
        
        for k = 1:length(x)
            currOverlap = currOverlap + ...
                edgesBright(y(k),x(k));%*slices(y(k),x(k));
        end
        
        % normalize for circumference
        currOverlap = currOverlap / peak(3);
        overlap = [overlap currOverlap];
        
    end
    
    cellNum
    areaInt(areaInt == 0) = ...
        max(areaInt)*ones(1,length(areaInt(areaInt == 0)));
    
    bestIO = find(overlap == max(overlap))
    bestA = find(areaInt == min(areaInt))
    areaInt
    
    [~, IOinds] = sort(overlap, 'descend');
    [~, Ainds] = sort(areaInt, 'ascend');
    
    areaList = areaInt(Ainds(1:2));
    overlapList = overlap(IOinds(1:2));
    
    % plot if circle found
    for i = 1:2
        if overlapList(i) > 0 && areaList(i) > 0
            
            [x, y] = circlepoints(peaks(3,Ainds(i)));
            x = x + peaks(1,Ainds(i));
            y = y + peaks(2,Ainds(i));
            
%             hold on;
%             plot(x, y, 'g-', 'LineWidth', 3);
%             pause;
            
            % see how well overlap agrees with area condition
            [x, y] = circlepoints(peaks(3,IOinds(i)));
            x = x + peaks(1,IOinds(i));
            y = y + peaks(2,IOinds(i));
            
            hold on;
            plot(x, y, 'g-', 'LineWidth', 3);
            pause;
            
        end
    end
        
end


