% Detect up to 2 pronuclei / cell in entire 3D volume
%
% Inputs: cell3D is the image stack
%         numCells is the total number of cells in the embryo
%         xCenter, yCenter are center coords of each cell
%         rList is interpolated list of cell radii at each slice in jList
%         jList is the list of z-slices at which to detect nuclei
%         subplotSize1 is the number of columns in the subplot array
%
% Outputs: radList is radius of 2 strongest peaks from each slice in jList
%          psRankList is relative rank of chosen peak in terms of strength
%          psList is value of peak strength
%          areaList is sum of brightness signal inside circle for each peak
%          xcList, ycList, zcList are (x,y,z) indices of possible nuclei
%
% Livia Zarnescu
% 1-15-14



function [radList, psRankList, psList, areaList, xcList, ycList, zcList] = ...
    detectNuclei(cell3D, numCells, xCenter, yCenter, rList, jList, subplotSize1)

if numCells == 1
    peaksPerCell = 2;
else
    peaksPerCell = 1;
end

radList = [];
psRankList = [];
psList = [];
areaList = [];
xcList = [];
ycList = [];
zcList = [];

subplotSize2 = ceil(length(jList) / subplotSize1);
figure(1);
figure(2);
% areaList = zeros(numCells,length(jList));

for j = jList
    
    j
    iterNum = find(jList == j);
    % close all;
    slices = cell3D(:,:,j);
    slicesG = imfilter(slices, fspecial('gaussian', [20 20], 8));
    slicesG = imadjust(slicesG);
    sliceEdge = edge(slicesG, 'canny');%, [.005 .05]);
    
    % dilate and multiply by local intensity
    edgesDilate = imdilate(sliceEdge, strel('disk', 5));
    edgesBright = edgesDilate .* slicesG;
    
    % keep only edges with moderate intensity values
    lowThresh = graythresh(slicesG(10:end-10,10:end-10));
    edgesBright(edgesBright < lowThresh) = 0;
    
    highThresh = graythresh(edgesBright(edgesBright > 0));
    edgesBright(edgesBright > highThresh) = 0;
    
    edgesBright = imdilate(edgesBright, strel('disk', 5));
    %     figure, imshow(edgesBright);
    
    for cellNum = 1:numCells
                
        if 1 % j > firstSlice(cellNum)
            
            overlap = [];
            [iC jC] = meshgrid(1:size(edgesBright,2), 1:size(edgesBright,1));
            areaInt = [];
            slicesM = medfilt2(slicesG, [10 10]);
            
            % get rid of everything close to or outside edge of cell
            % distance from center of current cell
            currDistCenter = sqrt((iC - xCenter(cellNum)).^2 + ...
                (jC - yCenter(cellNum)).^2);
            edgesBright(currDistCenter > rList(cellNum, iterNum) - 15) = 0;
            
            % Carry out a Hough Transform to detect cell body
            % select radii for each nucleus
            radii = 20:30;
            h = circle_hough(edgesBright, radii, 'same', 'normalise');
            
            % Find some peaks in the accumulator
            nHoodXY = 31;
            nHoodR = 31;
            
            if numCells == 1
                nPeaks = 4;
            else
                nPeaks = numCells*2;
            end
            
            peaks = circle_houghpeaks(h, radii, 'nhoodxy', nHoodXY, 'nhoodr',...
                nHoodR, 'npeaks', nPeaks, 'smoothxy', 4);
            
%             peaks
            
            % Look at the results
            if cellNum == 1
                figure(1);
                hold on;
                subplot(subplotSize1,subplotSize2,iterNum);
                imagesc(edgesBright); colormap gray;
                hold on;
%                 scatter(xCenter, yCenter, 'LineWidth', 2, ...
%                     'MarkerEdgeColor', [1 0 0], ...
%                     'Marker', '+');
%                 hold on;
                
                figure(2);
                hold on;
                subplot(subplotSize1,subplotSize2,iterNum);
                imagesc(slices); colormap gray;
                hold on;
%                 scatter(xCenter, yCenter, 'LineWidth', 2, ...
%                     'MarkerEdgeColor', [1 0 0], ...
%                     'Marker', '+');
%                 hold on;
            end
            
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
                
                % if any part of the potential nucleus is outside the cell, reject it
                if distFromCenter > rList(cellNum,iterNum) - peak(3)*1.5
                    currArea = 0;
                    currAreaBright = 0;
                end
                            
                areaInt = [areaInt currArea];

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
            
            % for circles too far from center, set area to be very high
            % (i.e. lots of bright signal inside them)
            if ~isempty(areaInt)
                areaInt(areaInt == 0) = ...
                    2*max(areaInt)*ones(1,length(areaInt(areaInt == 0)));
            end
            
            % only export candidate circles with high peak strengths and
            % with the lowest signal
            [~, areaIndList] =  sort(areaInt, 'ascend');
            bestA = areaIndList(1:min(peaksPerCell,length(areaInt)));
            peakStrengths = peaks(4,bestA);
            [~, peakInds] = sort(peaks(4,:), 'descend');
            
            radToSave = peaks(3,bestA);
            psRankToSave = peakInds(bestA);
            areaToSave = areaInt(bestA);
            xcToSave = peaks(1,bestA);
            ycToSave = peaks(2,bestA);
            
            peakCutOff = 3.5;
            peaks
                        
            radList = [radList radToSave(peakStrengths > peakCutOff)];
            psRankList = [psRankList psRankToSave(peakStrengths > peakCutOff)];
            psList = [psList peakStrengths(peakStrengths > peakCutOff)];
            areaList = [areaList areaToSave(peakStrengths > peakCutOff)];
            xcList = [xcList xcToSave(peakStrengths > peakCutOff)];
            ycList = [ycList ycToSave(peakStrengths > peakCutOff)];
            
            if ~isempty(areaToSave(peakStrengths > peakCutOff))
                zcList = [zcList j*ones(1,length(areaToSave(peakStrengths > peakCutOff)))];
            end
            
%             if j > 16 && j < 30
%                 1
%             end
            
            % plot if circle found
            if max(areaInt) > 0
                
                for ii = 1:min(peaksPerCell,length(areaInt))
                    
                    if peakStrengths(ii) > peakCutOff
                        % plot
                        [x, y] = circlepoints(peaks(3,bestA(ii)));
                        x = x + peaks(1,bestA(ii));
                        y = y + peaks(2,bestA(ii));
                        
                        figure(1);
                        hold on;
                        plot(x, y, 'b-', 'LineWidth', 3);
                        title(['j = ' num2str(j)]);
                        
                        figure(2);
                        hold on;
                        plot(x, y, 'b-', 'LineWidth', 3);
                        title(['j = ' num2str(j)]);

                    end
                    
                end
                

            end

            % plot cell boundary for reference
            [x, y] = circlepoints(rList(j));
            x = x + xCenter;
            y = y + yCenter;
            
            figure(1);
            hold on;
            plot(x, y, 'b-', 'LineWidth', 3);
            
            figure(2);
            hold on;
            plot(x, y, 'b-', 'LineWidth', 3);
            
            pause(.5);
            
        else
            % Look at the results anyway, even if cells havent been
            % detected yet
            if cellNum == 1
                figure(1);
                hold on;
                subplot(subplotSize1,subplotSize2,iterNum);
                imagesc(edgesBright); colormap gray;
                hold on;
%                 scatter(xCenter, yCenter, 'LineWidth', 2, ...
%                     'MarkerEdgeColor', [1 0 0], ...
%                     'Marker', '+');
%                 hold on;
%                 
                figure(2);
                hold on;
                subplot(subplotSize1,subplotSize2,iterNum);
                imagesc(slices); colormap gray;
                hold on;
%                 scatter(xCenter, yCenter, 'LineWidth', 2, ...
%                     'MarkerEdgeColor', [1 0 0], ...
%                     'Marker', '+');
%                 hold on;
            end
        end
    end
    
end
