% find centers of cells given number of cells

function peaks = findCellCentersByProjection(projectionNorm, slice2D, cellStage, numCells)

edges = edge(projectionNorm, 'canny');
% figure, imshow(edges);

edgesDilate = imdilate(edges, strel('disk', 5));
edgesBright = edgesDilate .* projectionNorm;

% keep only edges with moderate intensity values
lowThresh = graythresh(projectionNorm);
edgesBright(edgesBright < lowThresh) = 0;

highThresh = graythresh(edgesBright(edgesBright > 0));
edgesBright(edgesBright > highThresh) = 0;
% figure, imshow(edgesBright);

if cellStage == 1
    radii = 85:105;
elseif cellStage == 2
    radii = 65:85;
elseif cellStage == 3
    radii = 24:30;
end

edgesBright(edgesBright > 0) = 1;
h = circle_hough(edgesBright, radii, 'same', 'normalise');

% Find some peaks in the accumulator
nHoodXY=31;
nHoodR=31;
nPeaks=numCells;
peaks = circle_houghpeaks(h, radii, 'nhoodxy', nHoodXY, 'nhoodr',...
    nHoodR, 'npeaks', nPeaks, 'smoothxy', 4);

% Look at the results
figure, imshow(imadjust(slice2D, [.239 .6353], [0 .8])); colormap gray;
hold on;

% calculate average overlap of each circle with edges
overlap = [];

for i = 1:size(peaks,2)
    
    peak = peaks(:,i);
    [x, y] = circlepoints(peak(3));
    x = x+peak(1);
    y = y+peak(2);
    plot(x, y, 'b-', 'LineWidth', 3);
    hold on;
    
    outOfBounds = (x < 1) | (x > size(projectionNorm,2)) | (y < 1) | ...
        (y > size(projectionNorm,1));
    
    x(outOfBounds) = NaN;
    y(outOfBounds) = NaN;
    
    x = x(~isnan(x));
    y = y(~isnan(y));
    currOverlap = 0;
    
    for k = 1:length(x)
        %             currOverlap = currOverlap + im(y(k),x(k))*edgeLabel3(y(k),x(k));
        currOverlap = currOverlap + edgesBright(y(k),x(k));
    end
    %     currOverlap = currOverlap;% / peak(3);
    overlap = [overlap currOverlap];

end