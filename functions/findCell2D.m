function [ peaks ] = findCell2D( im,nCells)
debug=1;
%FINDCELL2D Summary of this function goes here
%   Detailed explanation goes here
imO=(im-min(min(im)))/(max(max(im))-min(min(im))); % scale
im= medfilt2(imO, [25 25]); % median filter
% threshold
level=graythresh(im);
BW=zeros(size(im)); 
BW(find(im>level))=im(find(im>level)); 
BW=imcomplement(BW); %
im=BW;
% canny edge detection
eZona = edge(im, 'canny',0.6);
if debug==1 imshow(eZona); end;
eAll = edge(im, 'canny',0.1);
if debug==1 imshow(eAll); end;
eWoZona=eAll-eZona;
if debug==1 imshow(eWoZona);end;
eWoZona=eAll;
if debug ==1 imshow(eWoZona); end;
%% Carry out the HT
% The circles round the cells have radii in the 100-150 pixels range. 
% We select the 'same' option to simplify later processing, and the
% 'normalise' option to avoid a bias towards finding larger circles.
radii = 100/2:1:250/2;
h = circle_hough(eWoZona, radii, 'same', 'normalise');
%% Find some peaks in the accumulator
% We use the neighbourhood-suppression method of peak finding to ensure
% that we find spatially separated circles. We select the 10 most prominent
% peaks, because as it happens we can see that there are 10 coins to find.
nHoodXY=50+1;
nHoodR=50+1;
nPeaks=2*nCells;
peaks = circle_houghpeaks(h, radii, 'nhoodxy', nHoodXY, 'nhoodr',...
        nHoodR, 'npeaks', nPeaks);

%% Look at the results
if debug==1
    % We draw the circles found on the image, using both the positions and the
    % radii stored in the |peaks| array. The |circlepoints| function is
    % convenient for this - italso used by |circle_hough| so comes with it.
    imshow(imO,[]);
    %chooseCell;
    hold on;
    for peak = peaks
        [x, y] = circlepoints(peak(3));
        plot(x+peak(1), y+peak(2), 'g-');
    end
    hold off
end
end

