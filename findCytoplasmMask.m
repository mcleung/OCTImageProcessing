% This function outputs a 3D mask giving the cytoplasm area to analyze
%
% Usage: [cellMask, minPNslice, maxPNslice] = findCytoplasmMask(cell3D, ...
%   Group1, Group2, groupNum, fileName, XYpixelsPerSlice, micronsPerSlice)

function [cellMask, minPNslice, maxPNslice] = findCytoplasmMask(cell3D, ...
    Group1, Group2, groupNum, fileName, XYpixelsPerSlice, micronsPerSlice)

cellMask = 0*cell3D;
cellMask = (cellMask == 1);

embryoStruct = eval(['Group' num2str(groupNum) '.E' fileName]);
rCell = embryoStruct.cellBody.r;

% zc is in slices
% r is in pixels
PN1zc = embryoStruct.PN1.zc;
PN1r = embryoStruct.PN1.r + 5;
PN2zc = embryoStruct.PN2.zc;
PN2r = embryoStruct.PN2.r + 5;

% find bottom of lower PN, and top of upper PN, units are in slices
if (PN1zc < PN2zc)
    minPNslice = round(PN1zc - PN1r/XYpixelsPerSlice - 5/micronsPerSlice);
    maxPNslice = round(PN2zc + PN2r/XYpixelsPerSlice + 5/micronsPerSlice);
else
    minPNslice = round(PN2zc - PN2r/XYpixelsPerSlice - 5/micronsPerSlice);
    maxPNslice = round(PN1zc + PN1r/XYpixelsPerSlice + 5/micronsPerSlice);
end

minPNslice = max(3,minPNslice);
maxPNslice = min(maxPNslice, size(cell3D,3)-2);

% we already have the radius of the cell body as a function of z
% now calculate expected size of each PN as a function of z
rList1 = zeros(1,size(cell3D,3));
rList2 = zeros(1,size(cell3D,3));

% units are in pixels, not slices
for i = 1:size(cell3D,3)
    
    % set values for PN1
    if abs(PN1zc - i) < round(PN1r/XYpixelsPerSlice)
        rList1(i) = ceil(sqrt((PN1r/XYpixelsPerSlice)^2 - (PN1zc - i)^2));
    end
    
    % set values for PN2
    if abs(PN2zc - i) < round(PN2r/XYpixelsPerSlice)
        rList2(i) = ceil(sqrt((PN2r/XYpixelsPerSlice)^2 - (PN2zc - i)^2));
    end
    
end

rList1 = rList1*XYpixelsPerSlice;
rList2 = rList2*XYpixelsPerSlice;

[iC jC] = meshgrid(1:size(cell3D,2), 1:size(cell3D,1));        

% D0 is distance form cell center
% D1 is distance from PN1
% D2 is distance from PN2
D0 = sqrt((iC - embryoStruct.cellBody.xc).^2 + ...
    (jC - embryoStruct.cellBody.yc).^2);
D1 = sqrt((iC - embryoStruct.PN1.xc).^2 + ...
    (jC - embryoStruct.PN1.yc).^2);
D2 = sqrt((iC - embryoStruct.PN2.xc).^2 + ...
    (jC - embryoStruct.PN2.yc).^2);

% figure(1);

% finally, go thru all slices in cell3D and make mask
for i = 1:size(cell3D,3)
    
    % only add true values to mask if between minPNslice and maxPNslice
    if i > minPNslice && i < maxPNslice
        
        % for curr slice, keep pixels inside rCell + 5, but not pixels
        % inside each PN

        currMask = zeros(size(cell3D,1), size(cell3D,2));
        currMask((D0 < (rCell(i)+5)) & (D1 > (rList1(i))) ...
            & (D2 > (rList2(i)))) = 1;
        
        cellMask(:,:,i) = (currMask == 1);
%         imshow(currMask);
%         pause(1);
        
    end
    
end










