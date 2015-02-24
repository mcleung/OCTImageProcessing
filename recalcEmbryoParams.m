% recalculate embryo parameters based on already processed .mat files

%% Step 1: load parameters

%clear all;

inFile = 'embryoParams4.mat';
load(inFile);

groupNumList = [2];
eNumList = [8];


for ii = 1:length(groupNumList)
    
    close all;
    groupNum = groupNumList(ii);
    fileName = ['E' num2str(eNumList(ii))];
    micronsPerSlice = 1;
    XYpixelsPerSlice = 4.5*micronsPerSlice;
    
    % Step 2: extract cell mask, etc
    % only non-NaN pixels are in cytoplasm
    
    cell3D = eval(['Group' num2str(groupNum) '.' fileName '.cell3D;']);
    
    [cellMask, minPNslice, maxPNslice] = findCytoplasmMask(cell3D, Group1, ...
        Group2, groupNum, fileName, XYpixelsPerSlice, micronsPerSlice);
    
    round(minPNslice/2 + maxPNslice/2)
    figure, imshow(cellMask(:,:,round(minPNslice/2 + maxPNslice/2)) .* ...
        cell3D(:,:,round(minPNslice/2 + maxPNslice/2)));
    
    % Step 3: Extract "clumpiness" parameters
    
    % right now just extract entropy and stdev for cytoplasm pixels in each
    % slice
    rFit = eval(['Group' num2str(groupNum) '.' fileName '.cellBody.r;']);
    xCenter = eval(['Group' num2str(groupNum) '.' fileName '.cellBody.xc;']);
    yCenter = eval(['Group' num2str(groupNum) '.' fileName '.cellBody.yc;']);
    
    paramList = extractCytoplasmParams(cell3D, cellMask, minPNslice, ...
        maxPNslice, rFit, xCenter, yCenter);
    
    % Step 4: Save new params
    
    eval(['Group' num2str(groupNum) '.' fileName '.stdev = paramList.stdev;']);
    eval(['Group' num2str(groupNum) '.' fileName '.entropy = paramList.entropy;']);
    eval(['Group' num2str(groupNum) '.' fileName '.entropyfilt = paramList.entropyfilt;']);
    eval(['Group' num2str(groupNum) '.' fileName '.homogeneity = paramList.homogeneity;']);
    eval(['Group' num2str(groupNum) '.' fileName '.energy = paramList.energy;']);
    eval(['Group' num2str(groupNum) '.' fileName '.rProfile = paramList.rProfile;']);
    eval(['Group' num2str(groupNum) '.' fileName '.sliceList = minPNslice+1:maxPNslice-1;']);
    
    
end

save(inFile, 'Group1', 'Group2');


