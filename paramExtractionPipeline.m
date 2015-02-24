% Parameter Extraction Pipeline for OCT images
% Livia Zarnescu
% 1-14-14

% Input: images of 2PN mouse embryos
% Steps: 1) Detect cell edges in middle 50% of cell
%        2) Detect pronuclei locations
%        3) Take slices from z-locations with pronuclei +/- 10 microns
%        4) Calculate "clumpiness" of all cytoplasm area
%
% Could be later converted into a function to process large numbers of
% image files

%% Step 1.1: Read in desired tif stack

clear all;
close all;

% add things to path
addpath('C:\Users\Livia\Desktop\CommonMatlabFiles');
addPathRecursive('.');

micronsPerSlice = 1;
XYpixelsPerSlice = 4.5*micronsPerSlice;

% *************** CHANGE THIS IF NECESSARY *****************************
filePath = 'C:\Users\Livia\Desktop\OCT\OCT data\SPIE_paper_data\group 2 (2PN vit)\GaussBlur 1-23-15';
fileName = 'E11_14'; % 1,5,8,10
embryoNum = num2str(12); 
maxNoVit = 14;
maxVit = 14;
% *********************************************************************

groupNum = filePath(strfind(filePath, '\group ') + 7);
currDate = filePath(strfind(filePath, '\GaussBlur ') + 11 : end);

s = strfind(filePath, ' ');
exptDate = filePath(s(end)+1:end);

% determine file to save in
if (strcmp(exptDate, '3-27-14'))
    fileToSave = 'embryoParams4.mat';
elseif (strcmp(exptDate, '5-29-14'))
    fileToSave = 'embryoParams5.mat';
elseif (strcmp(exptDate, '1-14-14'))
    fileToSave = 'embryoParams.mat';
elseif (strcmp(exptDate, '12-4-13'))
    fileToSave = 'embryoParams2.mat';
elseif (strcmp(exptDate, '9-25-14'))
    fileToSave = 'embryoParams6.mat';
elseif (strcmp(exptDate, '8-14-14'))
    fileToSave = 'embryoParams7.mat';
elseif (strcmp(exptDate, '1-23-15'))
    fileToSave = 'embryoParams8.mat';
else
    fileToSave = 'embryoParams3.mat';
end

dataPath = [filePath '\' fileName '.tif'];
%dataPath = [filePath '\E3.tif'];

imInfo = imfinfo(dataPath);
numImages = length(imInfo);
data = zeros(imInfo(1).Height, imInfo(1).Width, numImages);

% read in images
for i = 1:numImages
    data(:,:,i) = im2double(imread(dataPath, i));
end

% Step 1.1: Select box around a single cell (this will be automated later)

middleHeight = round(numImages/2);
figure(1);
imshow(data(:,:,middleHeight));
coord = getrect; % [xmin ymin width height]
coords = round([coord(1) coord(1)+coord(3) coord(2) coord(2)+coord(4)]);
cell3D = data(coords(3):coords(4), coords(1):coords(2), :);
close(1);

% Step 1.2: Detect cell center from 2D projection

numCells = 1;
cellStage = 1;

projection = sum(cell3D,3);
projectionNorm = (projection - min(min(projection))) / ...
    (max(max(projection)) - min(min(projection)));

peaks = findCellCentersByProjection(projectionNorm, ...
    cell3D(:,:,middleHeight), cellStage, numCells);
xCenter = peaks(1,:);
yCenter = peaks(2,:);
radiiBest = peaks(3,:);

% Step 1.3 Detect 3D shape of cell

jStep = 1; % jstep = 2 for stacks with 130 slices
jList = 5:jStep:(numImages-5); % 10:jStep:(numImages-10) for stacks with 130 slices
[rList, jListOut, ~, ~, ~] = detectCellShape3D(cell3D, jList, jStep, ...
    numCells, xCenter, yCenter, radiiBest, 75);

% rFit is the "expected" cell radius as a function of z over all slices
figure, plot(jListOut, rList);
p = polyfit(jListOut, rList,2);
rFit = round(polyval(p, 1:numImages));
rFit2 = round(polyval(p, jListOut));
hold on; plot(1:numImages, rFit, 'color', [1 0 0]);
rList(abs(rList - rFit2) > 5) = rFit2(abs(rList - rFit2) > 5);

% redo polynomial fit
figure, plot(jListOut, rList);
p = polyfit(jListOut, rList,2);
rFit = round(polyval(p, 1:numImages));
hold on; plot(1:numImages, rFit, 'color', [1 0 0]);

%% Step 2: Detect pronuclei in 3D

% rFit(jList) is the expected cell radius at each slice in jList

subplotSize1 = 5;

[radList, psRankList, psList, areaList, xcList, ycList, zcList] = ...
    detectNuclei(cell3D, numCells, ...
    xCenter, yCenter, rFit, jListOut, subplotSize1);

% p(1,1) and p(2,1) should be locations of pronuclei in xcList, ycList,
% zList, areaList, radList, etc.
[p h] = mspeaks(1:length(areaList), 1-areaList)

% figure, plot(1:length(areaList), 1-areaList);
% figure, plot(zcList, 1-areaList);

if size(p,1) == 1
    correctSlice = zcList(p(1));
    
    % if there is only one peak, that means both PNs are in same slice
    % find all circles from that z slice
    otherInd = find(zcList == correctSlice);
    
    % take just those 2 circles in case there are other circles found 
    areaOtherInd = 1 - areaList(otherInd); 
    [~ , otherIndSort] = sort(areaOtherInd, 'descend');
    p = repmat(p,2,1);
    p(1,1) = otherInd(otherIndSort(1));
    p(2,1) = otherInd(otherIndSort(2));
end

p = round(p);
xcList
ycList
zcList

% % manually change p if need be
% p = [8 .3; ...
%     12 .3];

% p = [9;16]

% use mspeaks to extract radii, xc, and yc of two strongest peaks
% verify manually

%% Step 2.1 : Load and save embryo data

if exist(fileToSave, 'file')
    load(fileToSave);
else
    
    Group1 = struct();
    
    for i = 1:maxNoVit
        fieldName = ['E' num2str(i)]
        Group1.(fieldName) = struct();
    end
    
    Group2 = struct();
    
    for i = 1:maxVit
        fieldName = ['E' num2str(i)]
        Group2.(fieldName) = struct();
    end
    
end

% set date of current embryo
% set (x,y,z,r) of pronuclei
% set (xc,yc,z,r) of cell body, where z is a list of z coords and r = r(z)
% (xc, yc) are in cartesian coords - to access matrix entry, use (yc, xc)
eval(['Group' num2str(groupNum) '.E' embryoNum '.date = currDate;']);
eval(['Group' num2str(groupNum) '.E' embryoNum '.PN1 = struct();']);
eval(['Group' num2str(groupNum) '.E' embryoNum '.PN2 = struct();']);
eval(['Group' num2str(groupNum) '.E' embryoNum '.cellBody = struct();']);

eval(['Group' num2str(groupNum) '.E' embryoNum '.PN1.xc = xcList(p(1,1));']);
eval(['Group' num2str(groupNum) '.E' embryoNum '.PN1.yc = ycList(p(1,1));']);
eval(['Group' num2str(groupNum) '.E' embryoNum '.PN1.zc = zcList(p(1,1));']);
eval(['Group' num2str(groupNum) '.E' embryoNum '.PN1.r = radList(p(1,1));']);

eval(['Group' num2str(groupNum) '.E' embryoNum '.PN2.xc = xcList(p(2,1));']);
eval(['Group' num2str(groupNum) '.E' embryoNum '.PN2.yc = ycList(p(2,1));']);
eval(['Group' num2str(groupNum) '.E' embryoNum '.PN2.zc = zcList(p(2,1));']);
eval(['Group' num2str(groupNum) '.E' embryoNum '.PN2.r = radList(p(2,1));']);

eval(['Group' num2str(groupNum) '.E' embryoNum '.cellBody.xc = xCenter;']);
eval(['Group' num2str(groupNum) '.E' embryoNum '.cellBody.yc = yCenter;']);
eval(['Group' num2str(groupNum) '.E' embryoNum '.cellBody.z = 1:numImages;']);
eval(['Group' num2str(groupNum) '.E' embryoNum '.cellBody.r = rFit;']);

save(fileToSave, 'Group1', 'Group2');

% Step 2.2 Return 3D mask of pixels to analyze for each slice

[cellMask, minPNslice, maxPNslice] = findCytoplasmMask(cell3D, Group1, ...
    Group2, groupNum, embryoNum, XYpixelsPerSlice, micronsPerSlice);


% Step 4: Extract "clumpiness" parameters

% right now just extract entropy and stdev for cytoplasm pixels in each
% slice
paramList = extractCytoplasmParams(cell3D, cellMask, minPNslice, ...
    maxPNslice, rFit, xCenter, yCenter);
% rProfile = extractRadialProfile(cell3D, cellMask, minPNslice, maxPNslice, rList, xCenter, yCenter);

eval(['Group' num2str(groupNum) '.E' embryoNum '.stdev = paramList.stdev;']);
eval(['Group' num2str(groupNum) '.E' embryoNum '.entropy = paramList.entropy;']);
eval(['Group' num2str(groupNum) '.E' embryoNum '.entropyfilt = paramList.entropyfilt;']);
eval(['Group' num2str(groupNum) '.E' embryoNum '.homogeneity = paramList.homogeneity;']);
eval(['Group' num2str(groupNum) '.E' embryoNum '.energy = paramList.energy;']);
eval(['Group' num2str(groupNum) '.E' embryoNum '.rProfile = paramList.rProfile;']);
eval(['Group' num2str(groupNum) '.E' embryoNum '.sliceList = minPNslice+1:maxPNslice-1;']);

% also save cell3D, cellMask in case I want to re-calc parameters
eval(['Group' num2str(groupNum) '.E' embryoNum '.cell3D = cell3D;']);
eval(['Group' num2str(groupNum) '.E' embryoNum '.cellMask = cellMask;']);

save(fileToSave, 'Group1', 'Group2');

close all;
% figure, plot(eval(['Group' num2str(groupNum) '.' fileName '.sliceList']), ...
%     eval(['100*Group' num2str(groupNum) '.' fileName '.stdev']));
% hold on; 
% plot(eval(['Group' num2str(groupNum) '.' fileName '.sliceList']), ...
%     eval(['Group' num2str(groupNum) '.' fileName '.entropy']), 'r');

round(minPNslice/2 + maxPNslice/2)
figure, imshow(cellMask(:,:,round(minPNslice/2 + maxPNslice/2)) .* ...
    cell3D(:,:,round(minPNslice/2 + maxPNslice/2)));







