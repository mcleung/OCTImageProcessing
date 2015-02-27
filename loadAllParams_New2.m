
% load params from all dates

function [paramsOut, groundTruthOut] = loadAllParams_New2(exptsToLoad, varToPlot, posClass)

if nargin < 2
    varToPlot = 'blastForm';
    posClass = 1;
end

% initialize all output parameters
blastFormOut = [];
cellClumpOut = [];
vitGroupOut = [];

stdevList = [];
diffStart = [];
diffMid = [];
diffEnd = [];
diffdiff1 = [];
diffdiff2 = [];
darkPercent = [];
pnDistFromCenter1 = [];
pnDistFromCenter2 = [];
pnSeparation = [];

rProfileOut = [];

% store actual ground truth data in cell arrays
vitGroup = cell(1,length(exptsToLoad));
cellClump = cell(1,length(exptsToLoad));
blastForm = cell(1,length(exptsToLoad));
fileNameAssoc = cell(1,length(exptsToLoad));

% 1 = vit, 0 = non-vit
vitGroup{1} = [zeros(1,4) ones(1,6)]; % 4 non-vit from 4-26-13 and 6 vit from 5-3-13
vitGroup{2} = [zeros(1,5) ones(1,5)]; % 5 non-vit and 5 vit from 12-4-13
vitGroup{3} = [zeros(1,5) ones(1,5)]; % 5 non-vit and 5 vit from 1-14-14
vitGroup{4} = [zeros(1,5) ones(1,10)]; % 5 non-vit and 10 vit from 3-27-14
vitGroup{5} = [zeros(1,5) ones(1,10)]; % 5 non-vit and 10 vit from 5-29-14
vitGroup{6} = [zeros(1,5) ones(1,12)]; % 5 non-vit and 12 vit from 9-25-14
vitGroup{7} = [zeros(1,5) ones(1,5)]; % 5 non-vit and 5 vit from 8-14-14
vitGroup{8} = [zeros(1,14) ones(1,14)]; % 1-14 non-vit (fresh, not slow frozen), and 1-14 post-vit
vitGroup{9} = [zeros(1,10) ones(1,10)]; % 1-10 non-vit (fresh, not slow frozen), and 1-10 post-vit

blastForm{1} = NaN*ones(1,10);
blastForm{2} = NaN*ones(1,10);
blastForm{3} = [0 0 1 1 1 0 0 0 1 1];
blastForm{4} = floor([1 1 0 1 1 1 1 0 0 .6 1 0 0 .6 0]); % .6 means very early blast or unsure
blastForm{5} = floor([1 1 1 0 1 .6 1 1 1 0 .6 0 .6 1 1]); % .6 means very early blast or unsure
blastForm{6} = floor([1 0 0 0 1 1 1 1 0 1 0 1 1 0 .6 1 0]); % .6 means very early blast or unsure
blastForm{7} = NaN*ones(1,10); % used for fluorescence imaging
blastForm{8} = NaN*ones(1,28); % frozen down for more vit cycles
blastForm{9} = NaN*ones(1,20);

% manual evaluation of clumping
% 1 = clump, 0 = no clump, 2 = can't tell
cellClump{1} = [0 0 0 0 1 1 1 1 1 1]; % 1-10
cellClump{2} = [0 0 0 0 1 1 0 1 1 1]; % 11-20
cellClump{3} = [1 0 1 0 0 1 1 1 1 1]; % 21-30
cellClump{4} = [1 0 0 0 1 1 1 1 1 1 0 1 1 1 0]; % 31-45 
cellClump{5} = [1 1 0 1 1 0 1 1 1 1 1 0 1 1 1]; % 46-60
cellClump{6} = [1 0 1 0 0 0 1 1 1 0 1 0 1 1 1 1 1]; % 61-77 
cellClump{7} = [1 0 1 0 0 0 1 0 1 1]; % 78-87 
cellClump{8} = [zeros(1,14) ones(1,14)]; % 88-101, and 102-115
cellClump{9} = [zeros(1,10) ones(1,10)];

fileNameAssoc{1} = 'embryoParams3.mat'; % 4-26-13 and 5-3-13
fileNameAssoc{2} = 'embryoParams2.mat'; % 12-4-13
fileNameAssoc{3} = 'embryoParams.mat'; % 1-14-14
fileNameAssoc{4} = 'embryoParams4.mat'; % 3-27-14
fileNameAssoc{5} = 'embryoParams5.mat'; % 5-29-14
fileNameAssoc{6} = 'embryoParams6.mat'; % 9-25-14
fileNameAssoc{7} = 'embryoParams7.mat'; % 8-14-14
fileNameAssoc{8} = 'embryoParams8.mat'; % 1-23-15
fileNameAssoc{9} = 'eParam1.mat'; % 2015-02-20

dh7Mat = [];

for i = 1:length(exptsToLoad)
    
    if (exptsToLoad(i) && exist(fileNameAssoc{i}, 'file'))
        
        load(fileNameAssoc{i});
        numNonVit = sum(1 - vitGroup{i});
        
        for j = 1:length(vitGroup{i})
            
            j
            
            if j < (numNonVit + 1)
                groupNum = 0;
                embryoNum = j;
            else
                groupNum = 1;
                embryoNum = j - numNonVit;
            end
            
            if eval([varToPlot '{i}(j) == posClass'])
                if groupNum == 1
                    currColor = [0 0 .6]; % non-vit in posClass (dark blue)
                else
                    currColor = [0 .6 .6]; % vit in posClass (light blue)
                end
            else
                if groupNum == 1
                    currColor = [1 0 0]; % non-vit in negative class (red)
                else
                    currColor = [.7 .7 0]; % vit in negative class (yellow)
                end
            end
    
            if isfield(eval(['Group' num2str(groupNum) ...
                    '.E' num2str(embryoNum)]), 'cellMask')
                
                
                % add ground truth data to corresponding output arrays
                vitGroupOut = [vitGroupOut vitGroup{i}(j)];
                cellClumpOut = [cellClumpOut cellClump{i}(j)];
                blastFormOut = [blastFormOut blastForm{i}(j)];
                
                figure(1);
                subplot(2,3,1);
                
                cellMask = eval(['Group' num2str(groupNum) '.E' num2str(embryoNum) ...
                    '.cellMask']);
                cell3D = eval(['Group' num2str(groupNum) '.E' num2str(embryoNum) ...
                    '.cell3D']);
                numZ = size(cell3D,3);
                imshow(cellMask(:,:,round(numZ/2)) .* cell3D(:,:,round(numZ/2)));
                
                
                % get proportion of dark pixels in cytoplasm
                
                numXY = sum(sum(cellMask(:,:,round(numZ/2))));
                cytoplasmPixels = cellMask(:,:,round(numZ/2)) .* cell3D(:,:,round(numZ/2));
                meanBright = mean(cytoplasmPixels(cytoplasmPixels > 0));
                stdBright = std(cytoplasmPixels(cytoplasmPixels > 0));
                darkPercent = [darkPercent sum(sum(cytoplasmPixels > 0 & ...
                    cytoplasmPixels < (meanBright - stdBright)))/numXY];
                
                %             cellClump{i}(j)
                %             mean(cytoplasmPixels(cytoplasmPixels > 0))
                %             std(cytoplasmPixels(cytoplasmPixels > 0))
                %             sum(sum(cytoplasmPixels > 0 & cytoplasmPixels < (meanBright - .05)))/numXY
                %             figure(3);
                %             imshow(cytoplasmPixels);
                %
                %             figure(4);
                %             imhist(cytoplasmPixels);
                
                stdevCurr = eval(['Group' num2str(groupNum) '.E' ...
                    num2str(embryoNum) '.stdev;']);
                
                stdevList = [stdevList mean(stdevCurr)];
                
                
                h1 = eval(['Group' num2str(groupNum) '.E' num2str(embryoNum) ...
                    '.rProfile']);
                
                for k = 1:size(h1,1)
                    h2 = h1(k,1:5);
                    h4 = h1(k,end-5:end);
                    h6 = h1(k,5:end-5);
                    h3(k) = mean(h2(~isnan(h2)));
                    h5(k) = mean(h4(~isnan(h4)));
                    h7(k) = mean(h6(~isnan(h6)));
                end
                
                dh7 = diff(h7);
                dh7Mat = [dh7Mat; dh7];
                
                diffStart = [diffStart dh7(round(length(dh7)*.25))];
                diffMid = [diffMid dh7(round(length(dh7)*.6))];
                diffEnd = [diffEnd dh7(round(length(dh7)*.8))];
                diffdiff1 = [diffdiff1 dh7(round(length(dh7)/2)) - dh7(round(length(dh7)*.25))];
                diffdiff2 = [diffdiff2 dh7(round(length(dh7)*.8)) - dh7(round(length(dh7)/2))];
                rProfileOut = [rProfileOut h7'];
                
                micronsPerSlice = 1;
                XYpixelsPerSlice = 4.5*micronsPerSlice;
                
                pnDist1 = sqrt((...
                    eval(['Group' num2str(groupNum) '.E' num2str(embryoNum) '.PN1.xc;']) - ...
                    eval(['Group' num2str(groupNum) '.E' num2str(embryoNum) '.cellBody.xc;']))^2 + ...
                    (eval(['Group' num2str(groupNum) '.E' num2str(embryoNum) '.PN1.yc;']) - ...
                    eval(['Group' num2str(groupNum) '.E' num2str(embryoNum) '.cellBody.yc;']))^2 + ...
                    (XYpixelsPerSlice*(eval(['Group' num2str(groupNum) '.E' num2str(embryoNum) '.PN1.zc;']) - ...
                    max(eval(['Group' num2str(groupNum) '.E' num2str(embryoNum) '.cellBody.z;']))/2))^2);
                
                pnDist2 = sqrt((...
                    eval(['Group' num2str(groupNum) '.E' num2str(embryoNum) '.PN2.xc;']) - ...
                    eval(['Group' num2str(groupNum) '.E' num2str(embryoNum) '.cellBody.xc;']))^2 + ...
                    (eval(['Group' num2str(groupNum) '.E' num2str(embryoNum) '.PN2.yc;']) - ...
                    eval(['Group' num2str(groupNum) '.E' num2str(embryoNum) '.cellBody.yc;']))^2 + ...
                    (XYpixelsPerSlice*(eval(['Group' num2str(groupNum) '.E' num2str(embryoNum) '.PN2.zc;']) - ...
                    max(eval(['Group' num2str(groupNum) '.E' num2str(embryoNum) '.cellBody.z;']))/2))^2);
                
                pnSep = sqrt((...
                    eval(['Group' num2str(groupNum) '.E' num2str(embryoNum) '.PN1.xc;']) - ...
                    eval(['Group' num2str(groupNum) '.E' num2str(embryoNum) '.PN2.xc;']))^2 + ...
                    (eval(['Group' num2str(groupNum) '.E' num2str(embryoNum) '.PN1.yc;']) - ...
                    eval(['Group' num2str(groupNum) '.E' num2str(embryoNum) '.PN2.yc;']))^2 + ...
                    (XYpixelsPerSlice*(eval(['Group' num2str(groupNum) '.E' num2str(embryoNum) '.PN1.zc;']) - ...
                    eval(['Group' num2str(groupNum) '.E' num2str(embryoNum) '.PN2.zc;'])))^2);
                
                
                pnDistFromCenter1 = [pnDistFromCenter1 pnDist1];
                pnDistFromCenter2 = [pnDistFromCenter2 pnDist2];
                pnSeparation = [pnSeparation pnSep];
                
                
                figure(1);
                subplot(2,3,2);
                hold on;
                set(gca, 'fontsize', 14);
                plot([1 10], [0 find(h7 == max(h7),1) - find(h7 == min(h7),1) + rand(1)], 'Color', currColor, 'LineWidth', 2);
                xlabel('Normalized distance from cell center');
                ylabel('Intensity');
                title('Intensity Profile Along Cell Radius');
                %             xlim([1 10]);
                
                figure(1);
                subplot(2,3,3);
                hold on;
                set(gca, 'fontsize', 14);
                plot(dh7, 'Color', currColor, 'LineWidth', 2);
                xlabel('Normalized distance from cell center');
                ylabel('Intensity');
                title('Intensity Profile Along Cell Radius');
                %             xlim([1 10]);
                
                figure(1);
                subplot(2,3,4);
                hold on;
                set(gca, 'fontsize',14);
                plot(h7, 'Color', ...
                    currColor, 'LineWidth', 2);
                xlabel('Slice');
                ylabel('Intensity Profile Along All Slices');
                title('Average Radial Intensity Profile');
                %             xlim([1 10]);
            
            else
                
                vitGroup{i}(j) = NaN;
                blastForm{i}(j) = NaN;
                cellClump{i}(j) = NaN;
                
                % add NaNs equivalent to number of embryos in experiment
                stdevList = [stdevList NaN];
                vitGroupOut = [vitGroupOut NaN];
                cellClumpOut = [cellClumpOut NaN];
                blastFormOut = [blastFormOut NaN];
                diffStart = [diffStart NaN];
                diffMid = [diffMid NaN];
                diffEnd = [diffEnd NaN];
                diffdiff1 = [diffdiff1 NaN];
                diffdiff2 = [diffdiff2 NaN];
                darkPercent = [darkPercent NaN];
                rProfileOut = [rProfileOut NaN*ones(20,1)];
                pnDistFromCenter1 = [pnDistFromCenter1 NaN];
                pnDistFromCenter2 = [pnDistFromCenter2 NaN];
                pnSeparation = [pnSeparation NaN];
                
            end
            
            
        end
        
    else
        
        % add NaNs equivalent to number of embryos in experiment
        stdevList = [stdevList NaN*ones(1,length(vitGroup{i}))];
        vitGroupOut = [vitGroupOut NaN*ones(1,length(vitGroup{i}))];
        cellClumpOut = [cellClumpOut NaN*ones(1,length(vitGroup{i}))];
        blastFormOut = [blastFormOut NaN*ones(1,length(vitGroup{i}))];
        diffStart = [diffStart NaN*ones(1,length(vitGroup{i}))];
        diffMid = [diffMid NaN*ones(1,length(vitGroup{i}))];
        diffEnd = [diffEnd NaN*ones(1,length(vitGroup{i}))];
        diffdiff1 = [diffdiff1 NaN*ones(1,length(vitGroup{i}))];
        diffdiff2 = [diffdiff2 NaN*ones(1,length(vitGroup{i}))];
        darkPercent = [darkPercent NaN*ones(1,length(vitGroup{i}))];
        rProfileOut = [rProfileOut NaN*ones(20,length(vitGroup{i}))];
        pnDistFromCenter1 = [pnDistFromCenter1 NaN*ones(1,length(vitGroup{i}))];
        pnDistFromCenter2 = [pnDistFromCenter2 NaN*ones(1,length(vitGroup{i}))];
        pnSeparation = [pnSeparation NaN*ones(1,length(vitGroup{i}))];

    end
end


% figure, plot(linspace(0,100,19),mean(dh7Mat(blastFormOut(end-14:end) == 1 & vitGroupOut(end-14:end) == 1,:),1), 'Color', [1 0 0], 'LineWidth', 2);
% hold on, plot(linspace(0,100,19), mean(dh7Mat(blastFormOut(end-14:end) == 0 & vitGroupOut(end-14:end) == 1,:),1), 'Color', [0 0 .6], 'LineWidth', 2);
% set(gca, 'fontsize', 14);
% ylim([-.03 .02]);
% legend('Blast', 'No Blast', 'location', 'southwest');
% title('Radial Intensity Profile');

paramsOut = struct();
groundTruthOut = struct();

paramsOut.stdevList = stdevList;
paramsOut.diffStart = diffStart;
paramsOut.diffMid = diffMid;
paramsOut.diffEnd = diffEnd;
paramsOut.diffdiff1 = diffdiff1;
paramsOut.diffdiff2 = diffdiff2;
paramsOut.darkPercent = darkPercent;
paramsOut.rProfileOut = rProfileOut;
paramsOut.pnDistFromCenter1 = pnDistFromCenter1;
paramsOut.pnDistFromCenter2 = pnDistFromCenter2;
paramsOut.pnSeparation = pnSeparation;

groundTruthOut.vitGroup = vitGroupOut;
groundTruthOut.cellClump = cellClumpOut;
groundTruthOut.blastForm = blastFormOut;




        
        
        
        
        
        
        
        
        
        
