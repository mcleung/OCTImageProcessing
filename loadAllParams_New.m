
% load params from all dates

function [paramsOut, groundTruthOut] = loadAllParams_New(exptsToLoad, varToPlot, posClass)

if nargin < 2
    varToPlot = 'blastForm';
    posClass = 1;
end

% initialize all output parameters
blastFormOut = [];
cellClumpOut = [];
vitGroupOut = [];

stdevList = [];
stdNorm = [];
entropyVals = [];
homogeneityAvg = [];
energyAvg = [];
homogeneityDiff = [];
energyDiff = [];
IntensityRatio1 = [];
IntensityRatio2 = [];
IntensityDiff = [];
IntensityAvg = [];
minMaxIndDiff = [];
diffStart = [];
diffMid = [];
diffEnd = [];
diffdiff1 = [];
diffdiff2 = [];
darkPercent = [];


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

blastForm{1} = NaN*ones(1,10);
blastForm{2} = NaN*ones(1,10);
blastForm{3} = [0 0 1 1 1 0 0 0 1 1];
blastForm{4} = floor([1 1 0 1 1 1 1 0 0 .6 1 0 0 .6 0]); % .6 means very early blast or unsure
blastForm{5} = floor([1 1 1 0 1 .6 1 1 1 0 .6 0 .6 1 1]); % .6 means very early blast or unsure
blastForm{6} = floor([1 0 0 0 1 1 1 1 0 1 0 1 1 0 .6 1 0]); % .6 means very early blast or unsure

% manual evaluation of clumping
cellClump{1} = [0 0 0 0 1 1 1 1 1 1]; % 1 = clump, 0 = no clump, 2 = can't tell
cellClump{2} = [0 0 0 0 1 1 0 1 1 1];
cellClump{3} = [1 0 1 0 0 1 1 1 1 1];
cellClump{4} = [1 0 0 0 1 1 1 1 1 1 0 1 1 1 0]; % 
cellClump{5} = [1 0 0 1 1 0 1 1 1 1 1 0 1 1 1]; % 
cellClump{6} = [1 0 1 0 0 0 1 1 1 0 1 0 1 1 1 1 1]; % 

fileNameAssoc{1} = 'embryoParams3.mat'; % 4-26-13 and 5-3-13
fileNameAssoc{2} = 'embryoParams2.mat'; % 12-4-13
fileNameAssoc{3} = 'embryoParams.mat'; % 1-14-14
fileNameAssoc{4} = 'embryoParams4.mat'; % 3-27-14
fileNameAssoc{5} = 'embryoParams5.mat'; % 5-29-14
fileNameAssoc{6} = 'embryoParams6.mat'; % 9-25-14

dh7Mat = [];

for i = 1:length(exptsToLoad)
    
    if (exptsToLoad(i) && exist(fileNameAssoc{i}, 'file'))
        
        load(fileNameAssoc{i});
        numNonVit = sum(1 - vitGroup{i});
        
        % add ground truth data to corresponding output arrays
        vitGroupOut = [vitGroupOut vitGroup{i}];
        cellClumpOut = [cellClumpOut cellClump{i}];
        blastFormOut = [blastFormOut blastForm{i}];
        
        for j = 1:length(vitGroup{i})
            
            j
            if j < (numNonVit + 1)
                groupNum = 1;
                embryoNum = j;
            else
                groupNum = 2;
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
            
            cellClump{i}(j)
            mean(cytoplasmPixels(cytoplasmPixels > 0))
            std(cytoplasmPixels(cytoplasmPixels > 0))
            sum(sum(cytoplasmPixels > 0 & cytoplasmPixels < (meanBright - .05)))/numXY
            figure(3);
            imshow(cytoplasmPixels);
            
            figure(4);
            imhist(cytoplasmPixels);
            
            entropyCurr = eval(['Group' num2str(groupNum) '.E' ...
                num2str(embryoNum) '.entropyfilt;']);
            stdevCurr = eval(['Group' num2str(groupNum) '.E' ...
                num2str(embryoNum) '.stdev;']);
            homogeneityCurr = mean(eval(['Group' num2str(groupNum) '.E' ...
                num2str(embryoNum) '.homogeneity;']),1);
            energyCurr = mean(eval(['Group' num2str(groupNum) '.E' ...
                num2str(embryoNum) '.energy;']),1);
            
            entropyVals = [entropyVals mean(entropyCurr)];
            stdevList = [stdevList mean(stdevCurr)];
            homogeneityAvg = [homogeneityAvg mean(mean(homogeneityCurr))];
            energyAvg = [energyAvg mean(mean(energyCurr))];
            
            hAvg = mean(homogeneityCurr,1);
            eAvg = mean(energyCurr,1);
            homogeneityDiff = [homogeneityDiff hAvg(end)-hAvg(1)];
            energyDiff = [energyDiff eAvg(end)-eAvg(1)];
            

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
            
            IntensityRatio1 = [IntensityRatio1 (max(h5)/min(h5))];
            IntensityRatio2 = [IntensityRatio2 (max(h3)/min(h3))];
            IntensityDiff = [IntensityDiff (max(h7) - min(h7))];
            IntensityAvg = [IntensityAvg; h7];
            minMaxIndDiff = [minMaxIndDiff find(h7 == max(h7),1) - find(h7 == min(h7),1)];
            diffStart = [diffStart dh7(round(length(dh7)*.25))];
            diffMid = [diffMid dh7(round(length(dh7)/2))];
            diffEnd = [diffEnd dh7(round(length(dh7)*.8))];
            diffdiff1 = [diffdiff1 dh7(round(length(dh7)/2)) - dh7(round(length(dh7)*.25))];
            diffdiff2 = [diffdiff2 dh7(round(length(dh7)*.8)) - dh7(round(length(dh7)/2))];
            stdNorm = [stdNorm mean(stdevCurr)/mean(h7)];
            
            
            
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
            
            figure(1); 
            subplot(2,3,5);
            hold on;
            set(gca, 'fontsize',14);
            plot(linspace(0,1,length(mean(homogeneityCurr,1))), ...
                mean(homogeneityCurr,1), 'Color', currColor, 'LineWidth', 2);
            xlabel('Slice');
            ylabel('homogeneity across slices');
            title('homogeneity');
            xlim([0 1]);
            
            figure(1); 
            subplot(2,3,6);
            hold on;
            set(gca, 'fontsize',14);
            plot(linspace(0,1,length(mean(energyCurr,1))), ...
                mean(energyCurr,1), 'Color', currColor, 'LineWidth', 2);
            xlabel('Slice');
            ylabel('Energy across slices');
            title('Energy');
            xlim([0 1]);
            
            1
            
        end
        
    else
        
        % add NaNs equivalent to number of embryos in experiment
        entropyVals = [entropyVals NaN*ones(1,length(vitGroup{i}))];
        stdevList = [stdevList NaN*ones(1,length(vitGroup{i}))];
        stdNorm = [stdNorm NaN*ones(1,length(vitGroup{i}))];
        homogeneityAvg = [homogeneityAvg NaN*ones(1,length(vitGroup{i}))];
        energyAvg = [energyAvg NaN*ones(1,length(vitGroup{i}))];
        homogeneityDiff = [homogeneityDiff NaN*ones(1,length(vitGroup{i}))];
        energyDiff = [energyDiff NaN*ones(1,length(vitGroup{i}))];
        IntensityRatio1 = [IntensityRatio1 NaN*ones(1,length(vitGroup{i}))];
        IntensityDiff = [IntensityDiff NaN*ones(1,length(vitGroup{i}))];
        IntensityRatio2 = [IntensityRatio2 NaN*ones(1,length(vitGroup{i}))];
        IntensityAvg = [IntensityAvg; NaN*ones(length(vitGroup{i}),20)];
        vitGroupOut = [vitGroupOut NaN*ones(1,length(vitGroup{i}))];
        cellClumpOut = [cellClumpOut NaN*ones(1,length(vitGroup{i}))];
        blastFormOut = [blastFormOut NaN*ones(1,length(vitGroup{i}))];
        minMaxIndDiff = [minMaxIndDiff NaN*ones(1,length(vitGroup{i}))];
        diffStart = [diffStart NaN*ones(1,length(vitGroup{i}))];
        diffMid = [diffMid NaN*ones(1,length(vitGroup{i}))];
        diffEnd = [diffEnd NaN*ones(1,length(vitGroup{i}))];
        diffdiff1 = [diffdiff1 NaN*ones(1,length(vitGroup{i}))];
        diffdiff2 = [diffdiff2 NaN*ones(1,length(vitGroup{i}))];
        darkPercent = [darkPercent NaN*ones(1,length(vitGroup{i}))];

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

paramsOut.entropyVals = entropyVals;
paramsOut.stdevList = stdevList;
paramsOut.stdNorm = stdNorm;
paramsOut.homogeneityAvg = homogeneityAvg;
paramsOut.energyAvg = energyAvg;
paramsOut.homogeneityDiff = homogeneityDiff;
paramsOut.energyDiff = energyDiff;
paramsOut.IntensityRatio1 = IntensityRatio1;
paramsOut.IntensityDiff = IntensityDiff;
paramsOut.IntensityRatio2 = IntensityRatio2;
paramsOut.IntensityAvg = IntensityAvg;
paramsOut.minMaxIndDiff = minMaxIndDiff;
paramsOut.diffStart = diffStart;
paramsOut.diffMid = diffMid;
paramsOut.diffEnd = diffEnd;
paramsOut.diffdiff1 = diffdiff1;
paramsOut.diffdiff2 = diffdiff2;
paramsOut.darkPercent = darkPercent;

groundTruthOut.vitGroup = vitGroupOut;
groundTruthOut.cellClump = cellClumpOut;
groundTruthOut.blastForm = blastFormOut;




        
        
        
        
        
        
        
        
        
        
