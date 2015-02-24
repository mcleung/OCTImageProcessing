% Plot and compare data from all embryos


%% Re-extract params from all embryos in groups 1 and 2

% clear all;
% close all;

varToPlot = 'blastForm'; % blastForm, cellClump, or vitGroup
posClass = 1;
negClass = 1-posClass;
exptsToLoad = [0 0 1 1 1 1];
embryosInExpt = [10 10 10 15 15 17 10 28];
[paramsOut, groundTruthOut] = loadAllParams_New2(exptsToLoad, varToPlot, posClass);

vitGroupAll = groundTruthOut.vitGroup;
blastFormAll = groundTruthOut.blastForm;
cellClumpAll = groundTruthOut.cellClump;
rProfileOut = paramsOut.rProfileOut;

%% Try new params

figure(1);
clf;
groundTruth = eval(['groundTruthOut.' varToPlot]);
allProfile = [];

testParam1 = NaN*zeros(1,length(groundTruth));
testParam2 = NaN*zeros(1,length(groundTruth));
testParam3 = NaN*zeros(1,length(groundTruth));

hProfile = zeros(1,20);
cProfile = zeros(1,20);

for i = 1:length(groundTruth)
    
    if ~isnan(groundTruth(i)) 
    
        currProfile = rProfileOut(:,i);
        currDiff = smooth(diff(currProfile));
        allProfile = [allProfile currProfile];
        rAvg = mean(currProfile);
        
        posPeaks = round(mspeaks(1:20, currProfile'));
        negPeaks = (mspeaks(1:20, 1-currProfile'));
        
        if length(negPeaks) < 1% || length(negPeaks) > 2
            negLocs = 16;
        elseif negPeaks(end,1) < 6
            negLocs = 10;
        else
            negLocs = round(negPeaks(end,1));
        end       
        
        if size(posPeaks,1) < 2 || size(posPeaks,1) > 3
            if size(posPeaks,1) == 1
                if posPeaks(1) < 10
                    posLocs = [posPeaks(1) min(negLocs+2,19)];
                else
                    posLocs = [6 posPeaks(1)];
                end
            else
                posLocs = [8 18];
            end
        else
            posLocs = round(posPeaks(end-1:end,1))';
        end
        
        if ~isempty(intersect(i, [56 72 76]))  
            posLocs = [12 19];
            negLogs = 16;
        end       
        
        testParam1(i) = currProfile(posLocs(1)) - currProfile(negLocs);
        testParam2(i) = testParam1(i)/(rAvg*(max(posLocs(1) - negLocs,1)));
%         testParam3(i) = currProfile(posLocs(2)) - currProfile(negLocs)/rAvg;
        testParam3(i) = currDiff(6);
        
        if groundTruth(i)
            currColor = [0 0 1];
            cProfile = [cProfile ; currProfile'];
        else
            currColor = [1 0 0];
            hProfile = [hProfile ; currProfile'];
        end
        
%         hold on;
%         plot(5:5:100, currProfile, 'color', [.7 .7 0], 'linewidth', 2);
%         ylim([.2 .8]);
         
        hold on;
        plot(10:5:100, smooth(diff(currProfile)), 'color', currColor, 'linewidth', 2);
%         ylim([.2 .8]);       
        
    end
end

figure(2);
clf;
plot(5:5:100, mean(cProfile,1), 'color', [1 0 0], 'linewidth', 3);
hold on;
plot(5:5:100, mean(hProfile,1), 'color', [0 0 1], 'linewidth', 3);

jbfill(5:5:100, mean(cProfile,1)+std(cProfile,1), mean(cProfile,1)-std(cProfile,1), ...
    [1 0 0], 'none', [], .3);
jbfill(5:5:100, mean(hProfile,1)+std(hProfile,1), mean(hProfile,1)-std(hProfile,1), ...
    [0 0 1], 'none', [], .3);
xlim([5 100]);

set(gca, 'fontsize', 14);
xlabel('% of max radius');
ylabel('average intensity');
title('Radial Profiles');
grid on;

%% 
fList = fieldnames(paramsOut);
startingEmbryoList = [1 cumsum(embryosInExpt)+1];
        
% combine batches
for i = 1:length(exptsToLoad)   
    if exptsToLoad(i)
        
        currFirst = startingEmbryoList(i);
        currLast = startingEmbryoList(i+1) - 1;
        
        testParam1(currFirst:currLast) = ...
            testParam1(currFirst:currLast) - ...
            mean(testParam1(currFirst:currLast));
        testParam2(currFirst:currLast) = ...
            testParam2(currFirst:currLast) - ...
            mean(testParam2(currFirst:currLast));
        testParam3(currFirst:currLast) = ...
            testParam3(currFirst:currLast) - ...
            mean(testParam3(currFirst:currLast));

        for j = 1:length(fList)
            
            % mean of each experiment is shifted to 0
            currParamList = paramsOut.(fList{j});
            currParamList(currFirst:currLast) = ...
                currParamList(currFirst:currLast) - ...
                mean(currParamList(currFirst:currLast));
            
            paramsOut.(fList{j}) = currParamList;

        end
    end
end

%% 

% p1 = testParam1;
% p2 = testParam2;
% p3 = testParam3;

p1 = paramsOut.diffdiff1;
p2 = paramsOut.pnDistFromCenter1 + paramsOut.pnDistFromCenter1;
p3 = paramsOut.pnSeparation;

groundTruth = eval(['groundTruthOut.' varToPlot]);
vitGroup = groundTruthOut.vitGroup;

% plot
% close all;
fig = figure; % & vitGroup == 1
% scatter3(p1(groundTruth == posClass), p2(groundTruth == posClass), ...
%     p3(groundTruth == posClass),  150, [0 0 .6], 'filled');
% hold on;
% scatter3(p1(groundTruth == negClass), p2(groundTruth == negClass), ...
%     p3(groundTruth == negClass), 150, [1 0 0], 'filled');
scatter3(p1(groundTruth == posClass), p2(groundTruth == posClass), ...
    p3(groundTruth == posClass),  150, [0 0 .6], 'filled');
hold on;
scatter3(p1(groundTruth == negClass), p2(groundTruth == negClass), ...
    p3(groundTruth == negClass), 150, [1 0 0], 'filled');
view(0,90);

a = (1:sum(exptsToLoad .* embryosInExpt))';
b = num2str(a);
c = cellstr(b);
dx = -0.01; dy = 0.01; dz = .002; % displacement so the text does not overlay the data points
% text(p1(~isnan(p1))+dx, p2(~isnan(p1))+dy, p3(~isnan(p1))+dz, c);
set(gca, 'fontsize', 14);
xlabel('Param1');
ylabel('Param2');
zlabel('Param3');
% axis([-.06 .04 -.1 .05 -.025 .015])

%% find statistically significant parameters
% For separating clump from non-clump: entropyVals, diffdiff2,
% diffEnd, IntensityDiff, stdevList, energyAvg, 
% For separating blast from non-blast: stdevList, IntensityDiff,
% homogeneityAvg, diffdiff1, diffdiff2, diffMid

p4 = paramsOut.pnSeparation;
[h p] = ttest2(p4(groundTruth == posClass), p4(groundTruth == negClass))



%% Maker kernel density plots

paramName = 'pnDistFromCenter1';
outcomeName = 'cellClump';

paramToEval = eval(['paramsOut.' paramName]);
groundTruth = eval([outcomeName 'All']);
posClass = 1;

figure(1);
colorMat = zeros(length(paramToEval),3);
colorMat(groundTruth == 1-posClass,:) = repmat([0 0 1], sum(1-groundTruth(~isnan(groundTruth))), 1);
colorMat(groundTruth == posClass,:) = repmat([1 0 0], sum(groundTruth(~isnan(groundTruth))), 1);

% plot([-.3 .3], [0 0], '--k');
% hold on;
% h1 = scatter(-1*paramToEval, rand(1,length(paramToEval)) - .5, 100, colorMat, 'filled');
% set(h1, 'Marker', 'o');
set(gca, 'FontSize', 14);
% title('Embryo Clumping');
xlabel(paramName);
ylabel('histogram');

[f0, x0] = ksdensity(-1*paramToEval(groundTruth == 1-posClass));
[f1, x1] = ksdensity(-1*paramToEval(groundTruth == posClass));

hold on;
v0 = plot(x0, f0, 'Color', [0 0 1], 'linewidth', 2);
plot(x0, -f0, 'Color', [0 0 1], 'linewidth', 2);
v1 = plot(x1, f1, 'Color', [1 0 0], 'linewidth', 2);
plot(x1, -f1, 'Color', [1 0 0], 'linewidth', 2);
h = area(x0, f0, 'EdgeColor', [0 0 1], 'FaceColor', [0 0 1]);
hc = get(h, 'Children');
set(hc, 'FaceAlpha', .25);
h = area(x0, -f0, 'EdgeColor', [0 0 1], 'FaceColor', [0 0 1]);
hc = get(h, 'Children');
set(hc, 'FaceAlpha', .25);
h = area(x1, f1, 'EdgeColor', [1 0 0], 'FaceColor', [1 0 0]);
hc = get(h, 'Children');
set(hc, 'FaceAlpha', .25);
h = area(x1, -f1, 'EdgeColor', [1 0 0], 'FaceColor', [1 0 0]);
hc = get(h, 'Children');
set(hc, 'FaceAlpha', .25);
% line([0 .25], [0 0], 'Color', 'k', 'linewidth', 2);
% xlim([-.08 .08]);
legend([v0, v1], ['No ' outcomeName], outcomeName);
grid on;

% [xR, yR] = perfcurve(groundTruth, paramToEval, 1-posClass);
% figure, plot(xR, yR, 'linewidth', 2);
% set(gca, 'fontsize', 14);
% grid on;
% xlabel('1 - specificity');
% ylabel('sensitivity');
% title('Automatic Clumping Classification');
% 
% xRD = diff(xR);
% AUC = sum(xRD.*yR(1:end-1))


%% Classification

groundTruth = cellClumpAll;
groundTruthN = groundTruth(~isnan(groundTruth));
posClass = 0;
negClass = 1-posClass;
cpEmbryos = classperf(groundTruthN', 'positive', posClass, 'negative', negClass);
% featureVec = [allProfile; testParam1(~isnan(groundTruth)); ...
%     testParam2(~isnan(groundTruth)); testParam3(~isnan(groundTruth))];
featureVec = [testParam1(~isnan(groundTruth)); ...
   testParam2(~isnan(groundTruth)); testParam3(~isnan(groundTruth))];

AUC_list = zeros(1,1);
   
for j = 1:1   
   
nGroups = 10;
C = cvpartition(length(groundTruthN), 'Kfold', nGroups);

predictOutAll = zeros(1,length(groundTruthN));
decDistAll = zeros(1,length(groundTruthN));

for i = 1:nGroups
    
    % first, data is partitioned into training and testing set
    test = C.test(i);
    train = C.training(i);
    
    VecTrain = featureVec(:, train);
    VecTest = featureVec(:, test);
    groupsTrain = groundTruthN(train);
    
    embryoClassifier = svmtrain(VecTrain', groupsTrain', ...
        'kernel_function', 'rbf', 'rbf_sigma', 2, 'boxconstraint', 1);
    [predictOut, decDist] = svmclassify(embryoClassifier, VecTest');
    predictOutAll(test) = predictOut;
    decDistAll(test) = decDist;

end

[xR, yR] = perfcurve(groundTruthN, decDistAll, posClass);
figure(1);
clf;
plot(xR, yR, 'linewidth', 2);
set(gca, 'fontsize', 14);
grid on;
xlabel('1 - specificity');
ylabel('sensitivity');
title('Automatic Clumping Classification');

xRD = diff(xR);
AUC_list(j) = sum(xRD.*yR(1:end-1));

end

mean(AUC_list)


% embryoClassifier = svmtrain(featureVec', groundTruthN', ...
%     'kernel_function', 'rbf', 'rbf_sigma', 2, 'boxconstraint', 1, ...
%     'showplot', 'true');
% 


%%
% plot
% close all;
fig = figure;
scatter3(entropyVals(predictOut == posClass), IntensityRatio2(predictOut == posClass), ...
    stdevList(predictOut == posClass),  150, [0 .6 .6], 'filled');
hold on;
scatter3(entropyVals(predictOut == negClass), IntensityRatio2(predictOut == negClass), ...
    stdevList(predictOut == negClass), 150, [.7 .7 0], 'filled');
% hold on;
% scatter3(iVal(cellClump == 2), iRatio(cellClump == 2), ...
%     stdevList(cellClump == 2), 150, [1 0 0], 'filled');

set(gca, 'fontsize', 14);
xlabel('Entropy');
ylabel('Intensity Ratio');
zlabel('Standard Deviation');
% title('Predicting Blastocyst Formation');
% title('Predicting Vitrification');
title('Predicted Class');
legend('No Clumping', 'Clumping', 'Location', 'North');
% legend('Not Vitrified', 'Vitrified', 'Location', 'North');
% legend('Blastocyst', 'No Blastocyst', 'Location', 'North');

% axis([1 1.8 1 1.8 .04 .075]);
% axis([1 2.2 1 2.2 .04 .095]);
axis([3.8 5.4 1 2 .04 .095]);

% axis([min(IntensityRatio1(groundTruth < 2 & ~isnan(groundTruth))) ...
%     max(IntensityRatio1(groundTruth < 2 & ~isnan(groundTruth))) ...
%     min(IntensityRatio2(groundTruth < 2 & ~isnan(groundTruth))) ...
%     max(IntensityRatio2(groundTruth < 2 & ~isnan(groundTruth))) ...
%     min(stdevList(groundTruth < 2 & ~isnan(groundTruth))) ...
%     max(stdevList(groundTruth < 2 & ~isnan(groundTruth)))]);
set(gca, 'view', [-113 50]);









%%


for i = 1:5

i
cellMask = eval(['Group1.E' num2str(i) '.cellMask;']);
cell3D = eval(['Group1.E' num2str(i) '.cell3D;']);
minPNslice = min(eval(['Group1.E' num2str(i) '.sliceList;']));
maxPNslice = max(eval(['Group1.E' num2str(i) '.sliceList;']));
rList = eval(['Group1.E' num2str(i) '.cellBody.r']);
xCenter = eval(['Group1.E' num2str(i) '.cellBody.xc']);
yCenter = eval(['Group1.E' num2str(i) '.cellBody.yc']);

%     paramList = extractCytoplasmParams(cell3D, cellMask, minPNslice, maxPNslice);
rProfile = extractRadialProfile(cell3D, cellMask, minPNslice, maxPNslice, rList, xCenter, yCenter);
eval(['Group1.E' num2str(i) '.rProfile = rProfile;']);

entropyVals(i) = mean(eval(['Group1.E' num2str(i) '.entropyfilt;']));
stdevList(i) = mean(eval(['Group1.E' num2str(i) '.stdev;']));
entropyList(i) = mean(eval(['Group1.E' num2str(i) '.entropy;']));

%     h1 = paramList.homogeneity(4,:) ./ paramList.homogeneity(7,:);
%     h2 = paramList.homogeneity(15,:) ./ paramList.homogeneity(18,:);
%
%     hRatio1(i) = mean(h1(~isnan(h1)));
%     hRatio2(i) = mean(h2(~isnan(h2)));
%
%     entropyVals(i) = mean(paramList.entropyfilt);
%     stdevList(i) = mean(paramList.stdev);
%     entropyList(i) = mean(paramList.entropy);
%
%     eval(['Group1.E' num2str(i) '.stdev = paramList.stdev;']);
%     eval(['Group1.E' num2str(i) '.entropy = paramList.entropy;']);
%     eval(['Group1.E' num2str(i) '.entropyfilt = paramList.entropyfilt;']);
%     eval(['Group1.E' num2str(i) '.homogeneity = paramList.homogeneity;']);
%     eval(['Group1.E' num2str(i) '.energy = paramList.energy;']);

%     entropyList(i) = entropy(cellMask(:,:,minPNslice+1:maxPNslice-1) .* ...
%         cell3D(:,:,minPNslice+1:maxPNslice-1));
end

for i = 1:5
    
    i
    cellMask = eval(['Group2.E' num2str(i) '.cellMask;']);
    cell3D = eval(['Group2.E' num2str(i) '.cell3D;']);
    minPNslice = min(eval(['Group2.E' num2str(i) '.sliceList;']));
    maxPNslice = max(eval(['Group2.E' num2str(i) '.sliceList;']));
    rList = eval(['Group2.E' num2str(i) '.cellBody.r']);
    xCenter = eval(['Group2.E' num2str(i) '.cellBody.xc']);
    yCenter = eval(['Group2.E' num2str(i) '.cellBody.yc']);
    
    rProfile = extractRadialProfile(cell3D, cellMask, minPNslice, maxPNslice, rList, xCenter, yCenter);    
    eval(['Group2.E' num2str(i) '.rProfile = rProfile;']);
    
    entropyVals(i+5) = mean(eval(['Group2.E' num2str(i) '.entropyfilt;']));
    stdevList(i+5) = mean(eval(['Group2.E' num2str(i) '.stdev;']));
    entropyList(i+5) = mean(eval(['Group2.E' num2str(i) '.entropy;']));    
    
%     paramList = extractCytoplasmParams(cell3D, cellMask, minPNslice, maxPNslice);
%     h1 = paramList.homogeneity(4,:) ./ paramList.homogeneity(7,:);
%     h2 = paramList.homogeneity(15,:) ./ paramList.homogeneity(18,:);
%     
%     hRatio1(i+5) = mean(h1(~isnan(h1)));
%     hRatio2(i+5) = mean(h2(~isnan(h2)));  
%     
%     entropyVals(i+5) = mean(paramList.entropyfilt);
%     entropyList(i+5) = mean(paramList.entropy);
%     stdevList(i+5) = mean(paramList.stdev);
%     
%     eval(['Group2.E' num2str(i) '.stdev = paramList.stdev;']);
%     eval(['Group2.E' num2str(i) '.entropy = paramList.entropy;']);
%     eval(['Group2.E' num2str(i) '.entropyfilt = paramList.entropyfilt;']);
%     eval(['Group2.E' num2str(i) '.homogeneity = paramList.homogeneity;']);
%     eval(['Group2.E' num2str(i) '.energy = paramList.energy;']);

    
    %     entropyList(i+5) = entropy(cellMask(:,:,minPNslice+1:maxPNslice-1) .* ...
%         cell3D(:,:,minPNslice+1:maxPNslice-1));
  
end

save('embryoParams.mat', 'Group1', 'Group2');

% figure, scatter3(hRatio1(~cellClump), hRatio2(~cellClump), ...
%     stdevList(~cellClump),  150, [1 0 0], 'filled');
% hold on;
% scatter3(hRatio1(cellClump), hRatio2(cellClump), ...
%     stdevList(cellClump), 150, [0 0 1], 'filled');


%%

close all;
figure(1);
figure(2);

for i = 1:10
    
    if blastForm(i) == 1
        currColor = [0 0 .6];
    elseif blastForm(i) == 0
        currColor = [1 0 0];
    else
        currColor = [1 0 0];
    end
    
    if i < 6
        groupNum = 1;
        embryoNum = i;
    else
        groupNum = 2;
        embryoNum = i - 5;
    end
        
    h1 = eval(['Group' num2str(groupNum) '.E' num2str(embryoNum) ...
        '.rProfile']);
    
    for j = 1:size(h1,1)
        h2 = h1(j,1:5);
        h4 = h1(j,end-5:end);
        h3(j) = mean(h2(~isnan(h2)));
        h5(j) = mean(h4(~isnan(h4)));
    end
        

    iRatio(i) = max(h5)/min(h5);
    iVal(i) = max(h3)/min(h3); % mean(h3(end-3:end));
%     eRatio(i) = h1(13)/h1(20);
%     hRatio1(i) = h1(19)/h1(20);
%     hRatio2(i) = h1(1)/h1(2);
    
    figure(1);
    hold on;
    plot(h3, 'Color', currColor);
    
    figure(2);
    hold on;
    plot(h5, 'Color', currColor);
    
end

% ylim([.45 .85]);

figure(3);
scatter3(iVal(blastForm == 0), iRatio(blastForm == 0), ...
    stdevList(blastForm == 0),  150, [1 0 0], 'filled');
hold on;
scatter3(iVal(blastForm == 1), iRatio(blastForm == 1), ...
    stdevList(blastForm == 1), 150, [0 0 .6], 'filled');
% hold on;
% scatter3(iVal(cellClump == 2), iRatio(cellClump == 2), ...
%     stdevList(cellClump == 2), 150, [1 0 0], 'filled');

set(gca, 'fontsize', 14);
xlabel('Intensity Profile 1');
ylabel('Intensity Profile 2');
zlabel('Intensity Stdev');
title('Predicting Blastocyst Formation');
legend('Blastocyst', 'No Blastocyst', 'Location', 'North');
xlim([1 1.6]);


















