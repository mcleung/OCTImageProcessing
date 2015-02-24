% Plot and compare data from all embryos


%% 1. Re-calc params for all embryos

clear all;
close all;

load('embryoParams3.mat');

%% Plot

figure(99);
figure(100);

colorList = [0 0 1; ...
             0 .6 0; ...
             0 .6 .6; ...
             .8 0 0; ...
             .6 .6 0];

for i = 1:5
    
    figure(99);
    hold on;
    plot(eval(['Group1.E' num2str(i) '.stdev;']), 'Color', colorList(i,:), ...
        'LineWidth', 2);
    
    figure(100);
    hold on;
    plot(eval(['Group1.E' num2str(i) '.entropy;']), 'Color', colorList(i,:), ...
        'LineWidth', 2);
    
end

figure(99);
legend('E1', 'E2', 'E3', 'E4', 'E5');

figure(100);
legend('E1', 'E2', 'E3', 'E4', 'E5');

%% 

cellMask = Group1.E2.cellMask;
cell3D = Group1.E2.cell3D;
rList = Group1.E2.cellBody.r;
minPNslice = min(Group1.E2.sliceList);
maxPNslice = max(Group1.E2.sliceList);
xCenter = Group1.E2.cellBody.xc;
yCenter = Group1.E2.cellBody.yc;

A = imadjust(cell3D(:,:,24) .* cellMask(:,:,24), [.25 .6], [0 1]);
A1 = medfilt2(A, [15 15]);
A2 = medfilt2(A1, [15 15]);

figure(10);
subplot(1,3,1); imshow(A);
subplot(1,3,2); imshow(A1);
subplot(1,3,3); imshow(A2);

std(A(A > 0))
std(A1(A1 > 0))
std(A2(A2 > 0))

%% Re-extract params from all embryos in groups 1 and 2

clear all;
% close all;
% 
embryoGroup = [1 1 1 1 1 2 2 2 2 2];
blastForm = [0 0 1 1 1 0 0 0 1 1];
cellClump = [1 0 1 0 0 1 1 1 1 0];
vitGroup = [0 0 0 0 0 1 1 1 1 1];

% embryoGroup = [1 1 1 1 1 2 2 2 2 2 ...
%                1 1 5 5 5 2 2 2 2 2 ...
%                1 1 1 1 2 2 2 2 2 2];
% blastForm = [0 0 1 1 1 0 0 0 1 1 ...
%              NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN ...
%              NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN];
% cellClump = [1 0 1 0 0 1 1 1 1 0 ...
%              0 0 NaN NaN NaN 1 0 1 1 1 ...
%              0 0 0 0 1 1 1 1 1 1]; % 1 = clump, 0 = no clump, 2 = can't tell

posClass = 0;
negClass = 1-posClass;
groundTruth = vitGroup;
[stdevList entropyVals IntensityRatio1 IntensityRatio2] = ...
    loadAllParams(embryoGroup, groundTruth, posClass);


% plot
% close all;
fig = figure;
scatter3(entropyVals(groundTruth == posClass), IntensityRatio2(groundTruth == posClass), ...
    stdevList(groundTruth == posClass),  150, [0 0 .6], 'filled');
hold on;
scatter3(entropyVals(groundTruth == negClass), IntensityRatio2(groundTruth == negClass), ...
    stdevList(groundTruth == negClass), 150, [1 0 0], 'filled');
% hold on;
% scatter3(iVal(cellClump == 2), iRatio(cellClump == 2), ...
%     stdevList(cellClump == 2), 150, [1 0 0], 'filled');

set(gca, 'fontsize', 14);
xlabel('Entropy');
ylabel('Intensity Ratio');
zlabel('Standard Deviation');
% title('Predicting Blastocyst Formation');
title('Predicting Vitrification');
% title('Parameters Predictive of Cytoplasmic Clumping');
% legend('No Clumping', 'Clumping', 'Location', 'NorthEast');
legend('Not Vitrified', 'Vitrified', 'Location', 'North');
% legend('Blastocyst', 'No Blastocyst', 'Location', 'NorthEast');

% axis([1 1.8 1 1.8 .04 .075]);
% axis([1 2.2 1 2.2 .04 .095]);
axis([3.8 5.4 1 2 .04 .095]);
% axis([4 4.8 1 1.7 .04 .07]);

% axis([min(IntensityRatio1(groundTruth < 2 & ~isnan(groundTruth))) ...
%     max(IntensityRatio1(groundTruth < 2 & ~isnan(groundTruth))) ...
%     min(IntensityRatio2(groundTruth < 2 & ~isnan(groundTruth))) ...
%     max(IntensityRatio2(groundTruth < 2 & ~isnan(groundTruth))) ...
%     min(stdevList(groundTruth < 2 & ~isnan(groundTruth))) ...
%     max(stdevList(groundTruth < 2 & ~isnan(groundTruth)))]);
set(gca, 'view', [-113 50]);


%% Classification

cpEmbryos = classperf(groundTruth(~isnan(groundTruth))', ...
    'positive', posClass, 'negative', negClass);

featureVec = [entropyVals' IntensityRatio2' IntensityRatio1', ...
    stdevList'];

embryoClassifier = svmtrain(featureVec(:,[1 2 3 4]), groundTruth');
[predictOut, decDist] = svmclassify(embryoClassifier, featureVec(:,[1 2 3 4]));

classperf(cpEmbryos, predictOut(~isnan(groundTruth)))
    
true_pos = length(groundTruth(groundTruth' == posClass & predictOut == posClass));
false_pos = length(groundTruth(groundTruth' == negClass & predictOut == posClass));
true_neg = length(groundTruth(groundTruth' == negClass & predictOut == negClass));
false_neg = length(groundTruth(groundTruth' == posClass & predictOut == negClass));

sensitivity = true_pos / (true_pos + false_neg)
specificity = true_neg / (true_neg + false_pos)
ppv = true_pos / length(predictOut(predictOut == posClass))
accuracy

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


















