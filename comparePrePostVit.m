% plot embryos to compare pre and post-vit
% for fresh embryos from 1-23-15


%% Re-extract params from all embryos in groups 1 and 2

% clear all;
% close all;
addpath('C:/Users/Livia/Desktop/IVF/Code/EmbryoProject/functions/')

varToPlot = 'vitGroup'; % blastForm, cellClump, or vitGroup
posClass = 1;
negClass = 1-posClass;
exptsToLoad = [0 0 0 0 0 0 0 1];
embryosInExpt = [10 10 10 15 15 17 10 28];
[paramsOut, groundTruthOut] = loadAllParams_New2(exptsToLoad, varToPlot, posClass);

vitGroupAll = groundTruthOut.vitGroup;
blastFormAll = groundTruthOut.blastForm;
cellClumpAll = groundTruthOut.cellClump;
rProfileOut = paramsOut.rProfileOut;

%% Compare slope param

figure(1);
clf;
groundTruth = vitGroupAll;
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
plot(5:5:100, mean(hProfile,1), 'color', [0 0 1], 'linewidth', 3);
hold on;
plot(5:5:100, mean(cProfile,1), 'color', [1 0 0], 'linewidth', 3);
hold on;

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
legend('Fresh', 'Vitrified')




%% Set parameters and plot

exptIndices = (repmat(cumsum(embryosInExpt(1:end-1))', 1, 2))';
exptIndices = exptIndices(:)';
offsetIndices = [zeros(length(exptsToLoad)-1,1) ones(length(exptsToLoad)-1,1)]'; 
offsetIndices = offsetIndices(:)';
exptIndices = [1, exptIndices + offsetIndices, sum(embryosInExpt)];

indicesToLoad = zeros(1,sum(embryosInExpt));

for i = 1:length(exptsToLoad)   
    if (exptsToLoad(i))
        indicesToLoad(exptIndices(2*i-1):exptIndices(2*i)) = 1;
    end
end

vitGroupAll = vitGroupAll(indicesToLoad == 1)
blastFormAll = blastFormAll(indicesToLoad == 1)
cellClumpAll = cellClumpAll(indicesToLoad == 1)
rProfileOut = rProfileOut(:,(indicesToLoad == 1))

%% 

p1 = testParam1;
% p1 = paramsOut.stdevList;
p1 = p1(indicesToLoad == 1)

figure;
h1 = bar(1, mean(p1(vitGroupAll == 0 & ~isnan(p1))), 'facecolor', [.3 .3 .9]);
h2 = bar(2, mean(p1(vitGroupAll == 1 & ~isnan(p1))), 'facecolor', [.9 .3 .3]);

hold on;
h1 = bar(1, mean(p1(vitGroupAll == 0 & ~isnan(p1))), 'facecolor', [.3 .3 .9]);
e1 = errorbar(1,mean(p1(vitGroupAll == 0 & ~isnan(p1))), ...
    std(p1(vitGroupAll == 0 & ~isnan(p1))),'color', 'k', 'linewidth', 2);

% make width same size
h1 = get(e1, 'children');
x1 = get(h1, 'xdata');
x1(2) = {[1 1 NaN 0.97 1.03 NaN 0.97 1.03 NaN]};
set(h1(2), 'xdata', x1{2});
set(e1, 'children', h1);

h2 = bar(2, mean(p1(vitGroupAll == 1 & ~isnan(p1))), 'facecolor', [.9 .3 .3]);
e2 = errorbar(2,mean(p1(vitGroupAll == 1 & ~isnan(p1))), ...
    std(p1(vitGroupAll == 1 & ~isnan(p1))),'color', 'k', 'linewidth', 2);

% make width same size
h2 = get(e2, 'children');
x2 = get(h2, 'xdata');
x2(2) = {[2 2 NaN 1.97 2.03 NaN 1.97 2.03 NaN]};
set(h2(2), 'xdata', x2{2});
set(e2, 'children', h2);

set(gca, 'xtick', [1 2])
set(gca, 'fontsize', 14);
set(gca, 'xticklabel', {'Fresh', 'Vitrified'});
ylabel('Clumping Parameter');
title('Clumping Increases Post-Vitrification');
xlim([0.5 2.5]);
grid on;
legend('Vitrified', 'Fresh')


%% Maker kernel density plots

paramName = 'Clumping Parameter';
paramToEval = p1(~isnan(p1));
groundTruth = vitGroupAll(~isnan(p1));
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
title('Vitrification Increases Embryo Clumping');
xlabel(paramName);
ylabel('parameter histogram');

[f0, x0] = ksdensity(-1*paramToEval(groundTruth == 1-posClass));
[f1, x1] = ksdensity(-1*paramToEval(groundTruth == posClass));

hold on;
v0 = plot(x0, f0, 'Color', [0 0 1], 'linewidth', 2);
v1 = plot(x1, f1, 'Color', [1 0 0], 'linewidth', 2);
h = area(x0, f0, 'EdgeColor', [0 0 1], 'FaceColor', [0 0 1]);
hc = get(h, 'Children');
set(hc, 'FaceAlpha', .25);
h = area(x1, f1, 'EdgeColor', [1 0 0], 'FaceColor', [1 0 0]);
hc = get(h, 'Children');
set(hc, 'FaceAlpha', .25);
% line([0 .25], [0 0], 'Color', 'k', 'linewidth', 2);
% xlim([-.08 .08]);
legend([v0, v1], 'Fresh', 'Vitrified');
grid on;
legend('boxoff')





