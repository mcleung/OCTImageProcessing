
% load params from all dates

function [stdevList entropyVals IntensityRatio1 IntensityRatio2] = ...
    loadAllParams(exptsToLoad, embryoGroup, cellClump, posClass)


stdevList = zeros(1,length(embryoGroup));
entropyVals = zeros(1,length(embryoGroup));
IntensityRatio1 = zeros(1,length(embryoGroup));
IntensityRatio2 = zeros(1,length(embryoGroup));

vitGroup = cell(1,exptsToLoad);
cellClump = cell(1,exptsToLoad);
blastForm = cell(1,exptsToLoad);
numEmbryos = cell(1,exptsToLoad); % number of embryos

% 1 = vit, 0 = non-vit
vitGroup{1} = [zeros(1,4) ones(1,6)]; % 4 non-vit from 4-26-13 and 6 vit from 5-3-13
vitGroup{2} = [zeros(1,5) ones(1,5)]; % 5 non-vit and 5 vit from 12-4-13
vitGroup{3} = [zeros(1,5) ones(1,5)]; % 5 non-vit and 5 vit from 1-14-14
vitGroup{4} = [zeros(1,5) ones(1,10)]; % 5 non-vit and 10 vit from 3-27-14

blastForm{1} = NaN*ones(1,10);
blastForm{2} = NaN*ones(1,10);
blastForm{3} = [0 0 1 1 1 0 0 0 1 1];
blastForm{4} = 


% load experiment from 1-14-14
if exist('embryoParams.mat', 'file')
    load('embryoParams.mat');
end

for i = 1:10
    
    i
    
    if cellClump(i) == posClass
        currColor = [0 0 .6];
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
        
    entropyVals(i) = mean(eval(['Group' num2str(groupNum) '.E' ...
        num2str(embryoNum) '.entropyfilt;']));
    stdevList(i) = mean(eval(['Group' num2str(groupNum) '.E' ...
        num2str(embryoNum) '.stdev;']));
    
    h1 = eval(['Group' num2str(groupNum) '.E' num2str(embryoNum) ...
        '.rProfile']);
    
    for j = 1:size(h1,1)
        h2 = h1(j,1:5);
        h4 = h1(j,end-5:end);
        h3(j) = mean(h2(~isnan(h2)));
        h5(j) = mean(h4(~isnan(h4)));
    end
        

    IntensityRatio1(i) = max(h5)/min(h5);
    IntensityRatio2(i) = max(h3)/min(h3);    
    
    figure(1);
    hold on;
    set(gca, 'fontsize', 14);
    plot(h3 - min(h3), 'Color', currColor, 'LineWidth', 2);
    xlabel('Normalized distance from cell center');
    ylabel('Intensity');
    title('Intensity Profile Along Cell Radius');
    xlim([1 10]);
    
    figure(2);
    hold on;
    set(gca, 'fontsize', 14);
    plot(h5, 'Color', currColor, 'LineWidth', 2);
    xlabel('Normalized distance from cell center');
    ylabel('Intensity');
    title('Intensity Profile Along Cell Radius');
    xlim([1 10]);

end

if length(embryoGroup) < 11
    return;
end

clearvars -except stdevList entropyVals IntensityRatio1 IntensityRatio2 ...
    embryoGroup cellClump posClass

% load experiment from 12-4-13
if exist('embryoParams2.mat', 'file')
    load('embryoParams2.mat');
end

for i = 1:10
    
    i
    
    if cellClump(i) == posClass
        currColor = [0 0 .6];
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
        
    entropyVals(i+10) = mean(eval(['Group' num2str(groupNum) '.E' ...
        num2str(embryoNum) '.entropyfilt;']));
    stdevList(i+10) = mean(eval(['Group' num2str(groupNum) '.E' ...
        num2str(embryoNum) '.stdev;']));
    
    h1 = eval(['Group' num2str(groupNum) '.E' num2str(embryoNum) ...
        '.rProfile']);
    
    for j = 1:size(h1,1)
        h2 = h1(j,1:5);
        h4 = h1(j,end-5:end);
        h3(j) = mean(h2(~isnan(h2)));
        h5(j) = mean(h4(~isnan(h4)));
    end
        

    IntensityRatio1(i+10) = max(h5)/min(h5);
    IntensityRatio2(i+10) = max(h3)/min(h3);     
    
    figure(1);
    hold on;
    plot(h3 - min(h3), 'Color', currColor, 'LineWidth', 2);
    
    figure(2);
    hold on;
    plot(h5, 'Color', currColor, 'LineWidth', 2);
    
end


if length(embryoGroup) < 21
    return;
end

clearvars -except stdevList entropyVals IntensityRatio1 IntensityRatio2 ...
    embryoGroup cellClump posClass

% load experiment from 4-26-13 and 5-3-13
if exist('embryoParams3.mat', 'file')
    load('embryoParams3.mat');
end

for i = 1:10
    
    i
    
    if cellClump(i) == posClass
        currColor = [0 0 .6];
    else
        currColor = [1 0 0];
    end
    
    if i < 5
        groupNum = 1;
        embryoNum = i;
    else
        groupNum = 2;
        embryoNum = i - 4;
    end
        
    entropyVals(i+20) = mean(eval(['Group' num2str(groupNum) '.E' ...
        num2str(embryoNum) '.entropyfilt;']));
    stdevList(i+20) = mean(eval(['Group' num2str(groupNum) '.E' ...
        num2str(embryoNum) '.stdev;']));
    
    h1 = eval(['Group' num2str(groupNum) '.E' num2str(embryoNum) ...
        '.rProfile']);
    
    for j = 1:size(h1,1)
        h2 = h1(j,1:5);
        h4 = h1(j,end-5:end);
        h3(j) = mean(h2(~isnan(h2)));
        h5(j) = mean(h4(~isnan(h4)));
    end
        

    IntensityRatio1(i+20) = max(h5)/min(h5);
    IntensityRatio2(i+20) = max(h3)/min(h3);     

    figure(1);
    hold on;
    set(gca, 'fontsize', 14);
    plot(h3 - min(h3), 'Color', currColor, 'LineWidth', 2);
    xlabel('Normalized distance from cell center');
    ylabel('Intensity');
    title('Intensity Profile Along Cell Radius');
    xlim([1 10]);
    
    figure(2);
    hold on;
    set(gca, 'fontsize', 14);
    plot(h5, 'Color', currColor, 'LineWidth', 2);
    xlabel('Normalized distance from cell center');
    ylabel('Intensity');
    title('Intensity Profile Along Cell Radius');
    xlim([1 10]);
    
end















