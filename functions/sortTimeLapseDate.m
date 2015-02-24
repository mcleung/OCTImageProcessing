% first sort indices based on date
function [mOrder bOrder] = sortTimeLapseDate(dirData, month)

% indices of dirData with "M" or "bright" in name
mInd = [];
bInd = [];

% time elapsed for each element of mInd or bInd
mDate = [];
bDate = [];

for i = 1:size(dirData,1)
    
    currItem = dirData(i,1);
    i
    
    dateStart = strfind(currItem.name, ['_' num2str(month) '_']);
    
    if numel(dateStart > 0) & (~currItem.isdir)
        dateAndTime = currItem.name(dateStart(1)+1:end);
        
        currMo = str2double(dateAndTime(1));
        
        % indices of underscores or spaces
        underScores = strfind(dateAndTime, '_');
        spaces = strfind(dateAndTime, ' ');
        underScores = [underScores spaces];
        
        restParams = zeros(1, length(underScores)+1);
        restParams(1) = currMo;
        
        for j = 1:length(underScores)-1
            
            nextString = dateAndTime(underScores(j)+1:underScores(j+1)-1);
            restParams(j+1) = str2double(nextString);
            
        end
        
        if numel(restParams) == 7
            
            restParams(end) = (dateAndTime(underScores(end)+1) == 'P');
            
            % date sum = 0*year + 0*month + 24*day + hour + 12*(PM) +
            % + minute/60 + seconds/3600
            dateSum = 24*restParams(2) + mod(restParams(4),12) + ...
                12*restParams(7) + restParams(5)/60 + restParams(6)/3600;
            
            % if it ends in "M" or "bright"
            if currItem.name(end) == 'M'
                mInd = [mInd i];
                mDate = [mDate dateSum];
            elseif isequal(currItem.name(end-5:end),'bright')
                bDate = [bDate dateSum];
                bInd = [bInd i];
            end
        end
        
    end
    
end

% now sort

[~, mI] = sort(mDate, 'ascend');
mOrder = mInd(mI);

[~, bI] = sort(bDate, 'ascend');
bOrder = bInd(bI);

