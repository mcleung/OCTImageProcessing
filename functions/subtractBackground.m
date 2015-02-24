% this function detects the dirt on the reference mirror and removes it

function dataSub2 = subtractBackground(data,ULB)

% maybe look at histogram to choose thresh
% Image sum up
dataMean = mean(data(:,:,end-5:end),3);
dataMeanScaledFilt = medfilt2(dataMean, [15 15]);

% compute background
dataBackground = abs(dataMean - dataMeanScaledFilt);
% dilate background image to get peaks
dataBackDilate = imdilate(dataBackground, ones(3,3));
% convert to BW to get neighborhoods of brightness peaks
dataBackThresh = zeros(size(dataBackDilate,1), size(dataBackDilate,2));
thresh = mean(mean(dataBackDilate))+ 3*std(std(dataBackDilate));
dataBackThresh(dataBackDilate > thresh) = 1;
dataSub = zeros(size(data,1),size(data,2),ULB(2)-ULB(1)+1);
%h.waitbar=waitbar(0);
counter=1;
for i = ULB(1):ULB(2)
    if mod(i,10)==0
        if i~=10
            fprintf('\b\b\b\b\b\b\b\b');
        end
        fprintf('\n %f percent complete' , ceil(i/size(data,3)*100));
    end
%    i 
    currFrame = data(:,:,i);
    currFrameR = reshape(currFrame, 1, size(currFrame,1)*size(currFrame,2));
    currFrameR = sort(currFrameR, 'ascend');
    
    % threshold at 90% of max intensity value
    currThresh = currFrameR(round(length(currFrameR)*.9));
    currMed = medfilt2(wiener2(currFrame),[10 10]);
    
    % take pixels that are brighter than 90% of the oi = ther pixels. If they
    % are close to the bright spots in the background image, set them to
    % the values of the median filtered image.
    currSub = currFrame;
    currSub(currFrame > currThresh & dataBackThresh == 1) = ...
        currMed(currFrame > currThresh & dataBackThresh == 1);
    dataSub(:,:,counter) = currSub;
    counter=counter+1;
    %if mod(i,10)==0 waitbar(i/size(data,3)); end; 
end

allSum = sum(dataSub(:,:,1:end),3);
dataSub2 = data(:,:,ULB(1):ULB(2));

% get rid of edge effects (corners are really dark but shouldnt be)
counter=1;
for i = ULB(1):ULB(2)
    
    currSub = dataSub(:,:,i-ULB(1)+1);
    currSub(allSum < 100) = mean(mean(allSum))/(size(dataSub,3));
    dataSub2(:,:,i-ULB(1)+1) = currSub;
    
end
%delete(h.waitbar);
