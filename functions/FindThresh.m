% find background threshold

function thresh = FindThresh(background, pct)

backgroundR = reshape(background, 1, size(background,1)*size(background,2));
backgroundSort = sort(backgroundR, 'ascend');
thresh = backgroundSort(round(pct*length(backgroundSort)));
    
    
    