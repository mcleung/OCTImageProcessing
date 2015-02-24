% implement 3D median filter
% windowSize is single odd number since window will be symmetric in 3D

function im_out = medfilt3D(im, windowSize)

offset = round((windowSize - 1) / 2);
im_out = im;

for i = 1+offset:size(im,1)-offset
    i
    
    for j = 1+offset:size(im,2)-offset
        
        for k = 1+offset:size(im,3)-offset
            
            currWindow = im(i-offset:i+offset, j-offset:j+offset, ...
                k-offset:k+offset);
            
            currWindowSort = sort(currWindow(:), 'ascend')';
            im_out(i,j,k) = currWindowSort(round((windowSize^2 + 1)/2));
            
        end
    end
end

