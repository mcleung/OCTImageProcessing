function Afilt = medfiltL(A, windowSize)

offset = (windowSize - 1)/2;
Afilt = A;

% loop thru all pixels of image
% i = 1:564
% j = 1:640

for i = 1:size(A,1)
   
    % determine offsets for pixels close to the edges
    xoLow = 0;
    xoHigh = 0;
    if i < offset+1
        xoLow = offset+1 - i;
    end
    if i > size(A,1) - offset
        xoHigh = i - size(A,1) + offset;
    end
    
    for j = 1:size(A,2)
        
        % determine offsets for pixels close to the edges
        yoLow = 0;
        yoHigh = 0;
        if j < offset+1
            yoLow = offset+1 - j;
        end
        if j > size(A,2) - offset
            yoHigh = j - size(A,2) + offset;
        end
        
        % determine window and weights matrices
        imWindow = A(i-offset+xoLow:i+offset-xoHigh,...
            j-offset+yoLow:j+offset-yoHigh);
        
        % populate vector of which to find median
        W = imWindow(:);
        
        Afilt(i,j) = median(W);
        
    end
end

1
