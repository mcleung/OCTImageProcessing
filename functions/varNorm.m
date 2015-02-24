function dataSumOut = varNorm(dataSum)

% scale values in matrix from 0 to 1
% works on matrices up to 3D

dataSumOut = (dataSum - min(min(min(dataSum)))) / ...
    (max(max(max(dataSum))) - min(min(min(dataSum))));