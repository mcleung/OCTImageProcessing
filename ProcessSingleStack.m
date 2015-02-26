function [output] = ProcessSingleStack(fileNameList, dirName)
% Process single stack into 16-bit TIF

close all;

nPixelX = 860;
nPixelY = 800;

if ~exist([dirName '\TIFS'],'dir')
    mkdir(dirName, 'TIFS');
    mkdir([dirName '\TIFS'], 'Filtered');
    mkdir([dirName '\TIFS\Filtered'],'GaussBlur');
elseif ~exist([dirName '\TIFS\Filtered'],'dir')
    mkdir([dirName '\TIFS'], 'Filtered');
    mkdir([dirName '\TIFS\Filtered'],'GaussBlur');
elseif ~exist([dirName '\TIFS\Filtered\GaussBlur'], 'dir')
    mkdir([dirName '\TIFS\Filtered'],'GaussBlur');
end

for ii = 1:length(fileNameList)

    fileName = fileNameList{ii};
    dataIn = load([dirName '\' fileName]);
    nPixelZ = size(dataIn,1)/nPixelX;
    
    % convert 2d table to 3d array
    data = zeros(nPixelX, nPixelY, nPixelZ);
    for j=1:nPixelZ
        data(:,:,j) = dataIn(((j-1)*nPixelX + 1):(j*nPixelX),:);
    end
    
    clear dataIn;
    data = uint16(data);
    
    % write as 16-bit tif stack
    n = 1;
    fileNameWrite = [dirName '\TIFS\' fileName];
    
    imwrite(data(:,:,n), [fileNameWrite '.tif']);
    
    for k = n+1:size(data,3)
        imwrite(data(:,:,k), [fileNameWrite '.tif'], ...
            'writemode', 'append');
    end
    
    % Or, just load 16-bit tif stack
    % post-processing
    
    % Step 0: Data preprocessing
    dataSub = data;
    dataSub2 = data;
    
    % Step 1: subtract background, convert to double between 0 and 1
    % gets rid of background = pillar like corruption
    disp('Step 1: Subtracting Background');
    dataBackSub = subtractBackground(data, [5 size(data,3)-5]);
    dataBackSub = im2double(uint16(dataBackSub));
    fprintf('\b'); disp(' - 100 % Complete');
    figure, imagesc(sum(dataBackSub(:,:,10:end),3));
    
    % Step 2: auto adjust contrast / brightness, filter out noise
    % define region in background to measure reference noise
    % coords = [xmin xmax ymin ymax]
    % use as data(ymin:ymax, xmin:xmax, i)
    coords = getNoiseLocation(data);

    
    disp('Step 2: Adjust contrast/brightness');
    dataAdj = AdjCBFilt(dataBackSub, coords);
    fprintf('\b'); disp(' - 100 % Complete');
    
    % Step 3:  3D gaussian blur
    % Use CONVN three times to filter your data with three 1D Gaussians...
    %    , one x-by-1-by-1, one 1-by-y-by-1, and one 1-by-1-by-z.
    % syntax: h = fspecial('gaussian', hsize, sigma)
    sigmaX=3; sigmaY=3; sigmaZ=3;
    disp('Step 3: Applying Gaussian Blur');
    dataGauss=gauss3d(dataAdj,sigmaX,sigmaY,sigmaZ);
    fprintf('\b'); disp(' - 100 % Complete');

    % %% Write in file
    % save(fNameSave,'dataGauss');
    % dirData(j).dataGauss=dataGauss;
    %
    % write as tif stack, leave out dish surface
    n = 3;
    dataAdj16 = uint16(dataAdj*65535);
    % Save Gauss blurred stack
    dataGauss16 = uint16(dataGauss*65535);
    
    imwrite(dataAdj16(:,:,n), [dirName '\TIFS\Filtered\' fileName '.tif']);
    
    for k = n+1:size(dataAdj,3)
        imwrite(dataAdj16(:,:,k), [dirName '\TIFS\Filtered\' fileName '.tif'], ...
            'writemode', 'append');
    end
    
    imwrite(dataGauss16(:,:,n), [dirName '\TIFS\Filtered\GaussBlur\' fileName '.tif']);
    
    for k = n+1:size(dataAdj,3)
        imwrite(dataGauss16(:,:,k), [dirName '\TIFS\Filtered\GaussBlur\' fileName '.tif'], ...
            'writemode', 'append');
    end
    
    fprintf('\nDONE\n');
    
end
end
