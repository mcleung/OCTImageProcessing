% Image Processing Pipeline for pre-processing OCT images
% Livia Zarnescu
% 1-9-14

%% Step 0: Load in Image Data

% add things to path
addPathRecursive('..');

% load in a single image stack
octDataPath = 'C:\Users\Livia\Desktop\OCT\OCT data\12-4-13 mouse thawed\TIFS';
imagePath = 'mouse_thawed_five_12_33_PM.tif';
imageData = imfinfo([octDataPath '\' imagePath]);
data = zeros(imageData(1,1).Height, imageData(1,1).Width, size(imageData,1)-5);

% read in all images in tiff stack
for i = 1:size(imageData,1)-5
    data(:,:,i) = imread([octDataPath '\' imagePath],i+5);
end

%% Step 1: Background Subtration

% gets rid of background = dirt on reference mirror
% input must be raw 16-bit data
disp('Step 1: Subtracting Background');
dataBackSub = subtractBackground(data, [1 size(data,3)]);
dataBackSub = im2double(uint16(dataBackSub));
fprintf('\b\n'); disp(' - 100 % Complete');

%% Step 2: auto adjust contrast / brightness, filter out noise

% define region in background to measure reference noise
% coords = [xmin xmax ymin ymax]
% use as data(ymin:ymax, xmin:xmax, i)
coords = getNoiseLocation(dataBackSub,0);

disp('Step 2: Adjust contrast/brightness');
dataAdj = AdjCBFilt(dataBackSub, coords);
fprintf('\b\n'); disp(' - 100 % Complete');

%% Step 3: Apply 3D Gaussian Blur

% 3D gaussian blur, sigma = 2
dataGauss=gauss3d(dataAdj,2,2,2);
figure, imshow(dataGauss(:,:,50));

%% Step 4: Get rid of remaining artifacts

% looks overly median filtered, but it's just to avoid
% confusing the edge detection 
dataFinal = removeLineArtifacts(data, coords);
figure, imshow(dataFinal(:,:,50))



















    
