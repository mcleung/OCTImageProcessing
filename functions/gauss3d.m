function [ dataOut ] = gauss3d( data,sigmaX,sigmaY,sigmaZ)
%GAUSS3D Summary of this function goes here
%   Detailed explanation goes here
    hsizeX=20; %sigmaX=2;
    hsizeY=20; %sigmaY=2;
    hsizeZ=20; %sigmaZ=2;
    gaussFilterX = fspecial('gaussian', [hsizeX 1], sigmaX);
    gaussFilterY = fspecial('gaussian', [hsizeY 1] , sigmaY);
    gaussFilterZ = fspecial('gaussian', [hsizeZ 1], sigmaZ);
    Xfilt=zeros(size(data));
    Yfilt=zeros(size(data));
    Zfilt=zeros(size(data));
    %total=sum(sum(size(data)));
    for iX=1:size(data,1);
        for iY=1:size(data,2);
            Zfilt(iX,iY,:) = conv(squeeze(data(iX,iY,:)), gaussFilterZ, 'same');
        end
    end
    for iY=1:size(data,2);
        for iZ=1:size(data,3);
            Xfilt(:,iY,iZ) = conv(squeeze(Zfilt(:,iY,iZ)), gaussFilterX, 'same');
        end
    end
    for iX=1:size(data,1);
        for iZ=1:size(data,3);
            Yfilt(iX,:,iZ) = conv(squeeze(Xfilt(iX,:,iZ)), gaussFilterY, 'same');
        end
    end
    dataOut=Yfilt;
end

