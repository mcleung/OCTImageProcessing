function [] = disp3d(Rstore)
%DISP3D Summary of this function goes here
%   Detailed explanation goes here
colormap(jet);
for i=0:180
    imagesc(squeeze(Rstore(:,:,i+1))); %,[]);
    %imshow(squeeze(Rstore(:,120:920,i+1)),[]);
    
    colorbar; %colormap(gray);
    pause(0.01);
end
end

