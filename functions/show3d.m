function [ Rstore ] = show3d( data )
%SHOW3DD Summary of this function goes here
%   Detailed explanation goes here
%% Projection
h.waitbar=waitbar(0);
theta=0:180;
%Rstore=zeros(115,1041,181);
for iStack=1:115
    R = radon(squeeze(data(:,:,iStack)),theta);
    if iStack==1
        Rstore=zeros(115,size(R,1),181);
    end
    Rstore(iStack,:,:)=R;
    if mod(iStack,20)==0 waitbar(iStack/115); end; 
end
end

