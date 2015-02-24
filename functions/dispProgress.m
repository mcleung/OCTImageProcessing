function [ ] = dispProgress(n,N)
%DISPPROGRESS Summary of this function goes here
%   Detailed explanation goes here
percent=ceil(n/N*100);
if percent<10 
        textD=['00' num2str(percent) ' %']; 
elseif percent>=10 & percent<100
        textD=['0' num2str(percent) ' %'];
else
        textD=[num2str(percent) ' %'];
end
fprintf('\b\b\b\b\b\b'); disp(textD);
end

