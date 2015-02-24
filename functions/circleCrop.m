function [ BWCircleCrop] = circleCrop(BW,peakSel)
%CIRCLECROP Summary of this function goes here
%   Detailed explanation goes here
[xN yN] = size(BW);
BWCircleCrop=BW;
r2=peakSel(3)^2;
for ix=1:xN
    for iy=1:yN
        cR2=(iy-peakSel(1))^2+(ix-peakSel(2))^2;
        if cR2>r2
            BWCircleCrop(ix,iy)=0;
        end
    end
end
  %[x, y] = circlepoints(peakSel(3));
  %imshow(BWCircleCrop);
  %plot(x+peakSel(1), y+peakSel(2), 'g-');
end

