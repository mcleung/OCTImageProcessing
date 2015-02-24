function [ metricStruct] = getMetric(BW,peakRef,peaks,nCell)
%GETMETRIC Summary of this function goes here
%   Detailed explanation goes here
nMetric=3;
metric=zeros(nMetric,size(peaks,2));  
for iRef=1:size(peakRef,2)
    peakRefSel=peakRef(:,iRef);
    for ipeak=1:size(peaks,2); %peaks(:,1:end)
        peakSel=peaks(:,ipeak);
        metric(1,ipeak)=peakRefSel(3,1)-peakSel(3,1);
        metric(2,ipeak)=sqrt((peakRefSel(1,1)-peakSel(1,1))^2+...
                        (peakRefSel(2,1)-peakSel(2,1))^2);
        BWcirclecrop=circleCrop(BW,peakSel);
        % cuts out edges
        
        BWSolidity=BWcirclecrop(max(peakSel(2)-peakSel(3),1):...
                                min(peakSel(2)+peakSel(3),size(BWcirclecrop,1)),...
                                max(peakSel(1)-peakSel(3),1):...
                                min(peakSel(1)+peakSel(3),size(BWcirclecrop,2)));
        solidit=(4/pi)*sum(sum(BWSolidity))/(size(BWSolidity,1)*size(BWSolidity,2));
        metric(3,ipeak)=solidit;
    end
    metricStruct(iRef).metric=metric;
end
end

