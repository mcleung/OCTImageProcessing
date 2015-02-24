function [ peakRef CellSel] = getPeakRef(nCell,CellSel,peakRefOld,peaks )
%GETPEAKREF Summary of this function goes here
%   Detailed explanation goes here
peakRef=peakRefOld;
for iRef=1:size(peakRefOld,2)
    peakRefSel=peakRefOld(:,iRef);
    peakRefSel=repmat(peakRefSel,1,size(peaks,2));
        proximity1=[sqrt(...
                         (peakRefSel(1,:)-peaks(1,:)).^2 ...
                       + (peakRefSel(2,:)-peaks(2,:)).^2 ...
                          )];
        proximity2=[abs(peakRefSel(3,:)-peaks(3,:))];
        proximity=proximity1+proximity2;
        [val Cellind]=min(proximity);
        peakRef(:,iRef)=peaks(:,Cellind);
        CellSel(iRef)=Cellind;
end
end

