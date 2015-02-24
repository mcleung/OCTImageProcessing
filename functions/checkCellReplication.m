function [ nCellF CellSel vals] = checkCellReplication( CellSel, nCell,metricStruct,peaks,TLVr )
%CHECKCELLREPLICATION Summary of this function goes here
%   Detailed explanation goes here
    %peaks=peaks(:,1:min(size(peaks,2),nCell*2));
    nCellF=nCell;
    nZeros=(unique(CellSel));
    nRef=length(metricStruct);
    CellSelProspect=zeros(nRef*size(peaks,2),1);
    CellSelProspect(1:(length(nZeros)-1))=nZeros(2:end);
for iRef=1:nRef
    metric=metricStruct.metric;
    for ipeak=1:size(peaks,2);
        dr=metric(1,ipeak);
        dd=metric(2,ipeak);
        Sdity=metric(3,ipeak);
        cond1=(abs(dr) < TLVr);
        cond2a=(dd>0.6*peaks(3,ipeak));
        cond2b=(dd<2*peaks(3,ipeak));
        cond3=Sdity>0.65;
        disp([cond1 cond2a cond2b cond3]);
        vals=[dr dd Sdity];
        if (cond1 & (cond2a & cond2b & cond3))
            if sum(ipeak==CellSel)==0
                nCellF=nCell+1;
                CellSel(nCellF)=ipeak;
            end
        else
            nCellF=nCellF;
        end
    end
    vals=[dr dd Sdity peaks(3,ipeak)];
    disp('-------');
end
%CellSel=flipu
end

