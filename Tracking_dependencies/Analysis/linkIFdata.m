function [IFdatalinked, gateIFMother]=linkIFdata(tracestats,samplecellsID,IFdata,gate)
gateIFMother=false(size(gate));
IFdatalinked = NaN * ones(size(IFdata));
for i=1:numel(samplecellsID)
    orgcellid=samplecellsID(i);
    cellid=orgcellid;
    condition=true;
    while condition
        IFdatalinked(cellid,:) = IFdata(orgcellid,:);
        if cellid ~= tracestats(cellid,4) & ~isnan(tracestats(cellid,4))
            cellid=tracestats(cellid,4);
            gateIFMother(cellid) = true;
        else
            condition = false;
        end
    end
   
end
