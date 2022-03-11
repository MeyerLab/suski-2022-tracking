function linkedtracedata=linkancestry_NR_full(tracedata,tracestats,samplecellsID)
linkedtracedata=tracedata;
for i=1:numel(samplecellsID)
    orgcellid=samplecellsID(i);
    cellid=orgcellid;
    condition=true;
    while condition
        goodframes=find(~isnan(tracedata(cellid,:,1)));
        linkedtracedata(orgcellid,goodframes,:)=tracedata(cellid,goodframes,:);
        if cellid ~= tracestats(cellid,4)
            cellid=tracestats(cellid,4);
            condition=~isnan(cellid);
        else
            condition = false;
        end
    end
   
end
