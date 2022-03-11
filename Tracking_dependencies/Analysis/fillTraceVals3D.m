function [ traceFilled, badtrace] = fillTraceVals3D( trace, traceStats, nanlim )
traceFilled = NaN*ones(size(trace));
badtrace = zeros(size(trace,1),1);
for i = 1:size(trace, 1)
        startFrame = traceStats(i,5);
        endFrame = traceStats(i,2);
        for j = 1:size(trace,3)
            if any(isnan(trace(i,startFrame:endFrame,j)))
                if isnan(trace(i,startFrame,j))
                    startFrame = find(~isnan(trace(i,:,j)),1,'first');
                end
                if isnan(trace(i,endFrame,j))
                    endFrame = find(~isnan(trace(i,:,j)),1,'last');
                end
                traceFilled(i, startFrame:endFrame,j) = fillmissing(trace(i,startFrame:endFrame,j),'linear');
                if sum(isnan(trace(i,startFrame:endFrame,j))) > nanlim
                    badtrace(i) = 1;
                end
            else
                traceFilled(i, startFrame:endFrame,j)  = trace(i,startFrame:endFrame,j);
            end
        end
end

end


