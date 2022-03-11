function [ traceFilled, badtrace] = fillTraceVals( trace, traceStats, nanlim )
traceFilled = NaN*ones(size(trace,1),size(trace,2));
badtrace = zeros(size(trace,1),1);
for i = 1:size(trace, 1)
        startFrame = traceStats(i,5);
        endFrame = traceStats(i,2);
    if any(isnan(trace(i,startFrame:endFrame)))
        if isnan(trace(i,startFrame))
            startFrame = find(~isnan(trace(i,:)),1,'first');
        end
        if isnan(trace(i,endFrame))
            endFrame = find(~isnan(trace(i,:)),1,'last');
        end
        traceFilled(i, startFrame:endFrame) = fillmissing(trace(i,startFrame:endFrame),'linear');
        if sum(isnan(trace(i,startFrame:endFrame))) > nanlim
            badtrace(i) = 1;
        end
    else
        traceFilled(i, startFrame:endFrame)  = trace(i,startFrame:endFrame);
    end
end

end


