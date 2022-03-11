function [ badtraces ] = gatenoisy( signal, trace_stats, daughter_gate, quiescent, noise_thresh, buff, end_buff)
badtraces = zeros(size(signal,1),1);
curve_collect = [];
for i=1:size(signal, 1)
    if daughter_gate == 1 || ~quiescent
        mitosisFrame = trace_stats(i, 1);
        startFrame = mitosisFrame + buff;
        endFrame = trace_stats(i,2);
        if startFrame < endFrame - 1
            try
                trace = signal(i,startFrame:endFrame);
            catch
                keyboard
            end
            if any(isnan(trace))
                trace = fillmissing(trace,'linear');
            end
            curvature = gradient(gradient(trace));
            curvature = curvature((1+end_buff):(end-end_buff));
        else
            curvature = NaN;
        end
    elseif quiescent
        startFrame = trace_stats(i,1);
        endFrame = trace_stats(i,2);
        trace = signal(i,startFrame:endFrame);
        if any(isnan(trace))
            trace = fillmissing(trace,'linear');
        end
        curvature = gradient(gradient(trace));
        curvature = curvature((1+end_buff):(end-end_buff));
        
    end
    curve_collect = [curve_collect curvature];  % for determining threshold
    badtraces(i) = any(abs(curvature)>noise_thresh);
end

%histogram(curve_collect)
end

