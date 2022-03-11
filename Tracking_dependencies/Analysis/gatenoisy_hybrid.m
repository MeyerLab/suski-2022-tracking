function [ badtraces ] = gatenoisy_hybrid( signal, trace_stats, isDaughter, noise_thresh, mitosis_buff, end_buff)
badtraces = zeros(size(signal,1),1);
curve_collect = [];
for i=1:size(signal, 1)
    
    startFrame = trace_stats(i,1);
    endFrame = trace_stats(i,2);
    trace = signal(i,startFrame:endFrame);
    if any(isnan(trace))
        trace = fillmissing(trace,'linear');
    end
    curvature = gradient(gradient(trace));
    gate = true(size(curvature));
    if length(gate) > end_buff
        gate(end-end_buff:end) = false;
    end
    if length(gate) > mitosis_buff & isDaughter(i)
        gate(1: 1+ mitosis_buff) = false;
    end
    curvature = curvature(gate);
    
    
    curve_collect = [curve_collect curvature];  % for determining threshold
    badtraces(i) = any(abs(curvature)>noise_thresh);
end

%histogram(curve_collect)
end

