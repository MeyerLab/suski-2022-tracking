function [ badtraces] = gateSigTraces( traceSig, traceStats, buff, postMitosis, minMaxThresh, minMaxNoise,maxPost)
badtraces = zeros(size(traceSig,1),1);
noisy = zeros(size(traceSig,1),1);
lowMax = zeros(size(traceSig,1),1);
highMin = zeros(size(traceSig,1),1);

for i=1:size(traceSig, 1)
    if postMitosis == 1
        mitosisFrame = traceStats(i, 1);
        startFrame = mitosisFrame + buff;
        endFrame = traceStats(i,2);
        noisePost = 0;
        highMinPost= 0;
        if startFrame < endFrame - 1
            tracePost = traceSig(i,startFrame:endFrame);
            noisePost = (diff(tracePost));
            highMinPost = min(tracePost) > minMaxThresh(1);
        end
        if startFrame < endFrame - 1 & maxPost
            lowMax(i) = max(tracePost) < minMaxThresh(2);
        else
            lowMax(i) = max(traceSig(i,:)) < minMaxThresh(2);
        end
        highMin(i)= highMinPost;
        noisy(i) = any(noisePost > minMaxNoise(2)) | any(noisePost < minMaxNoise(1));
    else
        startFrame = traceStats(i,5);
        endFrame = traceStats(i,2);
        trace = traceSig(i,startFrame:endFrame);
        noise = (diff(trace));
        lowMax(i) = max(trace) < minMaxThresh(2);
        highMin(i) = min(trace) > minMaxThresh(1);
        noisy(i) = any(noise > minMaxNoise(2)) | any(noise < minMaxNoise(1));
    end

    
end
badtraces = lowMax | highMin | noisy;


end

%         figure, plot(traceGem(35:40,:)')
%         figure, plot(sensor(i).nuc_mean(35:40,:)')
%histogram(real(log2((sensor(i).geminin_nuc(:,:)'))))