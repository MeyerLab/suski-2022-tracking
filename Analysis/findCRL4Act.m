function [ POI ] = findCRL4Act( CRLtrace, traceStats, settings )
%low, smooth, thresh, early, falseCall, buffer
[sampleSize,traceLength] = size(CRLtrace);
POI = ones(sampleSize,1)*NaN;
firstFrame = ones(sampleSize,1)*NaN;
lastFrame = ones(sampleSize,1)*NaN;
fullTraceStore = ones(sampleSize,traceLength)*NaN;
badTraces = false(sampleSize,1);

for i=1:size(CRLtrace,1)
    CRLact = NaN;
    if settings.cycling
        mitosisFrame = traceStats(i, 1);
        startFrame = mitosisFrame + settings.buff;
        endFrame = traceStats(i,2);
    else
        startFrame = traceStats(i,1);
        endFrame = traceStats(i,2);
    end
    
    firstFrame(i) = startFrame;
    lastFrame(i) = endFrame;
    fullTraceStore(i,:) = CRLtrace(i,:);
    
    if startFrame < endFrame - 5
        tracePost = CRLtrace(i,startFrame:endFrame);
        diffTest = diff(tracePost);
        diffTest = [diffTest diffTest(end)];
        
%         truncInd = find(tracePost < settings.low & diffTest < 0,1);
%         if ~isempty(truncInd) & truncInd > 1
%             tracePost = tracePost(1:truncInd);
%         end
%         tracePost = tracePost/max(tracePost); % renormalize
        
        if settings.smooth > 1
             %smoothTrace = nansmoothm(tracePost, settings.smooth,'sgolay');
              smoothTrace = smoothdata(tracePost, 'sgolay', [settings.smooth 0],'omitnan');
        else
            smoothTrace = tracePost;
        end
        
        
        
        firstD = gradient(smoothTrace);
        %firstD = [firstD firstD(end)];
        secD = gradient(firstD);
        %secD = [secD(1) secD];
        if settings.cycling
            filter = firstD < settings.firstD & secD < settings.secD & (1:length(smoothTrace)) > settings.early ...
                & tracePost > settings.falseCall;
        else
            filter = firstD < settings.firstD & secD < settings.secD & (1:length(smoothTrace)) + startFrame - 1 > settings.early ...
                & tracePost > settings.falseCall;
        end
        %%% secondary gate
        numFrames = size(tracePost,2);
        gate = zeros(1,numFrames);
        maxFutureHeight = zeros(1,numFrames);
        futureHeightDec = zeros(1,numFrames);
        allSlope = zeros(1,numFrames);
        allSec = zeros(1,numFrames);
        if length(smoothTrace) > settings.postBuffer & settings.postBuffer > 0
            for j = 1:numFrames - settings.postBuffer
                maxFutureHeight(j) = max(smoothTrace(j+1:j+settings.postBuffer)) < smoothTrace(j);
                futureHeightDec(j) = smoothTrace(j+settings.postBuffer) < smoothTrace(j) - settings.decrease; %default 0.1
                allSlope(j) = ~any(firstD(j:j+settings.postBuffer) > settings.firstD);
                allSec(j) = ~any(secD(j:j+settings.postBufferSec) > settings.secD);
                gate(j) =  maxFutureHeight(j) & futureHeightDec(j) & allSlope(j) & allSec(j);
            end
        elseif settings.postBuffer == 0
            gate = ones(1, numFrames);
        end
        
        CRLact= find(filter & gate,1,'first');
        if isempty(CRLact)
            POI(i) = NaN;
        else
            while CRLact > 2 && filter(CRLact-1)
               CRLact = CRLact-1;
            end
            
            POI(i) = CRLact + startFrame -1;
        end
        
    else
        tracePost = CRLtrace(i,:);
        firstD = 0;
        secD = 0;
        POI(i) = NaN;
    end
    
%    if any(i == [95 96 4])
%     figure(1),
%     subplot(3,1,1)
%     plot(tracePost), hold on
%     plot(smoothTrace,'--')
%     ylim([0 1]);
%     if ~isnan(CRLact)
%         scatter(CRLact, tracePost(CRLact),500,'r.');
%     end
%     hold off
%     subplot(3,1,2)
%     plot(firstD);
%     hline(settings.firstD);
%     
%     subplot(3,1,3)
%     plot(secD);
%     hline(settings.secD);
%      keyboard;
%    end
end
if settings.debug    
    %%% debug %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    traceIDs=[1:(96)];
    ylims=[0 1.2];
    POIs = {POI};
    POIdisplay(traceIDs,fullTraceStore, firstFrame, lastFrame,badTraces,ylims, POIs);  
    keyboard;
end

