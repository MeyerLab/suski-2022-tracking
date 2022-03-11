function [ POI ] = findCRL4Act_v2_reverse( data, settings )
%low, smooth, thresh, early, falseCall, buffer
CRLtrace = data.CRLtrace;
traceStats = data.traceStats;
if settings.cytoCheck
    cytoTrace = data.cyto;
end
nucTrace = data.nuc;

%%% Get start and end of frames based on tracestats
if settings.cycling
    mitosisFrame = traceStats(:, 1);
    startFrame = mitosisFrame + settings.buff;
    endFrame = traceStats(:,2);
else
    startFrame = traceStats(:,1);
    endFrame = traceStats(:,2);
end

%%% Mask traces, replacing values outside start and end with NaN
[numCells, time] = size(CRLtrace);
maskData= 1:time;
startMask = repmat(maskData, numCells,1) >= repmat(startFrame, 1, time);
endMask = repmat(maskData, numCells,1) <= repmat(endFrame, 1, time);
totalMask = startMask & endMask;
CRLtrace(~totalMask) = NaN;

% Smooth data
if settings.smooth > 1
    filtPoints = floor((settings.smooth - 1)/2);
    CRLsmooth = smoothdata(CRLtrace,2, 'sgolay',[filtPoints min([settings.postBuffer, filtPoints ])], 'omitnan','degree',3);
    CRLsmooth(~totalMask) = NaN;
else
    CRLsmooth = CRLtrace;
end

% Calculate derivatives 
firstD = [diff(CRLsmooth, [], 2) NaN(numCells,1)]; %forward 1st D
secD = [NaN(numCells,1) diff(firstD, [], 2)]; %Centered 2nd D
   
% Normalize
maxCRL = cummax(CRLsmooth,2);
normCRL = CRLsmooth ./ maxCRL;     
normFirstD = firstD ./ maxCRL;
normSecD = secD ./maxCRL;
     

% Calculate
POI = ones(numCells,2)*NaN;
fullTraceStore = ones(numCells,time)*NaN;
badTraces = false(numCells,1);

for i=1:size(CRLtrace,1)
    CRLact = NaN;
    fullTraceStore(i,:) = CRLtrace(i,:)/max(CRLtrace(i,:));
    
    if startFrame(i) < endFrame(i) - 5
        tracePost = CRLtrace(i,startFrame(i):endFrame(i));
        smoothTrace = CRLsmooth(i,startFrame(i):endFrame(i));
        normTrace = normCRL(i,startFrame(i):endFrame(i));
%         fD = firstD(i,startFrame(i):endFrame(i));
%         sD = secD(i,startFrame(i):endFrame(i));
        normFD = normFirstD(i,startFrame(i):endFrame(i));
        normSD = normSecD(i,startFrame(i):endFrame(i));
        nuc = nucTrace(i,startFrame(i):endFrame(i));
        
        filter = normFD < settings.firstD & normSD < settings.secD & normTrace > settings.falseCall &...
            (1:length(smoothTrace)) + ~settings.cycling * (startFrame(i) - 1) > settings.early & nuc > settings.minExp;
        
%         if settings.trunk & any(normTrace < settings.trunk)
%             filter = filter & (1:length(normTrace)) <= find(normTrace < settings.trunk,1,'first');
%         end
        if settings.cytoCheck
            cyto = cytoTrace(i,startFrame(i):endFrame(i));
            filter = filter & cyto < settings.cytoCheck;
        end
        
        
        
        %%% secondary gate
        numFrames = size(tracePost,2);
        gate = zeros(1,numFrames);
        maxFutureHeight = zeros(1,numFrames);
        futureHeightDec = zeros(1,numFrames);
        allSlope = zeros(1,numFrames);
        allSec = zeros(1,numFrames);
        allSec = zeros(1,numFrames);
        minSlope = zeros(1,numFrames);
        if length(smoothTrace) > settings.postBuffer & settings.postBuffer > 0
            for j = 1:numFrames - settings.postBuffer
                maxFutureHeight(j) = max(smoothTrace(j+1:j+settings.postBuffer)) < smoothTrace(j);
                futureHeightDec(j) = normTrace(j+settings.postBuffer) < normTrace(j) - settings.decrease; %default 0.1
                allSlope(j) = ~any(normFD(j:j+settings.postBuffer) > settings.firstD & normTrace(j:j+settings.postBuffer) > .25);
                allSec(j) = ~any(normFD(j:j+settings.postBufferSec) > settings.secD & normTrace(j:j+settings.postBufferSec) > .1);
                minSlope(j) = min(normFD(j:j+settings.postBuffer)) < settings.minSlope;
                gate(j) =  maxFutureHeight(j) & futureHeightDec(j) & allSlope(j) & allSec(j) & minSlope(j);
            end
        elseif settings.postBuffer == 0
            gate = ones(1, numFrames);
        end
        
        CRLact= find(filter & gate,1,'first');
        if isempty(CRLact)
            POI(i,1) = NaN;
        else
            while CRLact > 2 && filter(CRLact-1)
               CRLact = CRLact-1;
            end           
            POI(i,1) = CRLact + startFrame(i) -1;
            
            CRLactrev = CRLact;
            while CRLactrev > 2 && normFD(CRLactrev - 1) < 0
                CRLactrev = CRLactrev -1;
            end
            POI(i,2) = CRLactrev + startFrame(i) -1;
        end
        
    else
        tracePost = CRLtrace(i,:);
        firstD = 0;
        secD = 0;
        POI(i,1) = NaN;
    end
    
   if any(i == []) & settings.debug 
    figure(1),
    subplot(3,1,1)
    plot(tracePost), hold on
    plot(smoothTrace,'--')
%     ylim([0 1]);
    if ~isnan(CRLact)
        scatter([CRLact CRLactrev], tracePost([CRLact CRLactrev]),500,'r.');
    end
    hold off
    subplot(3,1,2)
    plot(normFD);
    hline(settings.firstD);
    
    subplot(3,1,3)
    plot(normSD);
    hline(settings.secD);
     keyboard;
   end
end
if settings.debug    
    %%% debug %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    numTraces = min([200 length(startFrame)]);
    traceIDs=[1:numTraces];
    ylims=[0 1.2];
    POIs = {POI(:,2)};
    POIdisplay(traceIDs,fullTraceStore, startFrame, endFrame,badTraces,ylims, POIs);  
    keyboard;
end

