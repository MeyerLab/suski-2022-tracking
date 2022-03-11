function [Sstart,badtraces]=getPCNAstart(sampletraces,samplestats,settings) 
initoffset = settings.initoffset;
punctathresh =settings.thresh;
punctalowthresh = settings.lowthresh;

%This function will return degstart assuming trace begins with mitosis
[samplesize,tracelength]=size(sampletraces);
firstframe=ones(samplesize,1)*NaN;
lastframe=ones(samplesize,1)*NaN;
Sstart=ones(samplesize,1)*NaN;
Sstartdb=Sstart;
badtraces=zeros(samplesize,1);
sigstore=ones(samplesize,tracelength)*NaN;

for i=1:samplesize
    signal=sampletraces(i,:);
    firstframe(i)=samplestats(i,1);
    lastframe(i)=samplestats(i,2);
    numframes=lastframe(i)-firstframe(i)+1;
%     signal=smooth(signal(firstframe(i):lastframe(i)),3)' ; %incoming signal always raw (unsmoothened)
%     signal = signal - min(signal);
    signal=signal(firstframe(i):lastframe(i)) - min(signal(firstframe(i):lastframe(i)));
    sigstore(i,firstframe(i):lastframe(i))=signal;
    
    %%%%%% general requirements %%%%%%%
    prebuffer=initoffset-1;
    gate=zeros(1,numframes);
    for j=initoffset:numframes
        pastmax=max(signal(j-prebuffer:j-1));
        presentgate=signal(j)>=punctathresh;
        gate(j)=signal(j)>1*pastmax & presentgate;
    end
    punctastart=find(gate,1,'first');
    lowerthresh=punctathresh*0.75;
    if isempty(punctastart) && signal(end)>lowerthresh
        badtraces(i)=1;
    elseif ~isempty(punctastart)
        %%% Find last zero before punctastart call %%%%%%%%%%%%%%%%%%%%%%%%
        earlypuncta=find(signal(1:punctastart-1)<=punctalowthresh,1,'last');
        if ~isempty(earlypuncta)
            punctastart=earlypuncta;
            Sstartdb(i)=punctastart;
            Sstart(i)=firstframe(i)+punctastart-1; %account for index
        end
    end
end

if settings.debug
traceIDs=1:96;
ylims=[0 200];
POIs = {Sstart};
POIdisplay(traceIDs,sigstore, firstframe, lastframe,badtraces,ylims, POIs);

keyboard;
end
