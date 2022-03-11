function [signals,badtraces]=gate_lengthandrange_noise_NR(trace,tracestats,minlength,maxthreshold,minthreshold,maxnoisethresh,minnoisethresh)
% From gate_lengthandrange_noise_rev05.m
%%% smooth traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numtraces=size(trace,1);
signals=trace;
signals_for_diff=trace;
for i=1:numtraces
    realframes=find(~isnan(signals(i,:)));
    signals(i,realframes)=smooth(signals(i,realframes));
end
diffsignals=diff(signals,1,2);

%%% gate traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% shorttraces=tracestats(:,3)<=minlength;
nosignal=max(signals,[],2)<=maxthreshold;
strange=min(signals,[],2)>=minthreshold;
noisy_signal=max(diffsignals,[],2)>=maxnoisethresh;
noisy_signal2=min(diffsignals,[],2)<=minnoisethresh;
% badtraces=shorttraces | nosignal | strange;
badtraces= nosignal | strange | noisy_signal | noisy_signal2;

signals=trace; %return raw signal (not smoothened)