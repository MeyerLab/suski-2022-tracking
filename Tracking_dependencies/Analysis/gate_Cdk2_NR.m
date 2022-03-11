function [CDK,badtracesCdk2]=gate_Cdk2_NR(nucsig,cytsig,maxthresh,noisethresh,smoothwindow)
numtraces=size(nucsig,1);
CDK=cytsig./nucsig;
% for i=1:numtraces
%     realframes=find(~isnan(nucsig(i,:)));
%     nucsigsmooth(i,realframes)=smooth(nucsig(i,realframes),smoothwindow);
%     CDKsmooth(i,realframes)=smooth(CDK(i,realframes),smoothwindow);
% end

%%% gate CDK2 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%hist(DHBnuc(:),0:10:2010); xlim([0 2000]);
%hist(max(DHBnuc,[],2),0:20:1020); xlim([0 1000]);

%%% remove traces where max nuc intensity is too low %%%%%%%%%%%%%%%%%%%%%%
lowsignal=max(nucsig,[],2)<=maxthresh; %as of 2014-01-30 was min
%%% remove ratio outliers %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sig2min=min(CDK,[],2);
sig2max=max(CDK,[],2);
%mu1=nanmean(DHBratio(:));
mumin=mean(sig2min);
mumax=mean(sig2max);
%iqr1=prctile(DHBratio(:),75)-prctile(DHBratio(:),25);
%outliers=sig2min<mu1-3*iqr1 | sig2max>mu1+3*iqr1;
iqrmin=prctile(sig2min,75)-prctile(sig2min,25);
iqrmax=prctile(sig2max,75)-prctile(sig2max,25);
outliers=sig2min<mumin-3*iqrmin | sig2max>mumax+3*iqrmax;
%%% remove noisy traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
noisy=max(diff(CDK,1,2),[],2)>noisethresh;
badtracesCdk2=lowsignal | outliers | noisy;
