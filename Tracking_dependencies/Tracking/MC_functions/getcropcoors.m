function [cropcoors]=getcropcoors(dims,jitx,jity)
numshifts=numel(jitx);
cropcoors=ones(numshifts+1,4)*NaN;
%%% get maximum shifts of each side of image relative to first image %%%%%%
maxtopshift=max([0;jity]);
maxbottomshift=min([0;jity]);
maxleftshift=max([0;jitx]);
maxrightshift=min([0;jitx]);
%%% get cropping coordinates for each shifted frame %%%%%%%%%%%%%%%%%%%%%%%
cropcoors(1,:)=[1+maxtopshift dims(2)+maxbottomshift 1+maxleftshift dims(1)+maxrightshift];
for i=1:numshifts
    topshift=maxtopshift-jity(i);
    bottomshift=jity(i)-maxbottomshift;
    cropcoors(i+1,1:2)=[1+topshift dims(2)-bottomshift];
    leftshift=maxleftshift-jitx(i);
    rightshift=jitx(i)-maxrightshift;
    cropcoors(i+1,3:4)=[1+leftshift dims(1)-rightshift];
end