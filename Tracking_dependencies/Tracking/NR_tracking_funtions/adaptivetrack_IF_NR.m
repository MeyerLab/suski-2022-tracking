function [tracked,nuc_label, areadiffrec]=adaptivetrack_IF_NR(prevdata,curdata,nuc_label,nucr,distthresh,arealow,areahigh,debugpackage)
%NOTE: same as adaptivetrack, but mass requirements loosened
%%% set up %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
winrad=2*nucr;          %farthest jump (pixels)
%masschangethreshold=0.20;
xprev=prevdata(:,1); yprev=prevdata(:,2); areaprev=prevdata(:,3); massprev=prevdata(:,4);
numprevorg=numel(xprev);     %this is the number of traces from beginning of movie
previd=find(~isnan(xprev));
numprev=numel(previd);
xprev=xprev(previd); yprev=yprev(previd); areaprev=areaprev(previd); massprev=massprev(previd);
xcur=curdata(:,1); ycur=curdata(:,2); areacur=curdata(:,3); masscur=curdata(:,4);
numcur=numel(xcur);
%%% match each cell from previous frame to a cell in the current frame %%%%
prevmatch=ones(numprev,1)*NaN;
curmatch=zeros(numcur,1); %number of times matched
areadiffrec=ones(numprev,1)*NaN;
noneigh=0;
for i=1:numprev
%     if previd(i) == 2037
%         keyboard;
%     end
    neighbors=find(abs(xcur-xprev(i))<winrad & abs(ycur-yprev(i))<winrad);
    %%% in case of zero neighbors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isempty(neighbors)
        noneigh=noneigh+1;
        continue;
    end
    %%% get distances and mass differences of neighbors %%%%%%%%%%%%%%%%%%%
    dist=sqrt((xcur(neighbors)-xprev(i)).^2+(ycur(neighbors)-yprev(i)).^2);
    areadiff=(areacur(neighbors)-areaprev(i))/areaprev(i);
    didx=find(dist==min(dist));
    matchid=neighbors(didx);
    %if dist(didx)<2*nucr && areadiff(didx)>-0.1 && areadiff(didx)<0.3
    %if dist(didx)<distthresh && abs(areadiff(didx))<areachangethreshold
    if dist(didx)<distthresh && areadiff(didx)>arealow && areadiff(didx)<areahigh
        prevmatch(i)=matchid;
        curmatch(matchid)=curmatch(matchid)+1;
    end
    areadiffrec(i)=areadiff(didx);
end
%%% check continuity of area %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%hist(areadiffrec,-2.01:0.01:2.01),xlim([-2 2]);
%%% remove conflicts %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
curconflict=find(curmatch>1);
conflicts = 0;
for c=curconflict'
    tempprevs=prevmatch==c;
    prevmatch(tempprevs)=NaN;
    curmatch(c)=0;
    conflicts = conflicts + 1;
end
%%%%%% add tracked cells %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tracked=ones(numprevorg,4)*NaN;
tempidx=find(~isnan(prevmatch));
matchidxprev=previd(tempidx);
matchidxcur=prevmatch(tempidx);
tracked(matchidxprev,:)=[xcur(matchidxcur) ycur(matchidxcur) areacur(matchidxcur) masscur(matchidxcur)];
%%%%%% re-label nuc_label %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
relabelidx=zeros(numcur,1);
relabelidx(matchidxcur)=matchidxprev;
nuc_info=regionprops(nuc_label,'PixelIdxList');
for i=1:numcur
    nuc_label(nuc_info(i).PixelIdxList)=relabelidx(i);
end

end
%%% debug %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
previdx=3767;
prevcellandcurrentneighborsIF(debugpackage,winrad,xprev,yprev,previdx);
%}