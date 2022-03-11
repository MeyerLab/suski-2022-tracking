function [ImageMask,th,bg]=ThreshImage(ImageOld)

TempSeries=ImageOld(:);

%%% get background value: search intensity distribution for most convex point %%%%%%%%%%%%%%%%%%%%%%%
tmax=max(TempSeries);
tmin=min(TempSeries);
nbin=200;tbin=(tmax-tmin)/nbin;
tmin=tmin+tbin/2;tmax=tmax-tbin/2;
[n,xout]=ksdensity(TempSeries,tmin:tbin:tmax);
gp=max([2,ceil(nbin/50)]);  %gp=4
ng=getcurvature(n,gp);      %returns the angle change of intensity histogram over one bin at each point

negpeak=regionprops(bwlabel(ng<0),'PixelIdxList','Area'); %returns indices and area (# of indices) for each contiguous set of bins with negative (convex) curvature
negth=prctile(ng,25)-1.5*iqr(ng);                       %threshold for low outliers
Ibg00=zeros(1,size(negpeak,1));AC=Ibg00;DD=Ibg00;
for cc=1:size(negpeak,1)
    AC(cc)=negpeak(cc).Area;
    tempeakvalues=ng(negpeak(cc).PixelIdxList);         %curvature for each index of negpeak group
    [DD(cc),I]=min(tempeakvalues);                      %most negative curvature
    Ibg00(cc)=negpeak(cc).PixelIdxList(I);              %index for this most negative curve for this negpeak group
end
neglo=((DD<=negth)|(AC>=gp));       %each peak legitimate if: curvature is low outlier or area > 4 bins
Ibg00=Ibg00(neglo);DD=DD(neglo);    %remove illegitimate peaks
if isempty(Ibg00)
    [DD,Ibg00]=min(ng);             %if set of peaks is empty, set to index with lowest global curvature
end
Ibg00=[Ibg00,size(ng,2)-gp];        %pad row vector with last intensity bin (for threshold calc contingency)
[dummy,I_DD]=min(DD);               %I_DD = negpeak group # with the lowest global curvature
bg0=xout(Ibg00(I_DD));              %intensity bin of peak with lowest global curvature

TempSeries2=TempSeries((TempSeries>(bg0-5*gp*tbin))&(TempSeries<(bg0+5*gp*tbin)));  %narrow search for peak
[n_bg,xout_bg]=ksdensity(TempSeries2);
[dummy,Ibg]=max(n_bg);              %index of peak
xfit=xout_bg(Ibg-1:Ibg+1)';yfit=n_bg(Ibg-1:Ibg+1)';
det_n=det([xfit.^2,xfit,[1;1;1]]);
a=det([yfit,xfit,[1;1;1]])/det_n;
b=det([xfit.^2,yfit,[1;1;1]])/det_n;
bg=(-b/a/2);                          %refined index of peak --> background intensity

%%% get threshold: search upstream of background for most concave point %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
upperbound = xout(Ibg00(I_DD+1))+10*gp*tbin;
lowerbound = bg-gp*tbin;
TempSeries3=TempSeries((TempSeries>lowerbound)&(TempSeries<upperbound));
[n_th,xout_th]=ksdensity(TempSeries3,bg:(upperbound-bg)/100:upperbound);    %shouldn't ng be lowerbound?
[dummy,Ibg]=min(abs(xout_th-bg));   %index at background
ng_th=getcurvature(n_th,gp);
partng=ng_th(Ibg:end);              %curvature array starting from background
pospeak=regionprops(bwlabel(partng>0),'PixelIdxList','Centroid','Area');    %same as negpeak, but for positive (concave) curvature
Ith00=zeros(1,size(pospeak,1));AC=Ith00;
for cc=1:size(pospeak,1)
    AC(cc)=sum(partng(pospeak(cc).PixelIdxList));   %sum of curvatures for each pospeak group
    %AC(cc)=pospeak(cc).Area;
    Ith00(cc)=round(pospeak(cc).Centroid(1));       %index of center of each pospeak group
    %Ith00(cc)=(sum((pospeak(cc).PixelIdxList)'.*(partng(pospeak(cc).PixelIdxList)))/AC(cc));
end
[dummy,I]=max(AC);Ith0=floor(Ith00(I));Ithr=(Ith00(I)-Ith0);
s1=xout_th(Ith0+Ibg-1);s2=xout_th(Ith0+Ibg);
th=s1+Ithr*(s2-s1);

%%% get mask %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ImageMask=single(ImageOld>th);