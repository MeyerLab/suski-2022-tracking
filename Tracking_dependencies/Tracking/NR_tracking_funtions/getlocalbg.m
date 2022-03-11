function [ localbg ] = getlocalbg( img_masked, center, region , percentile,threshold)
% Calculates local background in a box around cell.
%LOCALBG = GETLOCALBG(IMG_MASKED,CENTER,REGION,PERCENTILE) returns the
%value with percentile PERCENTILE in a square with side length REGION
%around a cell centered at CENTER. Takes already masked image with NaNs in
%foreground region. If box is clipped at the edge of an image the box is 
%bounded. If the number of background pixels is less than 20% (can adjust
%value) of box size, multiply box size by 1.5 (can adjust). 

fg = img_masked;

center = round(center);  % location of cell (rounded to pixel)

% define box corners
regr = round(region/2);
xadd = [-1 1; -1 1]*regr;
yadd = [1 1; -1 -1]*regr;
xcorners = xadd + center(1);
ycorners = yadd + center(2);

% clip box if outside image region
xcorners(xcorners(:)<1)=1;
ycorners(ycorners(:)<1)=1;
xcorners(xcorners(:)>size(fg,2))=size(fg,2);
ycorners(ycorners(:)>size(fg,1))=size(fg,1);

%% plot box around cell in image
% figure(1), imagesc(fg)
% hold on
% plot(center(1),center(2), '.r')
% plot(xcorners(:),ycorners(:), '.y')
%%
localmat=fg(ycorners(2,1):ycorners(1,1),xcorners(1,1):xcorners(1,2));  %crops image only around cell
% figure(2), imagesc(localmat) %visualize cropped image

vals=localmat(~isnan(localmat(:)));  %vector of background values
% figure(3), histogram(vals)

localbg=prctile(vals,percentile);   %calculate percentile pixel in background
%vline(localbg)
% 
% Alternate method of calculating background (estimate background from kernel smoothed histogram) 
% binmax=prctile(vals,95);
% binmin=prctile(vals,5);
% vals=vals(vals<binmax & vals>binmin);
% [kval,xval]=ksdensity(vals);
% localbg2=xval(find(kval==max(kval),1,'first'));


%% Increases box size if not enough pixels
% threshold: fraction of box pixels which are background threshold
multiplier = 2;   %multiply box size by
if length(vals)<threshold*region^2
    regr=regr*multiplier;
    xadd=[-1 1; -1 1]*regr;
    yadd=[1 1; -1 -1]*regr;

    xcorners=xadd+center(1);
    ycorners=yadd+center(2);

    xcorners(xcorners(:)<1)=1;
    ycorners(ycorners(:)<1)=1;

    xcorners(xcorners(:)>size(fg,2))=size(fg,2);
    ycorners(ycorners(:)>size(fg,1))=size(fg,1);
    
    localmat=fg(ycorners(2,1):ycorners(1,1),xcorners(1,1):xcorners(1,2));
    
    vals=localmat(~isnan(localmat(:)));
    
    if length(vals)<threshold*(multiplier*region)^2
        localbg=NaN;
    else
        localbg=prctile(vals,percentile);
    end
    
end






    


end

