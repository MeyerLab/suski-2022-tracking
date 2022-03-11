function thicklabelimg = labelthicken_better(labeledimg,thickenradius)
%Fixed from labelthicken.m by Yilin
[height,width]=size(labeledimg);
thickenedmask=bwmorph(labeledimg,'thicken',thickenradius);
% centers=round(squeeze(cell2mat(struct2cell(regionprops(thickenedmask,'pixelidx')')))');
% centers_idx=sub2ind([height,width],centers(:,2),centers(:,1));
% center_label=labeledimg(centers_idx);

pixels=regionprops(thickenedmask,'PixelIdxList');
numobjects=size(pixels,1);
thicklabelimg=zeros(height,width);
for i=1:numobjects
    tempidx = labeledimg(pixels(i).PixelIdxList);
    tempidx(tempidx==0) = nan;
    thicklabelimg(pixels(i).PixelIdxList) = nanmode(tempidx);%center_label(i);
end
end