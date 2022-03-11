function [image1,image2]=cropboth(image1,image2,reljitx,reljity)
if reljitx>0
    image1(:,1:reljitx)=[];
    image2(:,end-reljitx+1:end)=[];
elseif reljitx<0
    image1(:,end+reljitx+1:end)=[];
    image2(:,1:-reljitx)=[];
end
if reljity>0
    image1(1:reljity,:)=[];
    image2(end-reljity+1:end,:)=[];
elseif reljity<0
    image1(end+reljity+1:end,:)=[];
    image2(1:-reljity,:)=[];
end