function mask=threshmask_multi(image,blurradius,numthresh,toosmall)
blur=imfilter(image,fspecial('disk',blurradius),'symmetric'); %10x:3 20x:6
normlog=mat2gray(log(blur));
thresh=multithresh(normlog,numthresh);
% mask=imquantize(normlog,thresh);
% RGB = label2rgb(mask);
% figure;
% imshow(RGB)
mask=im2bw(normlog,thresh(end));
if sum(mask(:))/length(image(:))<toosmall
    mask = im2bw(normlog,thresh(end-1));
end
mask=imfill(mask,'holes');
end