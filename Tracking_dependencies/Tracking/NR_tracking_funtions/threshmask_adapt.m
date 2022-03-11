function mask=threshmask_adapt(image,blurradius)
blur=imfilter(image,fspecial('disk',blurradius),'symmetric'); %10x:3 20x:6
normlog=mat2gray(log(blur));

histeqsig=adapthisteq(normlog);
thresh=graythresh(histeqsig);
mask=im2bw(histeqsig,thresh);
mask=imfill(mask,'holes');

% thresh1=graythresh(normlog);
% mask1=im2bw(normlog,thresh1);
% mask1=imfill(mask1,'holes');
% 
% histeqsig=adapthisteq(normlog);
% thresh2=graythresh(histeqsig);
% mask2=im2bw(histeqsig,thresh2);
% mask2=imfill(mask2,'holes');
% figure,imagesc(image)
% figure,imagesc(histeqsig)
end