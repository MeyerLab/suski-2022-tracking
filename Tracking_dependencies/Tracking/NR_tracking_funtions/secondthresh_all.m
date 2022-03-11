function newmask=secondthresh_all(image,blurradius,mask)
blur=imfilter(log(image),fspecial('disk',blurradius),'symmetric');

maskvals=image(mask);
normlogvals=mat2gray(log(maskvals));
higherthresh=graythresh(normlogvals);
loghigherthresh=higherthresh*range(log(maskvals))+min(log(maskvals));
higher_mask=blur>loghigherthresh;

newmask=imfill(higher_mask,'holes');
end