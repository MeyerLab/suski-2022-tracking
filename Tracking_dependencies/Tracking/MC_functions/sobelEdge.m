function nuc_mask = sobelEdge(nuc_raw,nucr,debrisarea)
%nuc_filt = wiener2(nuc_raw, [5 5]);
nuc_mask = edge(nuc_raw,'Sobel');
nuc_mask = imdilate(nuc_mask,strel('disk',3));
nuc_mask=imfill(nuc_mask,'holes');
nuc_mask = bwareaopen(nuc_mask, debrisarea);
nuc_mask = imopen(nuc_mask,strel('disk',5));
nuc_mask = imerode(nuc_mask,strel('diamond',1));




%%% debugging %%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
extractmask=bwmorph(nuc_mask,'remove');
tempframe=imadjust(mat2gray(nuc_raw));
tempframe(:,:,2)=extractmask;
tempframe(:,:,3)=0;
figure,imshow(tempframe);
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%