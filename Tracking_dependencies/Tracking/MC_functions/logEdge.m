function nuc_mask = logEdge(nuc_raw,nucr,threshold, debrisarea)
sigma = 3;
%nuc_raw(nuc_raw < 1) = 1;
%nuc_log=log(nuc_raw);

%nuc_mask = edge(nuc_raw, 'log', .7, sigma);
 nuc_mask = edge(nuc_raw,'Sobel');
se90 = strel('line',3,90);
se45 = strel('line',3,45);
se0 = strel('line',3,0);
nuc_mask = imdilate(nuc_mask,[se90 se45 se0]);
nuc_mask=imfill(nuc_mask,'holes');
% nuc_mask=~bwmorph(~nuc_mask,'diag');
% nuc_mask=~bwmorph(~nuc_mask,'bridge');
nuc_mask = bwareaopen(nuc_mask, debrisarea);


%%% debugging %%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
extractmask=bwmorph(nuc_mask,'remove');
tempframe=imadjust(mat2gray(nuc_raw));
tempframe(:,:,2)=extractmask;
tempframe(:,:,3)=0;
figure,imshow(tempframe);
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%