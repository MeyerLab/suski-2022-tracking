function blockimg=blockpercentile_blockimage_NR_2(initimg,blocknum,prctilethresh,pixelthresh)
fun=@(block_struct) calcperc(block_struct.data,prctilethresh,pixelthresh);
[height,width]=size(initimg);
blockheight=ceil(height/blocknum);
blockwidth=ceil(width/blocknum);
blockimg=blockproc(initimg,[blockheight blockwidth],fun);
end

function perc=calcperc(img,thresh,pixelthresh)
%imagesc(img)
if sum(~isnan(img(:)))<pixelthresh
    perc=NaN;
else
    perc=prctile(img(:),thresh);
end
end