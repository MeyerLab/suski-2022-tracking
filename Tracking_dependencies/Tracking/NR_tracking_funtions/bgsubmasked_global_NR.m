function real=bgsubmasked_global_NR(raw,nanmask,numblocks,compression,sampleprctile)
[orgheight,orgwidth]=size(raw);
rawcomp=imresize(raw,1/compression,'nearest');
nanmaskcomp=imresize(nanmask,1/compression,'nearest');

fg=rawcomp; fg(nanmaskcomp)=NaN;
tfg=fg; tfg(nanmaskcomp)=[];
vals=tfg(:);

globalbg=prctile(vals,sampleprctile);
if numblocks>1
    [height width]=size(fg);
    blockimg=blockpercentile_blockimage_NR(fg,numblocks,sampleprctile);
    missingx=fillmissing(blockimg,'spline',1);
    missingy=fillmissing(blockimg,'spline',2);
    blockimg=(missingx+missingy)/2;
    bg=imresize(blockimg,[height width],'bicubic');
    bg=imresize(bg,[orgheight orgwidth],'bicubic');
    %figure,imagesc(bg);
else
    bg=ones(size(raw))*globalbg;
end
real=raw-bg;