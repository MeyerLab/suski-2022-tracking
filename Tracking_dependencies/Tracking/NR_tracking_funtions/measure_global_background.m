function [globalbg] = measure_global_background(raw,nanmask,compression,sampleprctile,option)
rawcomp=imresize(raw,1/compression,'nearest');
nanmaskcomp=imresize(nanmask,1/compression,'nearest');

fg=rawcomp; fg(nanmaskcomp)=NaN;
tfg=fg; tfg(nanmaskcomp)=[];
vals=tfg(:);
switch option
    case 'mode'
        binmax=prctile(vals,95);
        binmin=prctile(vals,5);
        vals=vals(vals<binmax & vals>binmin);
        [kval,xval]=ksdensity(vals);
        globalbg=xval(find(kval==max(kval),1,'first'));
    case 'percentile'       
        globalbg=prctile(vals,sampleprctile);
    otherwise
        error('global background option not available')
end
