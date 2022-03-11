function [frames, perctraces] = get_perc_traces( traces, POI, thresh, greater)

    numcells=size(traces,1);
    numframes=size(traces,2);
    alignmat=NaN*ones(numcells,numframes*2);
    
    for j=1:numcells
        if ~isnan(POI(j))
            alignmat(j,(numframes-POI(j)+1):(2*numframes-POI(j)))=traces(j,:);
        end
    end
    if greater
        perctraces = sum(alignmat > thresh, 1)./sum(~isnan(alignmat), 1);
    else
        perctraces = sum(alignmat < thresh, 1)./sum(~isnan(alignmat), 1);
    end
    nanmedian(alignmat,1);
    %semtrace=nanstd(alignmat,1);%./(2*sqrt(sum(~isnan(alignmat),1)));
    frames=-numframes:numframes-1;

    
end

