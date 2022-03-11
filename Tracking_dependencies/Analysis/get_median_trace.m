function [frames, meantrace, semtrace] = get_median_trace( varargin)
    traces = varargin{1};
    POI = varargin{2};
    if nargin == 3
        sem = varargin{3};
    else
        sem = 0;
    end
    

    numcells=size(traces,1);
    numframes=size(traces,2);
    alignmat=NaN*ones(numcells,numframes*2);
    
    for j=1:numcells
        if ~isnan(POI(j))
            alignmat(j,(numframes-POI(j)+2):(2*numframes-POI(j))+1)=traces(j,:);
        end
    end
    meantrace=nanmedian(alignmat,1);
    if sem
        semtrace=nanstd(alignmat,1)./(2*sqrt(sum(~isnan(alignmat),1)));
    else
        semtrace=nanstd(alignmat,1);%./(2*sqrt(sum(~isnan(alignmat),1)));
    end
    frames=-numframes:numframes-1;

    
end

