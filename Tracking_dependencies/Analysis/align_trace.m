function [frames, alignmat] = align_trace( varargin)
    traces = varargin{1};
    POI = varargin{2};
    
    numcells=size(traces,1);
    numframes=size(traces,2);
    alignmat=NaN*ones(numcells,numframes*2);
    
    for j=1:numcells
        if ~isnan(POI(j))
            alignmat(j,(numframes-POI(j)+2):(2*numframes-POI(j))+1)=traces(j,:);
        end
    end
    frames=-numframes:numframes-1;

    
end

