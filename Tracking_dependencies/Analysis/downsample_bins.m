function [keep] = downsample_bins(data,edges, varargin)
%DOWNSAMPLE_BINS downsample data bases on bin numbers
%
%   INPUTS:
%   DATA     - values to bin on
%   EDGES     - bin edges
%   
%   Optional name-value:
%   varargin(1)   - number to downsample by, otherwise will use median
%   bin size
%
%   OUTPUTS:
%   KEEP   - index of data to keep

    
[N,~,bins]=histcounts(data,edges);
keep= [];
if nargin == 3
    undersample = varargin{1};
else
    undersample = round(median(N(N ~= 0)));
end

for j = 1:length(N)
    if N(j) > undersample
        inds = find(bins==j);
        keep = [keep; randsample(inds,undersample,false)];
    else
        keep = [keep; find(bins==j)];
    end
end
end

