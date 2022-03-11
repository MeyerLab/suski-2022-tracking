function [binData, edges] = bin_xy(xdata, ydata, functions, varargin)
%BIN_XY bin data based on x data
%
%   INPUTS:
%   XDATA     - values to bin on
%   YDATA     - matrix of Y values to bin, rows are cells, cols are
%               measurement
%   FUNCTIONS - functions to apply to each bin
%               Options: 'mean', 'median', 'sem', 'std', 'perc'
%   Optional name-value:
%   'Threshs'   - cell array of cell arrays of threshs, including polarity
%                 (-1 for < thresh, 1 for > thresh) 
%                  e.g {{threshs, polarity},{threshs, polarity}...}
%   'Edges'     - Edges for binning, if not included use default bins from
%                 histcounts()
%   'Groups'    - Grouping variable
%   'MinGroupBin' - 

%   OUTPUTS:
%   BINDATA   - structure of bin values (array of each measurement)
%   EDGES     - bin edges

checkLength = @(x) size(x,2) == size(ydata,2);
p = inputParser;
addParameter(p, 'Threshs',[]);
addParameter(p, 'Edges',[]);
addParameter(p, 'MinGroupBin',1);
addParameter(p, 'Groups',[],checkLength);
addParameter(p, 'GroupFunc',@nanmedian);
parse(p,varargin{:})
GroupFunc = p.Results.GroupFunc;

ind = ~any(isnan(ydata),2) ;
xdata = xdata(ind);
ydata = ydata(ind,:);

if isempty(p.Results.Edges)
    [~,edges,bins]=histcounts(xdata);
else
    [~,edges,bins]=histcounts(xdata,p.Results.Edges);
end


        
for m = 1:size(ydata, 2) 
    xInd = [];
    groupData = [];
    
    for j = 1:length(edges)-1
        index = bins == j;
        binY = ydata(index, m);       
        if ~isempty(p.Results.Groups)
            binGroups = p.Results.Groups(index);
            unique_Groups = unique(binGroups);

            for g = 1:length(unique_Groups)
                g_filter = unique_Groups(g);
                inds = (contains(binGroups,g_filter));
                data = binY(inds);
                if length(data) >= p.Results.MinGroupBin
                    temp = GroupFunc(data);
                    groupData = [groupData temp];
                    xInd = [xInd repmat(j,1 ,length(temp))];
                end                
            end
        end
        
        binData(m).numCells(j) = sum(index);
        for f = 1:length(functions)
            switch functions{f}
                case 'mean'
                    binData(m).mean(j) = nanmean(binY);
                    binData(m).sem(j) = nanstd(binY)/sqrt(length(binY));
                    if ~isempty(groupData)
                        binData(m).g_mean(j) = nanmean(groupData(xInd == j));
                        binData(m).g_sem(j)= nanstd(groupData(xInd == j))/sqrt(length(groupData(xInd == j)));
                    end
                case 'median'
                    binData(m).median(j) = nanmedian(binY);
                    if ~isempty(groupData)
                        binData(m).g_median(j) = nanmedian(groupData(xInd == j));
                    end
                case 'std'
                    binData(m).std(j) = nanstd(binY);
                    if ~isempty(groupData)
                        binData(m).g_median(j) = nanstd(groupData(xInd == j));
                    end
                case 'perc'
                    Threshs =  p.Results.Threshs(m);
                    thresh = Threshs{1}{1};
                    polarity = Threshs{1}{2};
                    for level = 1:length(thresh)
                        binData(m).perc(level,j) = sum(polarity * binY > polarity * thresh(level))/length(binY);
                        binData(m).percerr(level,j) = sqrt(binData(m).perc(level,j) * (1-binData(m).perc(level,j))/length(binY));
                    end
            end
            
            
        end
    end
    binData(m).xInd = xInd;
    binData(m).g_ydata = groupData;
end

