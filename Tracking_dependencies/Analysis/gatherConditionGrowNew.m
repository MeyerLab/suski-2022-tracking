function [S] = gatherConditionGrowNew(settings)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%%% Find conditions
allNames = settings.conditions(:,1);
[~,uidx] = unique(allNames,'first');
uniqueNames = allNames(sort(uidx));
uniqueCondnum = numel(uniqueNames);

%% PROCESS CONDITION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S = struct();
for i = 1:uniqueCondnum
    cellStore = {};
    S(i).traceData = [];
    S(i).shot = [];
    idx = 1;
    %%% Load data from all sites of condition
    condRow = find(ismember(settings.conditions(:,1),uniqueNames{i}));
    for c = condRow'
        rowMat = cell2mat(settings.conditions(c,2));
        colMat = cell2mat(settings.conditions(c,3));
        siteMat = cell2mat(settings.conditions(c,4));
        for row = rowMat
            for col = colMat
                for site = siteMat
                    shot = [num2str(row),'_',num2str(col),'_',num2str(site)];
                    if exist(fullfile(settings.dataDir,[settings.liveLabel,shot,'.mat'])) && (~settings.IFoption || exist(fullfile(settings.dataDir,[settings.IFlabel,shot,'.mat'])))
                            loadLive = load(fullfile(settings.dataDir,['traceData_',shot,'.mat']),'tracedata','genealogy','jitters','bg');
                            S(i).traceData = [S(i).traceData; single(loadLive.tracedata)];
                            cellStore{idx,2} = loadLive.genealogy;
                            cellStore{idx,3} = uint16(repmat(reshape(loadLive.jitters, [1 size(loadLive.jitters,1) 2]),size(loadLive.tracedata,1),1));
                            cellStore{idx,4} = ones(size(loadLive.tracedata,1),3).*[row col site];
                            if settings.IFoption
                                loadIF = load(fullfile(settings.dataDir,[settings.IFlabel,shot,'.mat']),'IFdata');
                                cellStore{idx,5} = loadIF.IFdata;
                            end
                            
                            cellStore{idx,6} = repmat(reshape(loadLive.bg', [1 size(loadLive.bg,2)  size(loadLive.bg,1)]),size(loadLive.tracedata,1),1);
                            S(i).shot = [S(i).shot; repmat(string(shot),size(loadLive.genealogy,1),1)];
                            idx = idx + 1;
                    end
                end
            end
        end
    end
    %%% Add count to genealogy account for array concatenation
    cellCount = cellfun(@length,cellStore(:,2));
    cellCount = [0; cellCount(1:end-1)];
    cellCountAccum = cumsum(cellCount);
    cellIDoffsets = num2cell(cellCountAccum);
    cellStore(:,2) = cellfun(@(x,y) x + y, cellStore(:,2),cellIDoffsets,'UniformOutput',false);
    
    %%% Concatenate arrays
%     S(i).traceData = cell2mat(cellStore(:,1));    
    S(i).cellID = (1:size(S(i).traceData,1))';
    S(i).genealogy = cell2mat(cellStore(:,2));
    S(i).wellindex = cell2mat(cellStore(:,4));
    S(i).jitters =  cell2mat(cellStore(:,3));
    S(i).bg =  cell2mat(cellStore(:,6));
%     S(i).shot = string(strcat(num2str(S(i).wellindex(:,1)), '_',num2str(S(i).wellindex(:,2)),'_',num2str(S(i).wellindex(:,3))));
    S(i).pos = S(i).traceData(:,:,1:2);
    if settings.IFoption
        S(i).IFdata = cell2mat(cellStore(:,5));
    end
    clear cellStore loadLive
    
    
    %%% Get trace information
    mothers=S(i).genealogy(~isnan(S(i).genealogy));
    S(i).isMother = false(size(S(i).genealogy));
    S(i).isMother(unique(mothers)) = true;
    S(i).isDaughter = ~isnan(S(i).genealogy);
    S(i).traceStats = getstats(S(i).traceData,S(i).genealogy); %start,end,length,genealogy
    S(i).motherStats = getmotherstatsonly(S(i).traceData,S(i).traceStats,find(S(i).isDaughter));
    
    
    %%% Copy mother traces to daugher cells
    S(i).traceData=linkancestry_NR_full(S(i).traceData,S(i).traceStats,find(S(i).isDaughter));

end

end