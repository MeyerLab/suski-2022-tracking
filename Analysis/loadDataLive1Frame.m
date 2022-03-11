function sensor = loadDataLive1Frame(conditions, pth)

%% Combine wells%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
allnames = conditions(:,1);
[~,uidx] = unique(allnames,'first');
uniquenames = allnames(sort(uidx));
uniquecondnum = numel(uniquenames);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
condnum = size(conditions,1);
load([pth, 'settings_live.mat'],'names');
liveheader = names;
namesLive = liveheader(2,:);

load([pth, 'settings_live_IF.mat'],'header');
IFheader = header;
namesIF = IFheader(2,:);
for ii = 1:uniquecondnum
    condrow = find(ismember(conditions(:,1), uniquenames{ii}));
    tracedata = [];
    tracestats = [];
    motherstats = [];
    IFdata_combined = [];
    tracedata_combined = [];
    wellindex_combined = [];
    shotmat_combined = [];
    cc = 0;
    for c = condrow'
        rowmat = cell2mat(conditions(c,2));
        colmat = cell2mat(conditions(c,3));
        sitemat = cell2mat(conditions(c,4));
        for row = rowmat
            for col = colmat
                for site = sitemat
                    cc = cc+1;
                    shot = [num2str(row),'_',num2str(col), '_', num2str(site)];
                    if exist([pth,'IF_', shot, '.mat']) & exist([pth,'tracedata_', shot, '.mat']);
                        load([pth,'IF_', shot, '.mat'], 'IFdata');
                        shotmat = repmat({shot}, size(IFdata,1), 1);
                        %IFdata(:, find(ismember(names,'1_DAPI_mean'))) = IFdata(:, find(ismember(names,'1_DAPI_mean')))./median(IFdata(:, find(ismember(names,'1_DAPI_mean'))));
                        IFdata_combined = [IFdata_combined; IFdata];
                        shotmat_combined = [shotmat_combined; shotmat];
                        wellindextemp = ones(size(IFdata,1),3);
                        wellindextemp(:,1) = wellindextemp(:,1)*row;
                        wellindextemp(:,2) = wellindextemp(:,2)*col;
                        wellindextemp(:,3) = wellindextemp(:,3)*site;
                        wellindex_combined = [wellindex_combined; wellindextemp];
                        load([pth,'tracedata_', shot, '.mat'], 'tracedata');
                        tracedata_combined = [tracedata_combined; squeeze(tracedata)];
                        
                    end
                end
            end
        end
    end
    IFdata_combined(IFdata_combined(:)<0) = 1;
    
    %%% Extract Channels %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sensor(ii).area = IFdata_combined(:, find(ismember(namesIF,'nuclear area')));
    sensor(ii).DAPI1 = IFdata_combined(:, find(ismember(namesIF,'1_DAPI_mean')));
    sensor(ii).FarRed1 = IFdata_combined(:, find(ismember(namesIF,'1_Cy5_mean')));
  

    sensor(ii).YFPlive = tracedata_combined(:,find(ismember(namesLive,'YFP_mean')));
    sensor(ii).RFPlive = tracedata_combined(:,find(ismember(namesLive,'RFP_mean')));
    
    sensor(ii).shot = shotmat_combined(:);
    sensor(ii).pos = IFdata_combined(:, 1:2);
    sensor(ii).wellindex = wellindex_combined;
    
    
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
condnum = length(sensor);
for i = 1:condnum
    
    sensor(i).dna = (sensor(i).area.*sensor(i).DAPI1);
    ind = sensor(i).area > 1000;
    sensor(i)=gateout_all(sensor(i), ind);
    
    
    
end
end