function S = loadData_IF(conditions, pth)

%% Combine wells%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
allnames = conditions(:,1);
[~,uidx] = unique(allnames,'first');
uniquenames = allnames(sort(uidx));
uniquecondnum = numel(uniquenames);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
condnum = size(conditions,1);
load([pth, 'settings_IF.mat'],'header');

IFnames = header(2,:);
for i = 1:uniquecondnum
    condrow = find(ismember(conditions(:,1), uniquenames{i}));
    tracedata = [];
    tracestats = [];
    motherstats = [];
    IFdata = [];
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
                    if exist([pth,'IFdata_', shot, '.mat']);
                        temp = load([pth,'IFdata_', shot, '.mat']);
                        shotmat = repmat({shot}, size( temp.IFdata,1), 1);
                        %IFdata(:, find(ismember(names,'1_DAPI_mean'))) = IFdata(:, find(ismember(names,'1_DAPI_mean')))./median(IFdata(:, find(ismember(names,'1_DAPI_mean'))));
                        IFdata = [IFdata; temp.IFdata];
                        shotmat_combined = [shotmat_combined; shotmat];
                        wellindextemp = ones(size( temp.IFdata,1),3);
                        wellindextemp(:,1) = wellindextemp(:,1)*row;
                        wellindextemp(:,2) = wellindextemp(:,2)*col;
                        wellindextemp(:,3) = wellindextemp(:,3)*site;
                        wellindex_combined = [wellindex_combined; wellindextemp];
                        
                    end
                end
            end
        end
    end
    IFdata(IFdata(:)<0) = 1;
    
    %%% Extract Channels %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    S(i).area = IFdata(:, find(ismember(IFnames,'nuclear area')));
    
    
    S(i).IFarea = IFdata(:, find(ismember(IFnames,'nuclear area')));
    S(i).DAPI1 = IFdata(:, find(ismember(IFnames,'1_DAPI_mean')));
    S(i).YFP1 = IFdata(:, find(ismember(IFnames,'1_YFP_mean')));
    S(i).FarRed1 = IFdata(:, find(ismember(IFnames,'1_Cy5_mean')));
    
    S(i).DAPI2 = IFdata(:, find(ismember(IFnames,'2_DAPI_mean')));
    S(i).YFP2 = IFdata(:, find(ismember(IFnames,'2_YFP_mean')));
    S(i).FarRed2 = IFdata(:, find(ismember(IFnames,'2_Cy5_mean')));
 
    S(i).CFP3 = IFdata(:, find(ismember(IFnames,'3_CFP_mean')));
    S(i).YFP3 = IFdata(:, find(ismember(IFnames,'3_YFP_mean')));
    S(i).RFP3 = IFdata(:, find(ismember(IFnames,'3_RFP_mean')));
    S(i).YFP3cyt = IFdata(:, find(ismember(IFnames,'3_YFP_cyto ring')));
    S(i).RFP3cyt = IFdata(:, find(ismember(IFnames,'3_RFP_cyto ring')));
    S(i).iRFP3 = IFdata(:, find(ismember(IFnames,'3_iRFP_mean')));
    
    S(i).shot = shotmat_combined(:);
    S(i).pos = IFdata(:, 1:2);
    S(i).wellindex = wellindex_combined;
    
    
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
condnum = length(S);
for i = 1:condnum
    
    S(i).dna = (S(i).area.*S(i).DAPI1);
    ind = S(i).area > 200;
%     ind = S(i).area < 1000 & S(i).area > 250;
     S(i)=gateout_all(S(i), ind);
    
     S(i).CDK2 = S(i).YFP3cyt./S(i).YFP3;
    ind = S(i).CDK2 < 4 & S(i).CDK2 >.1;
    S(i)=gateout_all(S(i), ind);
    
    S(i).CDK4 = S(i).RFP3cyt./S(i).RFP3;
    ind = S(i).CDK4 < 4 & S(i).CDK4 >.1;
    S(i)=gateout_all(S(i), ind);
    
    ind = log2(S(i).iRFP3) > 7;
    S(i)=gateout_all(S(i), ind);
    
    S(i).normGem = S(i).CFP3 ./ S(i).iRFP3;
    
end
end