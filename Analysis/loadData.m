function [S] = loadData(conditions, dataDir)

%%% Analysis Options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
setting.motherOption = 1;    %0:no gating 1:mothers 2:no mothers
setting.daughterOption = 0;  %0:no gating 1:daughters 2:no daughters
setting.type = 'cycling';   %'release' or 'cycling'
setting.minTraceFrac = .25;
setting.nuc = 'IRFP';
setting.cdk = 'CFP';
setting.apc = 'YFP';
setting.crl = 'RFP';
setting.sig = '';
setting.PCNA = '';

setting.IFoption = 0;        %0:No IF 1:IF
setting.IFlabel = 'IF_';
setting.poiCdk = 0;
setting.poiApc = 1;
setting.poiCrl = 1;
setting.poiG2 = 1;
setting.poiPCNA = 0;

setting.saveName = 'D121_data.mat';
setting.dataDir = 'F:\Data\D-Replication Initiation\D121-live\Plot';

%% PROCESS DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Load data header
load([dataDir,'settings_live.mat'],'names');
names = names(2,:);

%%% Find conditions
allNames = conditions(:,1);
[~,uidx] = unique(allNames,'first');
uniqueNames = allNames(sort(uidx));
uniqueCondnum = numel(uniqueNames);
condNum = size(conditions,1);

%% PROCESS CONDITION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S = struct('traceData', [], 'traceStats', [], 'motherStats', [], 'wellindex', [], 'cellID', [], 'jitters', [], 'pos', [],'isMother',[],'isDaughter',[]);
for i = 1:uniqueCondnum
    IFdata = [];
    allIFdata = [];
    S(i).traceData = [];
    %%% Load data from all sites of condition
    condRow = find(ismember(conditions(:,1),uniqueNames{i}));   
    for c = condRow'
        rowMat = cell2mat(conditions(c,2));
        colMat = cell2mat(conditions(c,3));
        siteMat = cell2mat(conditions(c,4));
        for row = rowMat
            for col = colMat
                for site = siteMat
                    shot = [num2str(row),'_',num2str(col),'_',num2str(site)];
                    if exist([dataDir,'traceData_',shot,'.mat']) && (~setting.IFoption || exist([dataDir,setting.IFlabel,shot,'.mat']))
                        [traceDatatemp,traceStatstemp,motherStatstemp,IFdatatemp,jitters,samplecellsID, ismothertemp, allIFtemp] = ...
                            gathertracedata_NR(dataDir,shot,setting.motherOption,setting.daughterOption,setting.IFoption,setting.IFlabel);
                       try
                        S(i).traceData = [S(i).traceData; traceDatatemp];
                       catch
                           keyboard;
                       end
                        S(i).traceStats = [S(i).traceStats; traceStatstemp];
                        S(i).motherStats = [S(i).motherStats; motherStatstemp];
                        S(i).isMother = logical([S(i).isMother; ismothertemp]);
                        S(i).isDaughter = logical([S(i).isDaughter; ~isnan(traceStatstemp(:,4))]);

                        S(i).cellID = [S(i).cellID; samplecellsID];
                        wellindexTemp = ones(size(traceDatatemp,1),3);
                        wellindexTemp(:,1) = wellindexTemp(:,1)*row;wellindexTemp(:,2) = wellindexTemp(:,2)*col;wellindexTemp(:,3) = wellindexTemp(:,3)*site;
                        S(i).wellindex = [S(i).wellindex; wellindexTemp];
                        jittersTemp = repmat(reshape(jitters, [1 size(jitters,1) 2]),size(traceDatatemp,1),1);
                        S(i).jitters = [S(i).jitters; jittersTemp];
                        S(i).pos = [S(i).pos; traceDatatemp(:,:,1:2)];
                        IFdata = [IFdata; IFdatatemp];
                        allIFdata = [allIFdata; allIFtemp];
                    end
                end
            end
        end
    end
    %% Load IF data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if setting.IFoption
        load([dataDir, 'settings_live_IF.mat'],'header');
        IFnames = header(2,:);
        
        S(i).IFarea = IFdata(:, ismember(IFnames,'nuclear area'));
        S(i).DAPI1 = IFdata(:, ismember(IFnames,'1_DAPI_mean'));
        S(i).YFP1 = IFdata(:, ismember(IFnames,'1_YFP_mean'));
        S(i).FarRed1 = IFdata(:, ismember(IFnames,'1_Cy5_mean'));
        
        S(i).RFP2 = IFdata(:, ismember(IFnames,'2_RFP_mean'));

        S(i).dna = DAPI1 .* S(i).IFarea;
        S(i).x = IFdata(:, ismember(IFnames,'x'));
        S(i).y = IFdata(:, ismember(IFnames,'y'));
        
                
        S(i).jitterX = IFdata(:, ismember(IFnames,'jitterX'));
        S(i).jitterY = IFdata(:, ismember(IFnames,'jitterY'));
        S(i).fixedRow = IFdata(:, ismember(IFnames,'row'));
        S(i).fixedCol = IFdata(:, ismember(IFnames,'col'));
        S(i).fixedSite = IFdata(:, ismember(IFnames,'fixed site'));
        
    end
    %% Extract nuclear channels
    S(i).area = S(i).traceData(:,:,ismember(names,'nuclear area'));
    S(i).nucMean = S(i).traceData(:,:,ismember(names,[setting.nuc '_mean']));
    S(i).mass = S(i).area.*S(i).nucMean;
    S(i).massNorm = S(i).mass./repmat(max(S(i).mass,[],2),1,size(S(i).area,2));
    S(i).POI(:, 1) = S(i).traceStats(:,1);
    
    %% Gate on length
    numFrames = size(S(i).traceData,2);
    minLengthTrace = ceil(numFrames*setting.minTraceFrac);

    S(i).traceStats(:,5) = findFirstInMat(~isnan(S(i).area));
    S(i).traceStats(:,6) = findLastInMat(~isnan(S(i).area));
    badlengths = S(i).traceStats(:,6) - S(i).traceStats(:,5) < minLengthTrace; %| sensor(i).motherStats(:,3)<5;
    S(i)=gateout_all(S(i),~badlengths);
    
    %% Gate nuclear
    noiseThresh = .07; %(.07)
    badNoise = gatenoisy(S(i).massNorm, S(i).traceStats, setting.daughterOption, strcmp(setting.type,'release'), noiseThresh, 4,4);
    %     clear  highNoise;
    %     %%% Gate on absolute change
    %     for n = 1:size(S(i).massNorm,1)
    %         highNoise(n,1) = any(abs(diff(S(i).massNorm(n,:))) > .25);
    %     end
    S(i) = gateout_all(S(i),~(badNoise));
    
    %% Extract and gate cdk %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(setting.cdk)
        S(i).cdkNuc = S(i).traceData(:,:,ismember(names,[setting.cdk '_mean']));
        S(i).cdkCyt = S(i).traceData(:,:,ismember(names,[setting.cdk '_cyto ring']));
%         S(i).cdkLocalBg = S(i).traceData(:,:,ismember(names, [setting.cdk '_block bg']));
        
        maxThresh = 50; %threshold above which max of each trace must be  %150
        noiseThresh = .5;%0.20; %threshold rate of DHBratio change (not absolute) above which trace is considered noisy
        smoothWindow = 5;
        [S(i).cdk,badTracesCdk] = gate_Cdk2_NR(S(i).cdkNuc,S(i).cdkCyt,maxThresh,noiseThresh,smoothWindow);
        %sensor(i).cdk = sensor(i).cdkCyt./sensor(i).cdkNuc;
        S(i) = gateout_all(S(i),~badTracesCdk);
    end
    
    %% Extract and gate apc %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(setting.apc)
        S(i).apcNuc = S(i).traceData(:,:,ismember(names,[setting.apc '_mean']));% - S(i).apcLocalBg;
        %         S(i).apcNuc = correctBlankFrame(S(i).apcNuc, S(i).shot, 20:72, -5);
        %         S(i).apcNuc = fillTraceVals(S(i).apcNuc, S(i).traceStats, 5);
        
        %%% Gate cell traces
        switch setting.type
            case 'release'
                % Quiescent released cells
                meanWindow = 1:30;
                lowThresh = 20;
                highMean = nanmean(S(i).apcNuc(:,meanWindow),2) > lowThresh;
                %             noiseThresh = 200;
                %             maskTrace = [];
                %             noiseMask = ones(1,numFrames-1);
                %             noiseMask(maskTrace)=0;
                %             clear highTracestart highNoise;
                %             for n = 1:size(S(i).apcNuc,1)
                %                 highTracestart(n,1) = S(i).apcNuc(n,S(i).traceStats(n,1)) > lowThresh;
                %                 highNoise(n,1) = any(abs(diff(S(i).apcNuc(n,:))) > noiseThresh & noiseMask);
                %             end
                %S(i).apcBadTraces = highTracestart | highMean; % Don't gate out weird apc traces
                
                % Divided cells
                minMaxThresh = [15 200];
                minMaxNoise = [-50000 50000];
                traceGem = S(i).apcNuc;
                buff = 5;
                postMitosis = 0;
                [badAPCdiv] = gateSigTraces(traceGem, S(i).traceStats, buff, postMitosis, minMaxThresh, minMaxNoise);
                S(i) = gateout_all(S(i),~((badAPCdiv & S(i).isDaughter) | (highMean & ~S(i).isDaughter )));
            
            case 'cycling'
                % Cycling
                minMaxThresh = [100 1000];
                minMaxNoise = [-50000 50000];
                traceGem = S(i).apcNuc;
                buff = 5;
                postMitosis = 0;
                [badTracesApc] = gateSigTraces(traceGem, S(i).traceStats, buff, postMitosis, minMaxThresh, minMaxNoise);
                S(i) = gateout_all(S(i),~(badTracesApc));
        end
        
        %%% Transform data
        S(i).apcArea = (S(i).apcNuc).*S(i).area;
%         S(i).apcAreaSmooth = NaN*ones(size(S(i).apcNuc));
        if strcmp(setting.type,'cycling')
            S(i).apcNorm = (S(i).apcNuc-repmat(min(S(i).apcNuc,[],2),1,size(S(i).apcNuc,2)))...
                ./repmat(max(S(i).apcNuc,[],2),1,size(S(i).apcNuc,2));
            S(i).apcNormM = NaN*ones(size(S(i).apcNuc));
            S(i).apcAreaNormM = NaN*ones(size(S(i).apcArea));
            S(i).apcNorm = (S(i).apcNuc-repmat(min(S(i).apcNuc,[],2),1,size(S(i).apcNuc,2)));
            S(i).apcNorm = S(i).apcNorm./repmat(max(S(i).apcNorm,[],2),1,size(S(i).apcNorm,2));
            S(i).apcAreaNorm = (S(i).apcArea-repmat(min(S(i).apcArea,[],2),1,size(S(i).apcArea,2)));
            S(i).apcAreaNorm =    S(i).apcAreaNorm./repmat(max(S(i).apcArea,[],2),1,size(S(i).apcArea,2));
        end
        
        badTracesAPC = false(size(S(i).apcNuc,2),1);
        for n = 1:size(S(i).apcArea,1)
%             S(i).apcAreaSmooth(n,:) = nansmooth(S(i).apcArea(n,:),3);
            if strcmp(setting.type,'cycling')
                mitFrame = S(i).POI(n,1);                
                window = [-10:10];
                if ~isnan(mitFrame) & mitFrame + window(1) > 0 & mitFrame + window(end) <= size(S(i).apcNuc,2)
                    if max(S(i).apcNuc(n,mitFrame+window)) > minMaxThresh(2) & min(S(i).apcNuc(n,mitFrame+window)) < minMaxThresh(1)
                        S(i).apcNormM(n,:) =  S(i).apcNuc(n,:) - min(S(i).apcNuc(n,mitFrame+window(1):end));
                        S(i).apcNormM(n,:) = S(i).apcNormM(n,:)/max(S(i).apcNormM(n,mitFrame+window));
                        
                        S(i).apcAreaNormM(n,:) =  S(i).apcArea(n,:) - min(S(i).apcArea(n,mitFrame+window(1):end));
                        S(i).apcAreaNormM(n,:) = S(i).apcAreaNormM(n,:)/max(S(i).apcAreaNormM(n,mitFrame+window));
                    else
                        S(i).apcNormM(n,:) =  S(i).apcNorm(n,:);
                        S(i).apcAreaNormM(n,:) =  S(i).apcAreaNorm(n,:);
                        badTracesAPC(n) = true;
                    end
                else
                    S(i).apcNormM(n,:) =  S(i).apcNorm(n,:);
                    S(i).apcAreaNormM(n,:) =  S(i).apcAreaNorm(n,:);
                   badTracesAPC(n) = true;

                end
            end
        end
        %         S(i) = gateout_all(S(i),~(badTracesAPC));
        
        %%% Find APC inactivation
        if setting.poiApc
            switch setting.type
                case 'release'
                    %%% quiescent cells
                    subset = find(~S(i).isDaughter);
                    [S(i).POI(subset,2), badTracesAPC(subset)] = findAPCInact(S(i).apcNuc(subset,:), S(i).traceStats(subset,:), ...
                        struct('postBuffer',10,'smooth',5,'cycling',0,...
                        'buff',1,'thresh',.4,'preBuffer',1,'lowThresh',20, 'increase',10*.4,'trunc',50, 'medfilt',0,'early', 30,'debug',0));
                    %%% Divided cells
                    subset = find(S(i).isDaughter);
                    [S(i).POI(subset,2), badTracesAPC] = findAPCInact(S(i).apcNuc(subset,:), S(i).traceStats(subset,:), ...
                        struct('postBuffer',10,'smooth',5,'cycling',1,...
                        'buff',3,'thresh',.5,'preBuffer',5,'lowThresh',20, 'increase',10*1,'trunc',100, 'medfilt',0,'early', 5,'debug',0));
                    
                case 'cycling'
                    [S(i).POI(:,2), badTracesAPC] = findAPCInact(S(i).apcAreaNormM, S(i).traceStats, ...
                        struct('postBuffer',10,'smooth',7,'cycling',1,...
                        'buff',0,'thresh',.001,'preBuffer',5,'lowThresh',.05, 'increase',.001 * 10,'trunc',.5, 'medfilt',0,'early', 5,'debug',0));
            end
            %S(i) = gateout_all(S(i),~(badTracesAPC));
            %S(i).apcBadTraces = S(i).apcBadTraces | badTracesAPC;
        end
        
    end
    
    %% Extract and gate crl
    if ~isempty(setting.crl)
        S(i).crlNuc = S(i).traceData(:,:,ismember(names,[setting.crl '_mean']));% - S(i).crlLocalBg;
        %         S(i).crlNuc = correctBlankFrame(S(i).crlNuc, S(i).shot, 2:30, -50);
        %         S(i).crlNuc = fillTraceVals(S(i).crlNuc, S(i).traceStats, 5);
        S(i).crlArea = (S(i).crlNuc).*S(i).area;
        S(i).crlCyto = S(i).traceData(:,:,ismember(names,[setting.crl '_cyto ring']));
        
        %%% Gate out traces
        noiseThresh = 20000;
        maskTrace = [];
        %meanWindow = 20:30;
        lowThresh = 200;
        noiseMask = ones(1,numFrames-1);
        noiseMask(maskTrace)=0;
        clear lowTracestart highNoise;
        for n = 1:size(S(i).crlNuc,1)
            lowTracestart(n,1) = S(i).crlNuc(n,S(i).traceStats(n,1)) < lowThresh;
            highNoise(n,1) = any(abs(diff(S(i).crlNuc(n,:)))>noiseThresh & noiseMask);
        end
        %lowMean = nanmean(S(i).crlNuc(:,meanWindow),2) < lowThresh;
        lowMax = max(S(i).crlNuc,[],2) < lowThresh;
        %S(i).crlExpress = ~lowMax;
        S(i) = gateout_all(S(i),~( highNoise | lowMax ));
        
        %%% Transform traces
        S(i).crlNorm = NaN*ones(size(S(i).crlNuc));
        S(i).crlAreaNorm = NaN*ones(size(S(i).crlNuc));
        %S(i).crlAreaSmooth = NaN*ones(size(S(i).crlNuc));
        %S(i).crlDiff = NaN*ones(size(S(i).crlNuc));
        if strcmp(setting.type,'cycling')
            S(i).crlNormPostM = NaN*ones(size(S(i).crlNuc));
        end
        for n = 1:size(S(i).crlArea,1)
            S(i).crlNorm(n,:) =  S(i).crlNuc(n,:)./max(S(i).crlNuc(n,1:end-1));
            S(i).crlAreaNorm(n,:) =  S(i).crlArea(n,:)./max(S(i).crlArea(n,1:end-1));
            %S(i).crlNormSmooth(n,:) = nansmooth(S(i).crlNorm(n,:), 5);
            %S(i).crlDiff(n,:) = gradient(S(i).crlAreaSmooth(n,:));
            
            if strcmp(setting.type,'cycling')
                startFrame = S(i).POI(n,1);
                if ~isnan(startFrame) & startFrame +1 <= numFrames
                    S(i).crlNormPostM(n,:) =  S(i).crlArea(n,:)./max(S(i).crlArea(n,startFrame:end-1));
                else
                    S(i).crlNormPostM(n,:) =  S(i).crlAreaNorm(n,:);
                end
            end
        end
        
        %%% Find crl inactivation
        if setting.poiCrl
             switch setting.type
                 case 'release'
                     %%% quiescent cells
                     subset = find(~S(i).isDaughter);
                     S(i).POI(subset,3:4) = findCRL4Act_v2_reverse(struct('CRLtrace',S(i).crlArea(subset,:),'traceStats', S(i).traceStats(subset,:),'nuc',S(i).crlNuc(subset,:), 'cyto', S(i).crlCyto(subset,:)), ...
                         struct('cycling',0,'smooth',7 ,'buff',1 , 'postBuffer', 5 ,'postBufferSec',1,...
                         'firstD',0,'secD',-.005,'decrease', .1, 'minSlope', -.01, 'early', 3,'falseCall', .8, 'cytoCheck', 0,'minExp',75,...
                         'trunk',.4,'debug',0));
                     %%% Divided cells                     
                     subset = find(S(i).isDaughter);
                     S(i).POI(subset,3:4) = findCRL4Act_v2_reverse(struct('CRLtrace',S(i).crlArea(subset,:),'traceStats', S(i).traceStats(subset,:),'nuc',S(i).crlNuc(subset,:), 'cyto', S(i).crlCyto(subset,:)), ...
                         struct('cycling',1,'smooth',7 ,'buff',1 , 'postBuffer', 5 ,'postBufferSec',1,...
                         'firstD',0,'secD',-.005,'decrease', .05, 'minSlope', -.01, 'early', 3,'falseCall', .8, 'cytoCheck', 0,'minExp',75,...
                         'trunk',.4,'debug',0));
                     
                 case 'cycling'
                     S(i).POI(:,3:4) = findCRL4Act_v2_reverse(struct('CRLtrace',S(i).crlArea,'traceStats', S(i).traceStats,'nuc',S(i).crlNuc, 'cyto', S(i).crlCyto), ...
                         struct('cycling',1,'smooth',7 ,'buff',1 , 'postBuffer', 3 ,'postBufferSec',1,...
                         'firstD',0,'secD',-.005,'decrease', .05, 'minSlope', -.01, 'early', 3,'falseCall', .8, 'cytoCheck', 0,'minExp',75,...
                         'trunk',.4,'debug',0));
             end
             
             S(i).crlNormAct = NaN*ones(size(S(i).crlNuc));
             S(i).crlMaxAct = NaN*ones(size(S(i).crlNuc,1),1);
             for n = 1:size(S(i).crlArea,1)
                 if ~isnan(S(i).POI(n,3))
                     %%% change if use Nuc signal
                     S(i).crlNormAct(n,:) =  S(i).crlArea(n,:)/S(i).crlArea(n,S(i).POI(n,3));
                     S(i).crlMaxAct(n) = S(i).crlNuc(n,S(i).POI(n,3));
                 else
                     S(i).crlNormAct(n,:) = S(i).crlAreaNorm(n,:);
                 end
             end
             
            if setting.poiG2
                [S(i).POI(:,5), badTracesCRL] = findCRLInact(S(i).crlNormAct, S(i).traceStats, S(i).POI(:,4), ...
                    struct('postBuffer',10,'smooth',5,...
                    'buff',10,'thresh',.0025,'preBuffer',5,'lowThresh',.1, 'increase',.005*10,'trunc',.5, 'medfilt',0,'early', 15,'debug',0));
            end
        end
        
    end
    
    %% Extract and gate sig
    if ~isempty(setting.sig)
        S(i).sigNuc = S(i).traceData(:,:,ismember(names,[setting.sig '_mean']));
        S(i).sigArea = (S(i).sigNuc).*S(i).area;
        
        % Gate out traces on noise and misexpression
        noiseThresh = 20000;
        maskTrace = [];
        %meanWindow = 70:100;
        lowThresh = 20;
        noiseMask = ones(1,numFrames-1);
        noiseMask(maskTrace)=0;
        clear lowTracestart highNoise;
        for n = 1:size(S(i).sigNuc,1)
            highNoise(n,1) = any(abs(diff(S(i).sigNuc(n,:)))>noiseThresh & noiseMask);
        end
        %lowMean = nanmean(S(i).sigNuc(:,meanWindow),2) < lowThresh;
        lowMax = max(S(i).sigNuc,[],2) < lowThresh;
        S(i) = gateout_all(S(i),~( lowMax | highNoise));
        %S(i).sigExpress = ~lowMax;
        
        % Transform data
        for n = 1:size(S(i).sigArea,1)
            ignore_norm = 5;
            S(i).sigNorm(n,:) =  S(i).sigArea(n,:)./max(S(i).sigArea(n,ignore_norm:end-1));
            S(i).sigSmooth(n,:) = nansmooth(S(i).sigArea(n,:), 5);
            S(i).sigNormSmooth(n,:) = nansmooth(S(i).sigNorm(n,:), 5);
            S(i).sigDiff(n,:) = gradient(S(i).sigSmooth(n,:));
            
            
            if ~setting.quiescent
                startFrame = S(i).POI(n,1);
                if ~isnan(startFrame) & startFrame +1 <= numFrames & ~setting.quiescent
                    S(i).sigNormPostM(n,:) =  S(i).sigArea(n,:)./max(S(i).sigArea(n,startFrame:end-1));
                elseif ~setting.quiescent
                    S(i).sigNormPostM(n,:) =  S(i).sigNorm(n,:);
                end
            end
            S(i).sigNormSmoothPostM(n,:) = nansmoothm(S(i).sigNormPostM(n,:), 7, 'sgolay');
            S(i).sigNormNoise(n,:) = S(i).sigNormPostM(n,:) - S(i).sigNormSmoothPostM(n,:);
            
        end
        
        %%% Find crl inactivation
        if setting.poiCrl
            S(i).POI(:,3) = findCRL4Act(S(i).sigNormSmoothPostM,S(i).traceStats, ...
                struct('cycling',~setting.quiescent,'low',.1,'smooth',0,'smoothSec',3, ...
                'firstD',-.0,'secD',-.001,'early', 20,'falseCall', .8,'buff',10, 'postBuffer', 3,'decrease', 3*.005,'debug',0));
            
            S(i).crlNormAct = NaN*ones(size(S(i).sigNuc));
            S(i).sigMax = NaN*ones(size(S(i).sigNuc,1),1);
            
            for n = 1:size(S(i).sigArea,1)
                if ~isnan(S(i).POI(n,3))
                    S(i).sigNormAct(n,:) =  S(i).sigNuc(n,:)/S(i).sigNuc(n,S(i).POI(n,3));
                    %                     if(log2(S(i).YFP3(n)) > 9 | log2(S(i).dna) < 23.75)
                    %                         S(i).POI(n,3) = NaN;
                    %                     end
                    S(i).sigMax(n) = S(i).sigNuc(n,S(i).POI(n,3));
                    
                else
                    S(i).sigNormAct(n,:) = S(i).sigNorm(n,:);
                end
            end
        end
        
    end
    
    
    %% Extract and gate PCNA
    if ~isempty(setting.PCNA)
        S(i).pcnaNuc = S(i).traceData(:,:,ismember(names,'PCNA mean'));
        
        % Gate out traces on noise and misexpression
        noiseThresh = 20000;
        maskTrace = [];
        meanWindow = 70:100;
        lowThresh = 200;
        noiseMask = ones(1,numFrames-1);
        noiseMask(maskTrace)=0;
        clear lowTracestart highNoise;
        for n = 1:size(S(i).pcnaNuc,1)
            highNoise(n,1) = any(abs(diff(S(i).pcnaNuc(n,:)))>noiseThresh & noiseMask);
        end
        lowMean = nanmean(S(i).pcnaNuc(:,meanWindow),2) < lowThresh;
        %lowMax = max(S(i).pcnaNuc,[],2) < lowThresh;
        S(i) = gateout_all(S(i),~( lowMean | highNoise));
        %         S(i) = gateout_all(S(i),~( highNoise));
        %         S(i).pcnaExpress = ~lowMax;
        
        % load data
        S(i).filterIntensity = S(i).traceData(:,:,ismember(names,'Filter mean'));
        S(i).filterArea = S(i).traceData(:,:,contains(names,'Filter Masked area'));
        S(i).filterMaskIntensity = S(i).filterArea.*S(i).traceData(:,:,contains(names,'Filter Masked mean'));
        
        S(i).varIntensity = S(i).traceData(:,:,ismember(names,'Variance mean'));
        S(i).varStd = S(i).traceData(:,:,ismember(names,'Variance std'));
        %         S(i).varIntensityNorm = S(i).varIntensity./repmat(nanmean(S(i).pcnaNuc(:,meanWindow),2),[1 size(S(i).pcnaNuc,2)]);
        S(i).varIntensityNorm = S(i).varIntensity-repmat(min(S(i).varIntensity,[],2),[1 size(S(i).pcnaNuc,2)]);
        S(i).varArea = S(i).traceData(:,:,contains(names,'Variance Masked area'));
        S(i).varMaskIntensity = S(i).varArea.*S(i).traceData(:,:,contains(names,'Variance Masked mean'));
        %S(i).varMaskIntensityNorm = S(i).varMaskIntensity./repmat(nanmean(S(i).pcnaNuc(:,meanWindow),2),[1 size(S(i).pcnaNuc,2)]);
        
        if setting.poiPCNA
            S(i).filterPOI = getPCNAstart_C117(S(i).filterArea(:,:,2),S(i).traceStats,10,10);
            % S(i).filterPOI = getPCNAstart_C114_intensity(S(i).varIntensityNorm(:,:),S(i).traceStats,10,20);
            
        end
    end
    
    if setting.IFoption        
        S(i).allIFarea = allIFdata(:, find(ismember(IFnames,'nuclear area')));
        S(i).allDAPI1 = allIFdata(:, find(ismember(IFnames,'1_DAPI_mean')));
        S(i).allDNA = S(i).allDAPI1 .* S(i).allIFarea;
        S(i).allYFP1 = allIFdata(:, find(ismember(IFnames,'1_YFP_mean')));
        S(i).allFarRed1 = allIFdata(:, find(ismember(IFnames,'1_Cy5_mean')));
        S(i).allRFP2 = allIFdata(:, find(ismember(IFnames,'2_RFP_mean')));             
    end
end

%% Extra transformations
for i = 1:length(S)
    

end

S = rmfield(S,{'traceData'});

save(fullfile(setting.dataDir, setting.saveName),'S','-v7.3');


end

%% Extra code
%     S(i).traceData = S(i).traceData(:,1:end-1,:);
%     S(i).traceStats(S(i).traceStats == 156) = 155;
%     S(i).traceStats(:,3) = S(i).traceStats(:,2) - S(i).traceStats(:,1) + 1;
% Median filter to get rid of single frame noise
%     for numcell = 1:size(S(i).apcNuc,1)
%         S(i).apcNuc_filt(numcell,:) = S(i).apcNuc(numcell,:);
%         ind = ~isnan(S(i).apcNuc_filt(numcell,:));
%         S(i).apcNuc_filt(numcell,ind) = medfilt1(S(i).apcNuc_filt(numcell,ind),3);
%     end

%         S(i).crlLocalBg = S(i).traceData(:,:,ismember(names, [setting.crl '_block bg']));
%         [S(i).crlLocalBg, badBg] = fillTraceVals(S(i).crlLocalBg, S(i).traceStats, 5);
%         S(i).crlNuc = correctBlankFrame(S(i).crlNuc, S(i).shot, 10:40, -10);
%         S(i).crlNuc = fillTraceVals(S(i).crlNuc, S(i).traceStats, 5);
%         S(i).crlNuc(S(i).crlNuc <=0) = .01;
%         S(i).crl4Act = calcCRL4ActivityCycling_E1092(S(i).crlArea,S(i).traceStats,S(i).POI(:,3),...
%             struct('manualk_syn',NaN,'buffer',5,'k_deg', 0, 'k_mult',10, 'smooth',5));
%         S(i).crl4Act(S(i).crl4Act <= 0) = 1e-4;
%         S(i).logCrl4Act = log(S(i).crl4Act);
