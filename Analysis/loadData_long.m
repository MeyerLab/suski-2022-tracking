%% DEFINE PARAMETERS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
settings.conditions = {
        'siCtrl NR75 DMSO/DMSO',1:2,2:3,1,[0 0 0],[0 10]; %1    
%         'siFZR1 NR75 DMSO/DMSO',1:2,4:5,1,[0 0 0],[0 10]; %1    
%         'siCtrl NR75 DMSO/DOX',3:4,2:3,1,[0 0 0],[0 10]; %1    
%         'siFZR1 NR75 DMSO/DOX',3:4,4:5,1,[0 0 0],[0 10]; %1    
%         'siCtrl NR75 46i/DMSO',5:6,2:3,1,[0 0 0],[0 10]; %1    
%         'siFZR1 NR75 46i/DMSO',5:6,4:5,7,1,[0 0 0],[0 10]; %1   
%         'siCtrl NR75 46i/DOX',7:8,2:3,1,[0 0 0],[0 10]; %1    
%         'siFZR1 NR75 46i/DOX',7:8,4:5,1,[0 0 0],[0 10]; %1   
% 
%         
%         'siCtrl NR77 DMSO/DMSO',1:2,6:7,1,[0 0 0],[0 10]; %1
%         'siFZR1 NR77 DMSO/DMSO',1:2,8:9,1,[0 0 0],[0 10]; %1
%         'siCtrl NR77 DMSO/DOX',3:4,6:7,1,[0 0 0],[0 10]; %1
%         'siFZR1 NR77 DMSO/DOX',3:4,8:9,1,[0 0 0],[0 10]; %1
%         'siCtrl NR77 46i/DMSO',5:6,6:7,1,[0 0 0],[0 10]; %1
%         'siFZR1 NR77 46i/DMSO',5:6,8:9,7,1,[0 0 0],[0 10]; %1
%         'siCtrl NR77 46i/DOX',7:8,6:7,1,[0 0 0],[0 10]; %1
%         'siFZR1 NR77 46i/DOX',7:8,8:9,1,[0 0 0],[0 10]; %1
%         
%         'siCtrl NR78 DMSO/DMSO',1:2,10:11,1,[0 0 0],[0 10]; %1
%         'siCtrl NR78 DMSO/DOX',3:4,10:11,1,[0 0 0],[0 10]; %1
%         'siCtrl NR78 46i/DMSO',5:6,10:11,1,[0 0 0],[0 10]; %1
%         'siCtrl NR78 46i/DOX',7:8,10:11,1,[0 0 0],[0 10]; %1
        };

%%% Paths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
settings.dataDir = 'F:\Data\F-G1S\F023-live\Data';
settings.liveLabel = 'traceData_';
settings.liveNames = fullfile(settings.dataDir,'settings_live.mat');
settings.IFNames = fullfile(settings.dataDir, 'settings_live_IF.mat');
settings.saveName = 'F023_data.mat';
settings.saveDir = 'F:\Data\F-G1S\F023-live\Plot';
%%% Analysis Options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
settings.daughterGate = 0;  %0:no gating 1:daughters 2:no daughters
settings.type = 'release';   %'release', 'cycling' or 'hybrid
settings.minTraceFrac = .25;
settings.nuc = 'iRFP';
settings.cdk = 'CFP';
settings.apc = 'YFP';
settings.crl = 'RFP';
settings.sig = '';
settings.PCNA = '';

%%% POI analysis
settings.poiCdk = 0;
settings.poiApc = 1;
settings.poiCrl = 1;
settings.poiG2 = 1;
settings.poiPCNA = 0;

%%% Live-fixed analysis
settings.IFoption = 0;        %0:No IF 1:IF
settings.IFlabel = 'IF_';


%% PROCESS DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Load data header
load(settings.liveNames,'names');
names = names(2,:);
if settings.IFoption
    load(settings.IFNames,'header');
    IFnames = header(2,:);
end

%%% Load conditions
S = gatherConditionGrowNew(settings);
allIF = struct();
%% PROCESS CONDITION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(S)
    %%% Load IF data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if settings.IFoption
        S(i).IFarea = S(i).IFdata(:, ismember(IFnames,'nuclear area'));
        S(i).DAPI1 = S(i).IFdata(:, ismember(IFnames,'1_DAPI_mean'));
        S(i).FarRed1 = S(i).IFdata(:, ismember(IFnames,'1_FarRed_mean'));
        S(i).FarRed1cyt = S(i).IFdata(:, ismember(IFnames,'1_FarRed_cyto ring'));
        
        S(i).dna = S(i).DAPI1 .* S(i).IFarea;
        S(i).x = S(i).IFdata(:, ismember(IFnames,'x'));
        S(i).y = S(i).IFdata(:, ismember(IFnames,'y'));
        
        S(i).jitterX = S(i).IFdata(:, ismember(IFnames,'jitterX'));
        S(i).jitterY = S(i).IFdata(:, ismember(IFnames,'jitterY'));
        S(i).fixedRow = S(i).IFdata(:, ismember(IFnames,'row'));
        S(i).fixedCol = S(i).IFdata(:, ismember(IFnames,'col'));
        S(i).fixedSite = S(i).IFdata(:, ismember(IFnames,'fixed site'));
        
        
        %%% Ungated
        allIF(i).IFarea = S(i).IFdata(:, find(ismember(IFnames,'nuclear area')));
        allIF(i).DAPI1 = S(i).IFdata(:, find(ismember(IFnames,'1_DAPI_mean')));
        allIF(i).DNA = allIF(i).DAPI1 .* allIF(i).IFarea;      
        allIF(i).FarRed1 = S(i).IFdata(:, find(ismember(IFnames,'1_FarRed_mean')));
        allIF(i).FarRed1cyt = S(i).IFdata(:, find(ismember(IFnames,'1_FarRed_cyto ring')));
    end
    %% Gate based on division %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if settings.daughterGate == 1
        S(i) = gateout_all(S(i),S(i).isDaughter);
    elseif settings.daughterGate == 2
        S(i) = gateout_all(S(i),~S(i).isDaughter);
    end
    %% Extract nuclear channels%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    S(i).area = S(i).traceData(:,:,ismember(names,'nuclear area'));
    S(i).nucMean = S(i).traceData(:,:,ismember(names,[settings.nuc '_mean']));
    S(i).mass = S(i).area.*S(i).nucMean;
    S(i).massNorm = S(i).mass./repmat(max(S(i).mass,[],2),1,size(S(i).area,2));
    S(i).POI(:, 1) = S(i).traceStats(:,1);
    
    %% Gate on length%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    numFrames = size(S(i).traceData,2);
    minLengthTrace = ceil(numFrames*settings.minTraceFrac);
    
    S(i).traceStats(:,5) = findFirstInMat(~isnan(S(i).area));
    S(i).traceStats(:,6) = findLastInMat(~isnan(S(i).area));
    badlengths = (S(i).traceStats(:,6) - S(i).traceStats(:,5) < minLengthTrace) & ~S(i).isMother; %| sensor(i).motherStats(:,3)<5;
    S(i)=gateout_all(S(i),~badlengths);
    
    %% Gate nuclear %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    noiseThresh = .07; %(.07)
%     badNoise = gatenoisy(S(i).massNorm, S(i).traceStats, settings.daughterOption, strcmp(settings.type,'release'), noiseThresh, 4,4);
    badNoise = gatenoisy_hybrid(S(i).massNorm, S(i).traceStats, S(i).isDaughter, noiseThresh, 4,4);
    %     clear  highNoise;
    %     %%% Gate on absolute change
    %     for n = 1:size(S(i).massNorm,1)
    %         highNoise(n,1) = any(abs(diff(S(i).massNorm(n,:))) > .25);
    %     end
    S(i) = gateout_all(S(i),~(badNoise));
    
    %% Extract and gate cdk %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(settings.cdk)
        bg_ind = find(ismember(names,[settings.cdk '_mean']))- 4;
        bg_median = median(S(i).bg(:,:,bg_ind),1);
        bg_median = bg_median + 550;
        S(i).cdkNuc = S(i).traceData(:,:,ismember(names,[settings.cdk '_mean'])) + S(i).bg(:,:,bg_ind) - bg_median;
        S(i).cdkCyt = S(i).traceData(:,:,ismember(names,[settings.cdk '_cyto ring'])) + S(i).bg(:,:,bg_ind) - bg_median;
        %         S(i).cdkLocalBg = S(i).traceData(:,:,ismember(names, [setting.cdk '_block bg']));
        
        maxThresh = 200; %threshold above which max of each trace must be  %150
        noiseThresh = 100;%0.20; %threshold rate of DHBratio change (not absolute) above which trace is considered noisy
        smoothWindow = 5;
        [S(i).cdk,badTracesCdk] = gate_Cdk2_NR(S(i).cdkNuc,S(i).cdkCyt,maxThresh,noiseThresh,smoothWindow);
        lowCDK = min(S(i).cdkNuc,[],2) < 10;
        S(i).goodCdk = ~badTracesCdk & ~lowCDK;
%         S(i) = gateout_all(S(i),~badTracesCdk & ~lowCDK & S(i).isDaughter);
    end
    
    %% Extract and gate apc %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(settings.apc)
        S(i).apcNuc = S(i).traceData(:,:,ismember(names,[settings.apc '_mean']));% - S(i).apcLocalBg;
%         S(i).apcNuc = correctBlankFrame(S(i).apcNuc, S(i).shot, 3:10, -5);
%         S(i).apcNuc = fillTraceVals(S(i).apcNuc, S(i).traceStats, 5);
        
        %%% Gate cell traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        highMean = false(size(S(i).apcNuc,1),1);
        badAPCdiv = false(size(S(i).apcNuc,1),1);        
        if strcmp(settings.type, 'release') | strcmp(settings.type, 'hybrid')
            % Quiescent released cells
            meanWindow = 1:20;
            lowThresh = 2000;
            highMean(~S(i).isDaughter) = nanmean(S(i).apcNuc(~S(i).isDaughter,meanWindow),2) > lowThresh;
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
        end
        if strcmp(settings.type, 'cycling') | strcmp(settings.type, 'hybrid')
            % Divided cells
            minMaxThresh = [100 1000];
            minMaxNoise = [-50000 50000];
            traceGem = S(i).apcNuc;
            buff = 5;
            postMitosis = 0;
            badAPCdiv(S(i).isDaughter) = gateSigTraces(traceGem(S(i).isDaughter,:), S(i).traceStats(S(i).isDaughter,:),...
                buff, postMitosis, minMaxThresh, minMaxNoise);          
        end
        S(i) = gateout_all(S(i),~(badAPCdiv | highMean));
        
        
        %%% Transform data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        S(i).apcArea = (S(i).apcNuc).*S(i).area;
        %         S(i).apcAreaSmooth = NaN*ones(size(S(i).apcNuc));
                badTracesAPC = false(size(S(i).apcNuc,2),1);

        if strcmp(settings.type,'cycling') | strcmp(settings.type,'hybrid') 
            S(i).apcNorm = (S(i).apcNuc-repmat(min(S(i).apcNuc,[],2),1,size(S(i).apcNuc,2)))...
                ./repmat(max(S(i).apcNuc,[],2),1,size(S(i).apcNuc,2));
            S(i).apcNormM = NaN*ones(size(S(i).apcNuc));
            S(i).apcAreaNormM = NaN*ones(size(S(i).apcArea));
            S(i).apcNorm = (S(i).apcNuc-repmat(min(S(i).apcNuc,[],2),1,size(S(i).apcNuc,2)))...
                ./repmat(max(S(i).apcNuc,[],2),1,size(S(i).apcNuc,2));
            S(i).apcAreaNorm = (S(i).apcArea-repmat(min(S(i).apcArea,[],2),1,size(S(i).apcArea,2)))...
                ./repmat(max(S(i).apcArea,[],2),1,size(S(i).apcArea,2));
        end
        
        for n = 1:size(S(i).apcArea,1)
            %             S(i).apcAreaSmooth(n,:) = nansmooth(S(i).apcArea(n,:),3);
            if strcmp(settings.type,'cycling') | strcmp(settings.type,'hybrid') 
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
        
        %%% Find APC inactivation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if settings.poiApc            
            %%% quiescent cells
            subset = find(~S(i).isDaughter);
            [S(i).POI(subset,2), badTracesAPC(subset)] = findAPCInact(S(i).apcNuc(subset,:), S(i).traceStats(subset,:), ...
                struct('postBuffer',10,'smooth',5,'cycling',0,...
                'buff',1,'thresh',10,'preBuffer',1,'lowThresh',200, 'increase',10*10,'trunc',1000, ...
                'medfilt',0,'early', 0,'debug',0));
            %%% Divided cells
            if strcmp(settings.type,'cycling') | strcmp(settings.type,'hybrid')
                subset = find(S(i).isDaughter);
                traces = S(i).apcNuc;
                [S(i).POI(subset,2), badTracesAPC(subset)] = findAPCInact(S(i).apcAreaNormM(subset,:), S(i).traceStats(subset,:), ...
                    struct('postBuffer',10,'smooth',7,'cycling',1,...
                    'buff',0,'thresh',.002,'preBuffer',5,'lowThresh',.05, 'increase',.002 * 10,...
                    'trunc',.5, 'medfilt',0,'early', 3,'debug',0));
            end
            
            %S(i) = gateout_all(S(i),~(badTracesAPC));
            %S(i).apcBadTraces = S(i).apcBadTraces | badTracesAPC;
            
        end
        
    end
    
    %% Extract and gate crl %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(settings.crl)
        S(i).crlNuc = S(i).traceData(:,:,ismember(names,[settings.crl '_mean']));% - S(i).crlLocalBg;
%         S(i).crlNuc = correctBlankFrame(S(i).crlNuc, S(i).shot, 3:10, -5);
%         S(i).crlNuc = fillTraceVals(S(i).crlNuc, S(i).traceStats, 5);
        S(i).crlArea = (S(i).crlNuc).*S(i).area;
        S(i).crlCyto = S(i).traceData(:,:,ismember(names,[settings.crl '_cyto ring']));
        
        %%% Gate out traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        noiseThresh = 20000;
        maskTrace = [];
        %meanWindow = 20:30;
        lowThresh = 500;
        noiseMask = ones(1,numFrames-1);
        noiseMask(maskTrace)=0;
        clear lowTracestart highNoise;
        for n = 1:size(S(i).crlNuc,1)
            lowTracestart(n,1) = S(i).crlNuc(n,S(i).traceStats(n,1)) < lowThresh;
            highNoise(n,1) = any(abs(diff(S(i).crlNuc(n,:)))>noiseThresh & noiseMask);
%             inds = find(highNoise);
        end
        %lowMean = nanmean(S(i).crlNuc(:,meanWindow),2) < lowThresh;
        lowMax = max(S(i).crlNuc,[],2) < lowThresh;
        %S(i).crlExpress = ~lowMax;
        S(i) = gateout_all(S(i),~( highNoise | lowMax ));
        
        %%% Transform traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        S(i).crlNorm = NaN*ones(size(S(i).crlNuc));
        S(i).crlAreaNorm = NaN*ones(size(S(i).crlNuc));
        %S(i).crlAreaSmooth = NaN*ones(size(S(i).crlNuc));
        %S(i).crlDiff = NaN*ones(size(S(i).crlNuc));
        if strcmp(settings.type,'cycling')
            S(i).crlNormPostM = NaN*ones(size(S(i).crlNuc));
        end
        for n = 1:size(S(i).crlArea,1)
            S(i).crlNorm(n,:) =  S(i).crlNuc(n,:)./max(S(i).crlNuc(n,1:end-1));
            S(i).crlAreaNorm(n,:) =  S(i).crlArea(n,:)./max(S(i).crlArea(n,1:end-1));
            %S(i).crlNormSmooth(n,:) = nansmooth(S(i).crlNorm(n,:), 5);
            %S(i).crlDiff(n,:) = gradient(S(i).crlAreaSmooth(n,:));
            
            if strcmp(settings.type,'cycling')
                startFrame = S(i).POI(n,1);
                if ~isnan(startFrame) & startFrame +1 <= numFrames
                    S(i).crlNormPostM(n,:) =  S(i).crlArea(n,:)./max(S(i).crlArea(n,startFrame:end-1));
                else
                    S(i).crlNormPostM(n,:) =  S(i).crlAreaNorm(n,:);
                end
            end
        end
        
        %%% Find crl inactivation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if settings.poiCrl          
            %%% quiescent cells
            subset = find(~S(i).isDaughter);
            S(i).POI(subset,3:4) = findCRL4Act_v2_reverse(struct('CRLtrace',S(i).crlArea(subset,:),...
                'traceStats', S(i).traceStats(subset,:),'nuc',S(i).crlNuc(subset,:), 'cyto', S(i).crlCyto(subset,:)), ...
                struct('cycling',0,'smooth',7 ,'buff',1 , 'postBuffer', 7 ,'postBufferSec',1,...
                'firstD',0,'secD',-.005,'decrease', .5, 'minSlope', -.01, 'early', 5,'falseCall', .8, 'cytoCheck', 0,'minExp',75,...
                'trunk',.4,'debug',0));
            %%% Divided cells
            subset = find(S(i).isDaughter);
            S(i).POI(subset,3:4) = findCRL4Act_v2_reverse(struct('CRLtrace',S(i).crlArea(subset,:),...
                'traceStats', S(i).traceStats(subset,:),'nuc',S(i).crlNuc(subset,:), 'cyto', S(i).crlCyto(subset,:)), ...
                struct('cycling',1,'smooth',7 ,'buff',1 , 'postBuffer', 7 ,'postBufferSec',3,...
                'firstD',0,'secD',-.005,'decrease', .5, 'minSlope', -.01, 'early', 1,'falseCall', .8, 'cytoCheck', 0,'minExp',50,...
                'trunk',.4,'debug',0));
   
            
            %%% Transformations
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
            
            if settings.poiG2
                [S(i).POI(:,5), badTracesCRL] = findCRLInact(S(i).crlNormAct, S(i).traceStats, S(i).POI(:,4), ...
                    struct('postBuffer',10,'smooth',5,...
                    'buff',10,'thresh',.0025,'preBuffer',5,'lowThresh',.15,...
                    'increase',.005*10,'trunc',.75, 'medfilt',0,'early', 15,'debug',0));
            end
        end
        
    end
    
    %% Extract and gate sig %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(settings.sig)
        S(i).sigNuc = S(i).traceData(:,:,ismember(names,[settings.sig '_mean']));
        S(i).sigArea = (S(i).sigNuc).*S(i).area;
        
        %%% Gate out traces on noise and misexpression %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        
        %%% Transform data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for n = 1:size(S(i).sigArea,1)
            ignore_norm = 5;
            S(i).sigNorm(n,:) =  S(i).sigArea(n,:)./max(S(i).sigArea(n,ignore_norm:end-1));
            S(i).sigSmooth(n,:) = nansmooth(S(i).sigArea(n,:), 5);
            S(i).sigNormSmooth(n,:) = nansmooth(S(i).sigNorm(n,:), 5);
            S(i).sigDiff(n,:) = gradient(S(i).sigSmooth(n,:));
            
            
            if ~settings.quiescent
                startFrame = S(i).POI(n,1);
                if ~isnan(startFrame) & startFrame +1 <= numFrames & ~settings.quiescent
                    S(i).sigNormPostM(n,:) =  S(i).sigArea(n,:)./max(S(i).sigArea(n,startFrame:end-1));
                elseif ~settings.quiescent
                    S(i).sigNormPostM(n,:) =  S(i).sigNorm(n,:);
                end
            end
            S(i).sigNormSmoothPostM(n,:) = nansmoothm(S(i).sigNormPostM(n,:), 7, 'sgolay');
            S(i).sigNormNoise(n,:) = S(i).sigNormPostM(n,:) - S(i).sigNormSmoothPostM(n,:);
            
        end
        
        %%% Find crl inactivation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if settings.poiCrl
            S(i).POI(:,3) = findCRL4Act(S(i).sigNormSmoothPostM,S(i).traceStats, ...
                struct('cycling',~settings.quiescent,'low',.1,'smooth',0,'smoothSec',3, ...
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
    
    
    %% Extract and gate PCNA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(settings.PCNA)
        S(i).pcnaNuc = S(i).traceData(:,:,ismember(names,'PCNA mean'));
        
        %%% Gate out traces on noise and misexpression %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        
        %%% load data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        
        if settings.poiPCNA
            S(i).filterPOI = getPCNAstart_C117(S(i).filterArea(:,:,2),S(i).traceStats,10,10);
            % S(i).filterPOI = getPCNAstart_C114_intensity(S(i).varIntensityNorm(:,:),S(i).traceStats,10,20);
            
        end
    end
    
end

%% Extra transformations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(S)
    
    
end

S = rmfield(S,{'traceData'});

save(fullfile(settings.saveDir, settings.saveName),'S','allIF','settings','-v7.3');




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
