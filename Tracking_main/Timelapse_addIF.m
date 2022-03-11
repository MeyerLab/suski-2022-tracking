function Timelapse_addIF( settings, row,col,site, debug_mode, varargin)
%Timelapse_addIF adding IF imaging to live tracking analysis. Inputs parameter struct s, row,
%          col, site. debug_mode for setting up analysis, varargin contains
%          paths for error logging. Compatible with IX micro and nikon naming schemes. 
%% REQUIRED SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% Paths
% settings.data_path = 'F:\Data\D-Replication Initiation\D121-live\Data\';      % Output location of mat files
% settings.IF_imagesessions = {'L:\D121-IF\'};    
%                                                             % Cell array containing all imaging sessions to be matched
% settings.live_imagepath = 'G:\D121-live\Raw\';    % Live imaging raw image path
% settings.bgcmospath = '';
% settings.crop_save = '';           % Output cropped images, leave empty if no save
% settings.mask_save = 'L:\D121-IF\';           % Output mask, leave empty if no save
% settings.maskname = 'mask';
% settings.primaryMaskRound = 1;
% settings.maskIndex = [1];
% 
% %%% General parameters
% settings.microscope = 'nikon';
% % Formatting code for nikon files, 1$ is timepoint, 2$ is well, 3$ is point
% settings.formatCodeLive = 'Time%1$05d_Well%2$s_Point%2$s_%3$04d_*.tiff';
% settings.formatCodeFixed = 'Time%1$05d_Well%2$s_Point%2$s_%3$04d_*.tiff';
% settings.magnification = 20;     % 10 = 10x or 20 = 20x
% settings.binsize = 1;            % 1 = bin 1 or 2 = bin 2
% settings.postbin = [.5];           % 0 for no software bin, number for scaling factor (e.g. 0.5 for 2x2 bin);
% settings.match10x20x = 1;
% settings.numsites = 9;           % Number of sites imaged for 20x to 10x (double check order is correct)
% settings.scaleLive = 1;          % Scale live to have same pixel size as final IF image
% settings.signals = {{'DAPI','YFP','FarRed'}};  % Each imaging session in a cell inside cell array
% settings.nucLive = 'CFP_';       % Label for IX micro live nucleus (assume first tiff in stack for nikon)
% settings.liveChannel = 4;         % If nikon use this tiff in stack for nuclear live
% settings.manualLiveFrame = [];   % Manually set last frame for live imaging, otherwise use tracedata if empty
% 
% %%% Segmentation parameters
% settings.segmethod = {'concavity'};  % 'log' or 'single' or 'double'
% settings.nucr = 12;                           % 10x bin1 = 12 and 20x bin2 = 12
% settings.debrisarea = 100;                    % 100
% settings.boulderarea = 1500;                  % 1500 
% settings.blobthreshold = -0.03;               % For blobdector 'log' segmentation
% settings.blurradius = 3;                      % Blur for blobdetector
% settings.soliditythresh = 0.80;               % Minimum solidity of nucleus allowed
% settings.badFrameCheck = .25;                 % Min fraction of cells untracked before logging error, 0 for no check
% settings.split_mult = 1;
% 
% %%% Quantification parameters
% settings.bias = {[1 1 1]};
% settings.biasall = {[1 1 1]};
% settings.sigblur = {[0 0 0]};
% settings.signal_foreground = {[0 0 0]};
% settings.bleedthrough = {[0 0 0]};   % Calculate bleedthrough for channel
% settings.bleedthroughslope = {};              % Cell array of cell arrays
% settings.bleedthroughoff = {};                % Cell array of cell arrays
% settings.ringcalc = {[1 1 1]};
% settings.ringthresh = {[0 0 0 ]};    % Threshold for foreground ring pixels
% settings.punctacalc = {[0 0 0]};     % Analyze puncta in nucleus
% settings.punctaThresh = {{[],[],[]},{}};      % Threshold for tophat filter for puncta
% settings.punctatopsize = 2;                   % Top hat filter size
% settings.cytopuncta = {[0 0 0]};
% settings.thickenradius = 2*settings.nucr; 
% 
% settings.localbg = {[0 0 0]};        % Measure local background in box around cell
% settings.minringsize = 100;
% settings.bgcmoscorrection = 1;                % 1 = correct for shading; 0 = dont correct for shading;
% settings.compression = 4;
% settings.bgsubmethod = {{'global nuclear','global cyto', 'global cyto'}}; % 'global nuclear','global cyto','tophat','semi-local nuclear'
% settings.bgperctile = {[25 25 25]};  % Set bg as perctile for each channel
% settings.frameIF=1;
%  
% 
% 
% %%% Tracking parameters
% settings.distthresh = 3*settings.nucr;
% settings.arealowthresh = -.4;
% settings.areahighthresh = .5;
% settings.maxjit = 200;



%% PROCESS IMAGES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
shotLive=[num2str(row),'_',num2str(col),'_',num2str(site)];
if exist(fullfile(settings.data_path,['tracedata_',shotLive,'.mat']))% & ~exist(fullfile(settings.data_path,['IF_',shotLive,'.mat']))
    if debug_mode | nargin-5 == 0
        logFolder = 'test\';
        filePrefix = 'Worker1';
        mkdir(logFolder)
    else
        logFolder = varargin{1};
        filePrefix = varargin{2};
    end
    settings.logFile = [logFolder filePrefix '_log.txt'];
    settings.errorFile = [logFolder filePrefix '_error.txt'];
    
    %%% Setup time output and error log
    timetotal=tic;
    if ~exist(logFolder,'dir')
        mkdir(logFolder)
    end
    errorFileID = fopen(settings.errorFile,'a+');
    logFileID = fopen(settings.logFile,'a+');
    fprintf(logFileID, '%s : Shot %02d_%02d_%02d started ...\r\n',datestr(now,'HH:MM:SS'),row,col,site);
    
%     try
        %%% Load live data
        live = load(fullfile(settings.data_path,['tracedata_',shotLive,'.mat']),'tracedata','jitters');
        [live.totalcells,live.totalframes,~]=size(live.tracedata);
        if ~isempty(settings.manualLiveFrame)
            liveFrame = settings.manualLiveFrame;
        else
            liveFrame = live.totalframes;
        end
        
        %%% Load previous image and create live mask
        switch settings.microscope
            case 'IXM'
                rawLiveDir = fullfile(settings.live_imagepath,shotLive);                
                fileName = [shotLive '_' settings.nucLive,num2str(liveFrame),'.tif'];
                live.rawprev = single(imread(fullfile(rawLiveDir,fileName)));
            case 'nikon'
                rawLiveDir = fullfile(settings.live_imagepath,shotLive);                
                shotSearch = sprintf(settings.formatCodeLive,liveFrame-1,rowColumnTowellName(row,col),site-1);
                fileName = findFile(rawLiveDir,shotSearch);
                live.rawprev = single(imread(fullfile(rawLiveDir,fileName),settings.liveChannel)); 
        end
        
        [live.prevheight,live.prevwidth]=size(live.rawprev);
        live.nuc_mask_prev=threshmask(live.rawprev,3);
        
        
        %%% Process fixed image
        if settings.match10x20x
            [upperleft,upperright,lowerleft,lowerright]=siteChange(settings.numsites,site,settings.microscope);
            
            %%% Process each 20x image for the 10x live frame
            IFdata_upperleft = processImage(row, col, site, upperleft, live, settings, logFileID, errorFileID,debug_mode);
            IFdata_upperright = processImage(row, col, site, upperright, live, settings, logFileID, errorFileID,debug_mode);
            IFdata_lowerleft = processImage(row, col, site, lowerleft, live, settings, logFileID, errorFileID,debug_mode);
            IFdata_lowerright = processImage(row, col, site, lowerright, live, settings, logFileID, errorFileID,debug_mode);
            
            %%% Collect data
            IFdata = ones(size(IFdata_upperleft))*NaN;
            IFdata = addGoodRows(IFdata,IFdata_upperleft);
            IFdata = addGoodRows(IFdata,IFdata_upperright);
            IFdata = addGoodRows(IFdata,IFdata_lowerleft);
            IFdata = addGoodRows(IFdata,IFdata_lowerright);
        else
            %%% Process fixed image for live frame
            IFdata =  processImage(row, col, site, site, live, settings, logFileID, errorFileID,debug_mode);
        end
        
        %%% Save data and settings
        save(fullfile(settings.data_path,['IF_',shotLive,'.mat']),'IFdata');
        header = output_names_multiIF_puncta_livefixed(settings.signals,settings.localbg,settings.ringcalc,settings.punctacalc,settings.punctaThresh);
        save(fullfile(settings.data_path,'settings_live_IF.mat'),'settings','header'); %will overwrite with the most recent
        
        elapsedTime = toc(timetotal);
        fprintf(logFileID, '%s : Shot %02d_%02d_%02d finished in %05.1f sec\r\n\r\n',datestr(now,'HH:MM:SS'),row,col,site,elapsedTime);
%     catch ME
%         elapsedTime = toc(timetotal);
%         fprintf(errorFileID, '%s : ERRROR Shot %02d_%02d_%02d after %05.1f sec\r\n',datestr(now,'HH:MM:SS'),row,col,site,elapsedTime);
%         fprintf(errorFileID,'%s \r\n\r\n',ME.message);
%         fprintf( 'ERRROR Shot %02d_%02d_%02d after %05.1f sec\n',row,col,site,elapsedTime);
%         fprintf('%s\n',ME.message);
%         fprintf('Function:%s\nLine:%d\n',ME.stack(1).name,ME.stack(1).line);
%         
%     end
    fclose(logFileID);
    fclose(errorFileID);
end
end




%% PROCESSING FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [IFdata] = processImage(row, col, liveSite,  fixedSite,live, settings, logFileID, errorFileID,debug_mode)

%%% Print status to log
siteTime = tic;
fprintf(logFileID, '%1$s : Matching %2$02d_%3$02d_%4$02d fixed to %2$02d_%3$02d_%5$02d ...\r\n', ...
    datestr(now,'HH:MM:SS'),row,col,fixedSite,liveSite);

%% LOAD AND SEGMENT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Set paths
shotFixed=[num2str(row),'_',num2str(col),'_',num2str(fixedSite)];

%%% Set channnels
names = settings.signals;
numParams = 8 + 2*length([settings.signals{:}]) + sum([settings.localbg{:}])+ sum([settings.ringcalc{:}]) + 3* length(cell2mat([settings.punctaThresh{:}]));

%%% Set segmentation parameters
nucr = settings.nucr;
debrisarea = settings.debrisarea;
boulderarea = settings.boulderarea;
blobthreshold = settings.blobthreshold;
solidity = settings.soliditythresh;

%% Load auxilliary files %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Load bgcmos
if ~isempty(settings.bgcmospath)
    load([settings.bgcmospath],'cmosoffset');
elseif strcmp(settings.microscope,'nikon')
    cmosoffset = 100;
end

%%% Load bias
for i = 1:length(names)
    for j = 1:length(names{i})
        if settings.biasall{i}(j)
            biasdir = fullfile(settings.IF_imagesessions{i}, 'Bias_all');
        else
            biasdir = fullfile(settings.IF_imagesessions{i}, 'Bias');
        end
        if settings.bias{i}(j)
            load(fullfile(biasdir,[names{i}{j},'_',num2str(fixedSite),'.mat']));
            bias_cell{i}{j} = bias;
        end
    end
end

%% Load images and segment %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(names)
    %%% Load data
    for j = 1:length(names{i})
         switch settings.microscope
            case 'IXM'
               rawFixedDir = fullfile(settings.IF_imagesessions{i},'Raw',shotFixed);                
               fileName = [shotFixed '_' names{i}{j} '_' num2str(settings.frameIF),'.tif'];
               raw{i}{j} = single(imread(fullfile(rawFixedDir,fileName)));
            case 'nikon'
                rawFixedDir = fullfile(settings.IF_imagesessions{i},'Raw',shotFixed);                
                shotSearch = sprintf(settings.formatCodeFixed,settings.frameIF-1,rowColumnTowellName(row,col),fixedSite-1);
                fileName = findFile(rawFixedDir,shotSearch);
               raw{i}{j} = single(imread(fullfile(rawFixedDir,fileName),j));
         end
        
        [height,width]=size(raw{i}{j});
        
        if settings.bgcmoscorrection
            raw{i}{j} = (raw{i}{j}-cmosoffset);
            raw{i}{j}(raw{i}{j}<1) = 1;
        end
        if settings.bias{i}(j)
            raw{i}{j} = raw{i}{j}./bias_cell{i}{j};
        end
        if settings.postbin(i)> 0
            raw{i}{j} = imresize(raw{i}{j},settings.postbin(i));
            [height,width]=size(raw{i}{j});
        end
    end
    
    %% Segment nuclei %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    blurradius = settings.blurradius;
    switch settings.segmethod{i}
        case 'log'
            nuc_mask{i} = blobdetector_4(log(raw{i}{settings.maskIndex(i)}),nucr,blobthreshold,debrisarea);
        case 'thresh'            
            nuc_mask{i} = threshmask(raw{i}{settings.maskIndex(i)},blurradius);
        case 'single'
            nuc_mask{i} = threshmask(raw{i}{settings.maskIndex(i)},blurradius);
            nuc_mask{i} = markershed_filter(nuc_mask{i},round(nucr*2/3),4);
        case 'double'
            nuc_mask{i} = threshmask(raw{i}{settings.maskIndex(i)},blurradius);
            nuc_mask{i} = markershed_filter(nuc_mask{i},round(nucr*2/3),4);
            nuc_mask{i} = secondthresh(raw{i}{settings.maskIndex(i)},blurradius,nuc_mask{i},boulderarea*2);
        case 'double marker'
            nuc_mask{i} = threshmask(raw{i}{settings.maskIndex(i)},blurradius);
            nuc_mask{i} = secondthresh_all(raw{i}{settings.maskIndex(i)},blurradius,nuc_mask{i});
            nuc_mask{i} = markershed_filter(nuc_mask{i},round(nucr*.4),4);
        case 'multithresh'
            nuc_mask{i} = threshmask_multi(raw{i}{settings.maskIndex(i)},blurradius,3, 1e-02);
            nuc_mask{i} = markershed_filter(nuc_mask{i},round(nucr*.4),4);
            %nuc_mask{i} = secondthresh(raw{i}{settings.maskIndex(i)},blurradius,nuc_mask{i},boulderarea);
        case 'concavity'
            imblur = imfilter(raw{i}{settings.maskIndex(i)},fspecial('gaussian',blurradius),'symmetric');
            nuc_mask{i} = ThreshImage(imblur);
            nuc_mask{i} = logical(imfill(nuc_mask{i},'holes'));
    end  
    %%% Split and filter out nuclei only for primary mask
    if i==settings.primaryMaskRound
        whole_mask = logical(nuc_mask{i});
        nuc_mask{i} = bwareaopen(nuc_mask{i},debrisarea);
        
        %%% adaptive bridging settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [bounds,label_mask]=bwboundaries(nuc_mask{i},'noholes');
        mask_props = regionprops(label_mask,'Area');%,'Solidity');
        sizes = [mask_props(:).Area];
        %         solid = [mask_props(:).Solidity];%
        medSize = prctile(sizes,50);
        gate = sizes >= medSize*settings.split_mult;
        
        nuc_mask{i} = segmentdeflections_gated(nuc_mask{i},nucr, bounds,gate);
        nuc_mask{i} = bwareaopen(nuc_mask{i},debrisarea);
        
        clean_mask = excludelargeandwarped_3(nuc_mask{i},boulderarea,solidity);
        clean_mask = imclearborder(clean_mask);
        exclude_mask = nuc_mask{i} - clean_mask;
        nuc_mask{i} = clean_mask;
        
        %% Check for bad frame %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        IF_foreground = sum(nuc_mask{i}(:))/length(nuc_mask{i}(:));
        live_foreground = sum(live.nuc_mask_prev(:))/length(live.nuc_mask_prev(:));
        if  IF_foreground/live_foreground < settings.badFrameCheck
            IFdata = NaN * ones(live.totalcells,numParams);
            
            %% output to log
            fprintf(logFileID, '%1$s : ERROR matching %2$02d_%3$02d_%4$02d fixed to %2$02d_%3$02d_%5$02d. Too many lost cells\r\n\r\n', ...
                datestr(now,'HH:MM:SS'),row,col,fixedSite,liveSite);
            fprintf(errorFileID, '%1$s : ERROR matching %2$02d_%3$02d_%4$02d fixed to %2$02d_%3$02d_%5$02d. Too many lost cells\r\n\r\n', ...
                datestr(now,'HH:MM:SS'),row,col,fixedSite,liveSite);
            return;
        end
    else
        nuc_mask{i} = bwareaopen(nuc_mask{i},debrisarea);
        nuc_mask{i} = imclearborder(nuc_mask{i});
    end
    
end

%%% Clear unnecessary variables
clear bias bias_cell cmosoffset imblur clean_mask

%% Check  segmentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if debug_mode
    keyboard;
end
%{
%%% Check nuclear mask
channeltemp = 1;
extractmask = bwmorph(nuc_mask{channeltemp},'remove');
tempframe = imadjust(mat2gray(raw{channeltemp}{1}));
tempframe(:,:,2) = extractmask;
tempframe(:,:,3) = 0;
figure,imshow(tempframe);

nuc_info = struct2cell(regionprops(nuc_mask{channeltemp},'Area')');
nuc_area = squeeze(cell2mat(nuc_info(1,1,:)));
figure, hist(nuc_area,100);
%}

%% Align images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate jitter between IF rounds
jitmatx = [];
jitmaty = [];

if length(names) > 1
    regheight = round(0.25*height:1:.75*height);
    regwidth = round(0.25*width:1:.75*width);
    for i = 1:length(names)
        [reljitx,reljity] = registerimages(nuc_mask{settings.primaryMaskRound}(regheight,regwidth),nuc_mask{i}(regheight,regwidth));
        jitmatx = [jitmatx; reljitx];
        jitmaty = [jitmaty; reljity];
    end
    
    %%% Error if large offset
    if max(jitmatx) > settings.maxjit | max(jitmaty) > settings.maxjit
        fprintf(errorFileID, '%1$s : WARNING matching %2$02d_%3$02d_%4$02d fixed to %2$02d_%3$02d_%5$02d max jitter too high (%6$d, %7$d)\r\n', ...
            datestr(now,'HH:MM:SS'),row,col,fixedSite,liveSite,max(jitmatx), max(jitmaty));
        regheight = 1:height;
        regwidth = 1:width;
        jitmatx = [];
        jitmaty = [];
        for i = 1:length(names)
            [reljitx,reljity] = registerimages(nuc_mask{settings.primaryMaskRound}(regheight,regwidth),nuc_mask{i}(regheight,regwidth));
            jitmatx = [jitmatx; reljitx];
            jitmaty = [jitmaty; reljity];
        end
        fprintf(errorFileID, '%1$s : WARNING matching %2$02d_%3$02d_%4$02d fixed to %2$02d_%3$02d_%5$02d Fixed to (6$d, 7$d)\r\n', ...
            datestr(now,'HH:MM:SS'),row,col,fixedSite,liveSite,max(jitmatx), max(jitmaty));
        
    end
else
    jitmatx = [jitmatx 0];
    jitmaty = [jitmaty 0];
end



%%% Align and crop IF rounds
cropcoors = getcropcoors([width height ], jitmatx, jitmaty);
for i = 1:length(names)
    for j = 1:length(names{i})
        raw{i}{j} = raw{i}{j}(cropcoors(i+1,1):cropcoors(i+1,2),cropcoors(i+1,3):cropcoors(i+1,4));
        
    end
end
nuc_mask_aligned = nuc_mask{settings.primaryMaskRound}(cropcoors(1,1):cropcoors(1,2),cropcoors(1,3):cropcoors(1,4));
whole_mask_aligned = whole_mask(cropcoors(1,1):cropcoors(1,2),cropcoors(1,3):cropcoors(1,4));
exclude_mask_aligned = exclude_mask(cropcoors(1,1):cropcoors(1,2),cropcoors(1,3):cropcoors(1,4));
[height,width]=size(nuc_mask_aligned);

%%% Clear unnecessary variables
clear nuc_mask normnuc normraw sumraw exclude_mask

%% Subtract background %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
compression = settings.compression;
for i=1:length(names)
    for j = 1:length(names{i})
        %%% Set foreground mask
        if settings.signal_foreground{i}(j) 
            %calculate foreground based on sum of signal
            normnuc = (raw{settings.primaryMaskRound}{settings.maskIndex(settings.primaryMaskRound)}-prctile(raw{settings.primaryMaskRound}{settings.maskIndex(settings.primaryMaskRound)}(:),1));
            normnuc(normnuc < 1) = 1;
            normnuc = normnuc/prctile(normnuc(:),99);
            normraw = (raw{i}{j}-prctile(raw{i}{j}(:),1));
            normraw(normraw < 1) = 1;
            normraw = normraw/prctile(normraw(:),99);
            sumraw = normnuc + normraw;
            foreground_aligned = threshmask(sumraw,blurradius);
            % account for large foregrounds
            if sum(foreground_aligned(:))/length(foreground_aligned(:)) > .9
                foreground_aligned = zeros(height, width);
            end
        else
            foreground_aligned = zeros(height, width);
        end
        
        %Mask background
        mask_total  =  whole_mask_aligned | foreground_aligned;
        nanmask_aligned{i}{j}  = imdilate(mask_total,strel('disk',nucr));
        nanmaskcyto_aligned{i}{j} = imdilate(mask_total,strel('disk',nucr*2));
        
        
        %%% Blur signal
        if settings.sigblur{i}(j) > 0
            blur = imfilter(raw{i}{j},fspecial('disk',settings.sigblur{i}(j)),'symmetric');
        else
            blur  =  raw{i}{j};
        end
        
        switch settings.bgsubmethod{i}{j}
            case 'global nuclear'
                real{i}{j} = bgsubmasked_global_NR(blur,nanmask_aligned{i}{j},1,compression,settings.bgperctile{i}(j));
                if settings.punctacalc{i}(j)
                    unblurred{i}{j} = bgsubmasked_global_NR(raw{i}{j},nanmask_aligned{i}{j},1,compression,settings.bgperctile{i}(j));
                end
            case 'global cyto'
                real{i}{j} = bgsubmasked_global_2(blur,nanmaskcyto_aligned{i}{j},1,compression,settings.bgperctile{i}(j));
                if settings.punctacalc{i}(j)
                    unblurred{i}{j} = bgsubmasked_global_2(raw{i}{j},nanmaskcyto_aligned{i}{j},1,compression,settings.bgperctile{i}(j));
                end
            case 'tophat'
                real{i}{j} = imtophat(blur,strel('disk',nucr,0));
                if settings.punctacalc{i}(j)
                    unblurred{i}{j} = imtophat(raw{i}{j},strel('disk',nucr,0));
                end
            case 'semi-local nuclear'
                real{i}{j} = bgsubmasked_global_2(blur,nanmask_aligned{i}{j},11,compression,settings.bgperctile{i}(j));
                if settings.punctacalc{i}(j)
                    unblurred{i}{j} = bgsubmasked_global_2(raw{i}{j},nanmask_aligned{i}{j},11,compression,settings.bgperctile{i}(j));
                end
            case 'none'
                real{i}{j} = blur;
                if settings.punctacalc{i}(j)
                    unblurred{i}{j} = raw{i}{j};
                end
        end
    end

end

if ~any(cell2mat(settings.localbg))
    clear nanmask_aligned nanmaskcyto_aligned
end
%% correct IF for bleedthrough
for i=1:length(names)
    for j = 1:length(names{i})
        if settings.bleedthrough{i}(j)
            realbleedthrough = zeros(size(real{i}{j}));
            %unblurbleedthrough = zeros(size(unblurred{i}{j}));
            for k = 1:length(names)
                for l = 1:length(names{k})
                    if settings.bleedthroughslope{i}{j}{k}(l) > 0
                        bleedreal = real{k}{l}*settings.bleedthroughslope{i}{j}{k}(l) + settings.bleedthroughoff{i}{j}{k}(l);
                        %bleedunblurred = unblurred{k}{l}*settings.bleedthroughslope{i}{j}{k}(l) + settings.bleedthroughoff{i}{j}{k}(l);
                        realbleedthrough = realbleedthrough + bleedreal;
                        %unblurbleedthrough = unblurbleedthrough + bleedunblurred;
                        clear bleedreal bleedunblurred;
                    end
                end
            end
            real{i}{j} = real{i}{j}-realbleedthrough;
            %unblurred{i}{j} = unblurred{i}{j}-unblurbleedthrough;
            unblurred{i}{j} = unblurred{i}{j}-realbleedthrough;
            
            clear realbleedthrough unblurbleedthrough;
        end
    end
end


%% Check alignment segmentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if debug_mode & length(names) > 1
    keyboard;
end
%{
%%% Check nuclear mask with primary nuclear channel
extractmask = bwmorph(nuc_mask_aligned,'remove');
tempframe = imadjust(mat2gray(raw{settings.primaryMaskRound}{settings.maskIndex(settings.primaryMaskRound)}));
tempframe(:,:,2) = extractmask;
tempframe(:,:,3) = 0;
figure,imshow(tempframe);

%%% Check nuclear areas
nuc_info = struct2cell(regionprops(nuc_mask_aligned,'Area')');
nuc_area = squeeze(cell2mat(nuc_info(1,1,:)));
hist(nuc_area,100);

%%% Check masks of each session
for session = 1:length(names)
    extractmask = bwmorph(nuc_mask_aligned,'remove');
    tempframe = imadjust(mat2gray(raw{session}{settings.maskIndex(session)}));
    tempframe(:,:,2) = extractmask;
    tempframe(:,:,3) = 0;
    figure,imshow(tempframe);
end

%%% Check foreground
round = 2;
channel = 2;
extractmask = bwmorph(foreground_aligned{round}{channel},'remove');
tempframe = imadjust(mat2gray(raw{round}{channel}));
tempframe(:,:,2) = extractmask;
tempframe(:,:,3) = 0;
figure,imshow(tempframe);

%%% Overlay nuclear images
tempframe=imadjust(mat2gray(real{1}{settings.maskIndex(1)}));
tempframe(:,:,2)=imadjust(mat2gray(real{2}{settings.maskIndex(2)}));
tempframe(:,:,3)=0;
figure,imshow(tempframe)
%}

%%% Clear unnecessary variables
clear nuc_mask blur mask_total raw whole_mask foreground_aligned 

%% Puncta processing
punctaFilt = {};
for i = 1:length(names)
    for j = 1:length(names{i})
        if settings.punctacalc{i}(j)
            punctaFilt{i}{j} = imtophat(unblurred{i}{j}, strel('disk',settings.punctatopsize,0));
            if debug_mode
                keyboard;
                %{
                imtool(unblurred{i}{j})
                imtool(punctaFilt{i}{j}.*erodeMask)
                %}
            end

        end
    end
end

%% Save cropped/corrected images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(settings.crop_save) & ~debug_mode
    save_dir = fullfile(settings.crop_save,shotFixed);
    for i = 1:length(names)
        for j = 1:length(names{i})
            if ~exist(save_dir)
                mkdir(save_dir)
            end
            save_name = fullfile(save_dir, [shotFixed,'_' names{i}{j} '_round' num2str(i) '.tif']);
            imwrite(uint16(real{i}{j}/4), save_name);
        end
    end
    
end

if ~isempty(settings.mask_save)
    save_dir = fullfile(settings.mask_save,shotFixed);
    if ~exist(save_dir)
        mkdir(save_dir)
    end
    %%% Save mask
    extractmask = bwmorph(nuc_mask_aligned,'remove');
    save_name = fullfile(save_dir, [shotFixed,'_' settings.maskname '.tif']);
    imwrite(uint16(extractmask),save_name);
end

%% MATCH LIVE AND FIXED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pad IF images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[cropheight,cropwidth] = size(real{1}{1});
diff_y = live.prevheight*settings.scaleLive-cropheight;
tophalf = floor(diff_y/2);
bottomhalf = diff_y - tophalf;
diff_x = live.prevwidth*settings.scaleLive-cropwidth;
lefthalf = floor(diff_x/2);
righthalf = diff_x - lefthalf;
% [cropheight,cropwidth] = size(real{1}{1});
% tophalf = floor((live.prevheight*settings.scaleLive-cropheight)/2);
% if mod(cropheight,2) == 1
%     bottomhalf = tophalf+1;
% else
%     bottomhalf = tophalf;
% end
% lefthalf = floor((live.prevwidth*settings.scaleLive-cropwidth)/2);
% if mod(cropwidth,2) == 1
%     righthalf = lefthalf + 1;
% else
%     righthalf = lefthalf;
% end

%%% Pad masks and images
nuc_mask_padded = padarray(nuc_mask_aligned,[tophalf lefthalf],'pre');
nuc_mask_padded = padarray(nuc_mask_padded,[bottomhalf righthalf],'post');
exclude_mask_padded = padarray(exclude_mask_aligned,[tophalf lefthalf],'pre');
exclude_mask_padded = padarray(exclude_mask_padded,[bottomhalf righthalf],'post');

for i=1:length(names)
    for j = 1:length(names{i})
        real_padded{i}{j} = padarray(real{i}{j},[tophalf lefthalf],'pre');
        real_padded{i}{j} = padarray( real_padded{i}{j},[bottomhalf righthalf],'post');
        if settings.localbg{i}(j)
            nanmask_padded{i}{j} = padarray(nanmask_aligned{i}{j},[tophalf lefthalf],'pre');
            nanmask_padded{i}{j} = padarray( nanmask_padded{i}{j},[bottomhalf righthalf],'post');
        end
        if settings.punctacalc{i}(j)
            punctaFilt_padded{i}{j} = padarray(punctaFilt{i}{j},[tophalf lefthalf],'pre');
            punctaFilt_padded{i}{j} = padarray(punctaFilt_padded{i}{j},[bottomhalf righthalf],'post');
        end
     
    end
end

%%% Clear unnecessary variables
clear nuc_mask_aligned nanmask_aligned real unblurred punctaMask punctaIntensity nanmask_aligned nanmaskcyto_aligned exclude_mask_aligned
 
%% Align live and fixed iamges %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Extract IF cell features to match to live data
nuc_label = bwlabel(nuc_mask_padded);
nuc_info = struct2cell(regionprops(nuc_mask_padded,real_padded{settings.primaryMaskRound}{settings.maskIndex(settings.primaryMaskRound)},'Area','Centroid','MeanIntensity')');
nuc_area = squeeze(cell2mat(nuc_info(1,1,:)));
nuc_center = squeeze(cell2mat(nuc_info(2,1,:)))';
nuc_density = squeeze(cell2mat(nuc_info(3,1,:)));
nuc_mass = nuc_density.*nuc_area;

%%% Calculate jitter between IF and live
[newheight, newwidth] = size(nuc_mask_padded);
newheight = newheight/settings.scaleLive;
newwidth = newwidth/settings.scaleLive;
if settings.match10x20x
    regheight = 1:1*newheight;
    regwidth = 1:1*newwidth;
else
    regheight = round(0.25*newheight:1:.75*newheight);
    regwidth = round(0.25*newwidth:1:.75*newwidth);
end
reg_nuc_mask_padded = imresize(nuc_mask_padded,1/settings.scaleLive);
[reljitx, reljity] = registerimages(live.nuc_mask_prev(regheight,regwidth),reg_nuc_mask_padded(regheight,regwidth));
reljitx = reljitx *settings.scaleLive;
reljity = reljity *settings.scaleLive;
[padheight, padwidth] = size(nuc_mask_padded);
jitcoors = getcropcoors([padwidth padheight],reljitx,reljity);

%%% Align IF masks and images to live coordinates
reljitter = [reljitx, reljity];
prevjitter = live.jitters(live.totalframes,:)*settings.scaleLive;
IFjitter = prevjitter + reljitter;
nuc_center(:,1) = nuc_center(:,1) + IFjitter(1);
nuc_center(:,2) = nuc_center(:,2) + IFjitter(2);

%% Match live and fixed cells and correct labeling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
initdata = [nuc_center(:,1), nuc_center(:,2), nuc_area, nuc_mass];
debugpackage = {live.rawprev, live.nuc_mask_prev, nuc_mask_padded, prevjitter, reljitter};
livetrack = squeeze(live.tracedata(:,live.totalframes,1:4));
livetrack(:,1:2) = livetrack(:,1:2)*settings.scaleLive;
livetrack(:,3:4) = livetrack(:,3:4)*settings.scaleLive^2;

[tracked,nuc_label_track, areas] = adaptivetrack_IF_NR(livetrack, initdata, nuc_label, nucr,...
    settings.distthresh, settings.arealowthresh, settings.areahighthresh, debugpackage);

%%% Calculate number lost cells and record in log file
if settings.match10x20x
    approx_perc_lost_cells = 100*(sum(~isnan(live.tracedata(:,live.totalframes,1)))/4 - sum(~isnan(tracked(:,1))))/...
        sum(~isnan(live.tracedata(:,live.totalframes,1)));
else
    approx_perc_lost_cells = 100*(sum(~isnan(live.tracedata(:,live.totalframes,1))) - sum(~isnan(tracked(:,1))))/...
        sum(~isnan(live.tracedata(:,live.totalframes,1)));   
end
fprintf(logFileID, '%1$s : Lost ~%2$4.1f perc of cells\r\n', ...
    datestr(now,'HH:MM:SS'),approx_perc_lost_cells);
if approx_perc_lost_cells > 10
    fprintf(errorFileID, '%1$s : WARNING matching %2$02d_%3$02d_%4$02d fixed to %2$02d_%3$02d_%5$02d Lost ~%6$4.1f perc of cells\r\n', ...
        datestr(now,'HH:MM:SS'),row,col,fixedSite,liveSite, approx_perc_lost_cells);
end

numcells = numel(tracked(:,1));
nuc_center = tracked(:,[1 2]);
nuc_area = tracked(:,3);
nuc_mass = tracked(:,4);
cellid = find(~isnan(tracked(:,1)));
numlivecells = numel(cellid);
nuc_info = regionprops(nuc_label_track, 'PixelIdxList', 'Area', 'Centroid');

tracked_obj = nuc_label_track > 0;
untracked_obj = nuc_mask_padded & ~tracked_obj | exclude_mask_padded;
untracked_label = bwlabel(untracked_obj);
untracked_label(untracked_label > 0) = untracked_label(untracked_label > 0) + numcells;
all_obj = nuc_label_track + untracked_label;

%% Check tracking quality %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if debug_mode
    keyboard;
end
%{
%Create visualization variables
prev_mask_jit = imresize(live.nuc_mask_prev,settings.scaleLive);
prev_mask_jit = prev_mask_jit(jitcoors(1,1):jitcoors(1,2),jitcoors(1,3):jitcoors(1,4));
nuc_mask_jit = nuc_mask_padded(jitcoors(2,1):jitcoors(2,2),jitcoors(2,3):jitcoors(2,4));
for i=1:length(names)
    for j = 1:length(names{i})
        real_padded_jit{i}{j} = real_padded{i}{j}(jitcoors(2,1):jitcoors(2,2),jitcoors(2,3):jitcoors(2,4));
    end
end
tracked_cells = length(unique(nuc_label));
nuc_label_track_jit = nuc_label_track(jitcoors(2,1):jitcoors(2,2),jitcoors(2,3):jitcoors(2,4));

%Check tracking (dots are previously recorded, red are untracked, yellow tracked
extractmask1=bwmorph(nuc_mask_padded,'remove');
extractmask2=bwmorph(nuc_label_track,'remove');
[tempheight,tempwidth]=size(extractmask1);
tempframe=zeros(tempheight,tempwidth,3);
tempframe(:,:,1)=extractmask1;
tempframe(:,:,2)=extractmask2;
tempframe(:,:,3)=imadjust(mat2gray(real_padded{1}{1}));
figure,imshow(tempframe);
hold on
scatter(livetrack(:,1)-IFjitter(1),livetrack(:,2)-IFjitter(2),'w')

%Check jitter
extractmask=bwmorph(nuc_mask_jit,'remove');
[tempheight,tempwidth]=size(extractmask);
tempframe=zeros(tempheight,tempwidth,3);
tempframe(:,:,1)=prev_mask_jit;
tempframe(:,:,2)=nuc_mask_jit;
figure,imshow(tempframe);
hold on

%%% overlay nuclear images %%%%%%%%%%%%%%%%%
[rawprevcrop,raw1crop]=cropboth(imresize(live.rawprev,settings.scaleLive),real_padded{settings.primaryMaskRound}{settings.maskIndex(settings.primaryMaskRound)},reljitx,reljity);
tempframe=imadjust(mat2gray(rawprevcrop));
tempframe(:,:,2)=imadjust(mat2gray(raw1crop));
tempframe(:,:,3)=0;
figure,imshow(tempframe)

%Histogram of area changes
figure,
histogram(areas)

%Check tracking (dots previously recorded)
extractmask=bwmorph(nuc_mask_padded,'remove');
[tempheight,tempwidth]=size(extractmask);
tempframe=zeros(tempheight,tempwidth,3);
tempframe(:,:,1)=imadjust(mat2gray(real_padded{settings.primaryMaskRound}{settings.maskIndex(settings.primaryMaskRound)}));
tempframe(:,:,2)=extractmask;
%tempframe(:,:,3)=marker_mask;
figure,imshow(tempframe);
hold on
scatter(livetrack(:,1)-IFjitter(1),livetrack(:,2)-IFjitter(2),'w')

%Double check alignment of live and previously tracked cells
figure, imshow(imresize(live.rawprev,settings.scaleLive),[]);
hold on
scatter((live.tracedata(:,live.totalframes,1)*settings.scaleLive-prevjitter(1)),(live.tracedata(:,live.totalframes,2)*settings.scaleLive-prevjitter(2)),'w')

%}

%%% Clear unnecessary variables
clear debugpackage nuc_label


%% Create secondary masks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Cyto ring
if any(cell2mat(settings.ringcalc))
    if (settings.magnification == 10 & settings.binsize == 1) | ...
            (settings.magnification == 20 & (settings.binsize == 2 | settings.postbin))
        innerrad = 1; outerrad = 5;
    else
        innerrad = 2; outerrad = 10;
    end
    ring_label = getcytoring_thicken(all_obj,innerrad,outerrad,real_padded{1}{1});
    ring_info = regionprops(ring_label,'PixelIdxList');
    
    %%% Check cytoring
    if debug_mode
        keyboard;
    end
    %{
            %%% debugging: view images %%%%%%%%%%
            temp_ring  =  ring_label > 0;
            extractmask = bwmorph(temp_ring,'remove');
            signal  =  real_padded{1}{2};
            intensities  =  signal.*temp_ring;
            [height,width] = size(signal);
            tempframe = zeros(height,width,3);
            tempframe(:,:,1) = imadjust(mat2gray(signal));
            tempframe(:,:,2) = extractmask;;
            %tempframe(:,:,3) = marker_mask;
            figure,imshow(tempframe);
    %}
end

%%% Whole "cytoplasm" for FISH
if any(cell2mat(settings.punctacalc))
    if any(cell2mat(settings.cytopuncta))
        cyto_label = labelthicken_better(all_obj,settings.thickenradius);
        cyto_info = regionprops(cyto_label, 'PixelIdxList','Area');
    end
    
    %%% Check cyto mask
    if debug_mode
        keyboard;
    end
    %{
            %%% debugging: view images %%%%%%%%%%
            temp_im  =  cyto_label > 0;
            extractmask = bwmorph(temp_im,'remove');
            signal  =  real_padded{1}{1};
            intensities  =  signal.*temp_ring;
            [height,width] = size(signal);
            tempframe = zeros(height,width,3);
            tempframe(:,:,1) = imadjust(mat2gray(signal));
            tempframe(:,:,2) = extractmask;;
            %tempframe(:,:,3) = marker_mask;
            figure,imshow(tempframe);
    
            output_im = cyto_label .* punctaMask_padded{1}{2}{1};
            figure, imagesc(output_im);
    %}
end

%%% Clear unnecessary variables
clear  nuc_label tracked_obj untracked_obj untracked_label exclude_mask_padded

%% MEASSURE FEATURES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize storage variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nanvec = ones(numcells,1)*NaN;
sigringmedian{1} = [];
localbg{1} = [];
punctaArea{1}{1} = [];
punctaIntensity{1}{1} = [];
punctaNumber{1}{1} = [];
IF_position  = [nanvec nanvec];

for i = 1:length(names)
    for j = 1:length(names{i})
        sigmean{i}{j} = nanvec;
        sigmedian{i}{j} = nanvec;
        if settings.ringcalc{i}(j)
            sigringmedian{i}{j} = nanvec;
        end
        if settings.localbg{i}(j)
            localbg{i}{j} = nanvec;
        end
        if settings.punctacalc{i}(j)
            punctaArea{i}{j} = repmat(nanvec,1,length(settings.punctaThresh{i}{j}));
            punctaIntensity{i}{j} = repmat(nanvec,1,length(settings.punctaThresh{i}{j}));
            punctaNumber{i}{j} = repmat(nanvec,1,length(settings.punctaThresh{i}{j}));
        end
    end
end
%% Measure cell features%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(names)
    for j = 1:length(names{i})
        real_temp = real_padded{i}{j};
        if settings.localbg{i}(j)
            real_temp_masked = real_temp;
            real_temp_masked(nanmask_padded{i}{j}) = NaN;
        end       
        for n = 1:numlivecells
            cc=cellid(n);
            sigmean{i}{j}(cc) = mean(real_temp(nuc_info(cc).PixelIdxList));
            sigmedian{i}{j}(cc) = median(real_temp(nuc_info(cc).PixelIdxList));
            IF_position(cc,:) = nuc_info(cc).Centroid -  [lefthalf tophalf];
            if settings.localbg{i}(j)
                localbg{i}{j}(cc) = getlocalbg(real_temp_masked, nuc_info(cc).Centroid,100, 50,.15);
            end
            if settings.ringcalc{i}(j)
                if cc<numel(ring_info)
                    ringall = real_temp(ring_info(cc).PixelIdxList);
                    ringall(ringall>prctile(ringall,98)) = [];
                    ringforeground = ringall(ringall>settings.ringthresh{i}(j));
                    if numel(ringforeground)<settings.minringsize
                        ringforeground = ringall;
                    end
                    if numel(ringall)>=settings.minringsize
                        sigringmedian{i}{j}(cc) = nanmedian(ringforeground);
                    end
                end
            end
        end
    end
end

%% Measure puncta features%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(names)
    for j = 1:length(names{i})
        if settings.punctacalc{i}(j)
            if(settings.cytopuncta{i}(j))
                mask_info = cyto_info;
            else
                mask_info = nuc_info;
            end
            
            for t = 1:length(settings.punctaThresh{i}{j})
                punctaMask = punctaFilt_padded{i}{j} > settings.punctaThresh{i}{j}(t);
                punctaIntensity_Masked = punctaFilt_padded{i}{j} .* punctaMask;
                puncta_label = bwlabel(punctaMask);
                
                %%% Insert code to output puncta masks here %%%
                for n = 1:numlivecells
                    cc=cellid(n);
                    punctaArea{i}{j}(cc,t) = sum(punctaMask(mask_info(cc).PixelIdxList));
                    punctaIntensity{i}{j}(cc,t) = sum(punctaIntensity_Masked(mask_info(cc).PixelIdxList));
                    numPuncta = length(unique(puncta_label(mask_info(cc).PixelIdxList))) - 1;
                    punctaNumber{i}{j}(cc,t) = numPuncta;
                    
                end
            end
        end
    end
end

%%% Compile data
IFdata = [IF_position(:,1),IF_position(:,2),nuc_area,repmat(IFjitter,size(IF_position,1),1),...
    repmat(row,size(IF_position(:,1),1),1),...
    repmat(col,size(IF_position(:,1),1),1),...
    repmat(fixedSite,size(IF_position(:,1),1),1),...
    cell2mat([sigmean{:}]),cell2mat([sigmedian{:}]),cell2mat([localbg{:}]),...
    cell2mat([sigringmedian{:}]),...
    cell2mat([punctaArea{:}]), cell2mat([punctaIntensity{:}]), cell2mat([punctaNumber{:}])];

%%% Output to log file
elapsedSiteTime = toc(siteTime);
fprintf(logFileID, '%1$s : Finished successfully, %2$05.1f sec elapsed\r\n', ...
    datestr(now,'HH:MM:SS'),elapsedSiteTime);
end




function [upperleft,upperright,lowerleft,lowerright]=siteChange(numsites,site, microscope)
if numsites==4 & strcmp(microscope,'IXM')
    switch site
        case 1
            upperleft=1;
            upperright=2;
            lowerleft=5;
            lowerright=6;
        case 2
            upperleft=3;
            upperright=4;
            lowerleft=7;
            lowerright=8;
        case 3
            upperleft=9;
            upperright=10;
            lowerleft=13;
            lowerright=14;
        case 4
            upperleft=11;
            upperright=12;
            lowerleft=15;
            lowerright=16;
    end
end

if numsites==4 & strcmp(microscope,'nikon')
    switch site
        case 1
            upperleft=1;
            upperright=2;
            lowerleft=8;
            lowerright=7;
        case 2
            upperleft=3;
            upperright=4;
            lowerleft=6;
            lowerright=5;
        case 3
            upperleft=11;
            upperright=12;
            lowerleft=14;
            lowerright=13;
        case 4
            upperleft=9;
            upperright=10;
            lowerleft=16;
            lowerright=15;
    end
end


if numsites==9 & strcmp(microscope,'IXM')
    switch site
        case 1
            upperleft=1;
            upperright=2;
            lowerleft=7;
            lowerright=8;
        case 2
            upperleft=3;
            upperright=4;
            lowerleft=9;
            lowerright=10;
        case 3
            upperleft=5;
            upperright=6;
            lowerleft=11;
            lowerright=12;
        case 4
            upperleft=13;
            upperright=14;
            lowerleft=19;
            lowerright=20;
        case 5
            upperleft=15;
            upperright=16;
            lowerleft=21;
            lowerright=22;
        case 6
            upperleft=17;
            upperright=18;
            lowerleft=23;
            lowerright=24;
        case 7
            upperleft=25;
            upperright=26;
            lowerleft=31;
            lowerright=32;
        case 8
            upperleft=27;
            upperright=28;
            lowerleft=33;
            lowerright=34;
        case 9
            upperleft=29;
            upperright=30;
            lowerleft=35;
            lowerright=36;
    end
end

if numsites==9 & strcmp(microscope,'nikon')
    switch site
        case 1
            upperleft=1;
            upperright=2;
            lowerleft=12;
            lowerright=11;
        case 2
            upperleft=3;
            upperright=4;
            lowerleft=10;
            lowerright=9;
        case 3
            upperleft=5;
            upperright=6;
            lowerleft=8;
            lowerright=7;
        case 4
            upperleft=17;
            upperright=18;
            lowerleft=20;
            lowerright=19;
        case 5
            upperleft=15;
            upperright=16;
            lowerleft=22;
            lowerright=21;
        case 6
            upperleft=13;
            upperright=14;
            lowerleft=24;
            lowerright=23;
        case 7
            upperleft=25;
            upperright=26;
            lowerleft=36;
            lowerright=35;
        case 8
            upperleft=27;
            upperright=28;
            lowerleft=34;
            lowerright=33;
        case 9
            upperleft=29;
            upperright=30;
            lowerleft=32;
            lowerright=31;
    end
end


if numsites==6 & strcmp(microscope,'IXM')
    switch site
        case 1
            upperleft=1;
            upperright=2;
            lowerleft=7;
            lowerright=8;
        case 2
            upperleft=3;
            upperright=4;
            lowerleft=9;
            lowerright=10;
        case 3
            upperleft=5;
            upperright=6;
            lowerleft=11;
            lowerright=12;
        case 4
            upperleft=13;
            upperright=14;
            lowerleft=19;
            lowerright=20;
        case 5
            upperleft=15;
            upperright=16;
            lowerleft=21;
            lowerright=22;
        case 6
            upperleft=17;
            upperright=18;
            lowerleft=23;
            lowerright=24;
    end
end

if numsites==6 & strcmp(microscope,'nikon')
    switch site
        case 1
            upperleft=1;
            upperright=2;
            lowerleft=12;
            lowerright=11;
        case 2
            upperleft=3;
            upperright=4;
            lowerleft=10;
            lowerright=9;
        case 3
            upperleft=5;
            upperright=6;
            lowerleft=8;
            lowerright=7;
        case 4
            upperleft=17;
            upperright=18;
            lowerleft=20;
            lowerright=19;
        case 5
            upperleft=15;
            upperright=16;
            lowerleft=22;
            lowerright=21;
        case 6
            upperleft=13;
            upperright=14;
            lowerleft=24;
            lowerright=23;
    end
end

end




function updateddata=addGoodRows(orgdata,newdata)
goodrows=~isnan(newdata(:,1));
updateddata=orgdata;
updateddata(goodrows,:)=newdata(goodrows,:);
end