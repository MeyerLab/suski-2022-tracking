function Timelapse_singleFrame(s,row,col,site,debug_mode,varargin)
%Timelapse_singleFrame timelapse single frame analysis. Inputs parameter struct s, row,
%          col, site. debug_mode for setting up analysis, varargin is
%          ignored. Compatible with IX micro and nikon naming schemes.
%% REQUIRED SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Paths
% %%% Paths
% s.experiment_name='C196-control';                     % Experiment name
% s.image_drive = 'J:\Nikon\C196_two_stage_control\';   % Location of images base folder
% s.savepath=fullfile('F:\Data\C-Cdt1\',s.experiment_name,'Data'); % Output location of mat files
% s.imagepath=fullfile(s.image_drive,'Raw');                       % All sites in same folder for nikon
% s.biaspath=fullfile(s.image_drive,'Bias');                       % Location of bias files
% s.maskpath=fullfile(s.image_drive,'Mask');                       % Future location of mask saving
% s.bgcmospath='';
% % Formatting code for nikon files, 1$ is timepoint, 2$ is well, 3$ is point
% s.formatCode = 'Time%1$05d_Well%2$s_Point%2$s_%3$04d_*.tiff';
% s.microscope = 'nikon';
% 
% %%% General parameters
% s.frame = 1; 
% s.magnification=20;              % 10=10x or 20=20x
% s.binsize=1;                     % 1=bin 1 or 2=bin 2
% s.postbin=0;                     % 0 for no software bin, number for scaling factor (e.g. 0.5 for 2x2 bin);
% s.signals={'CFP_','RFP_'};
% s.maskwrite = 1;                 % 1 = save an image with the nucmask; 0 = dont save an image with the nucmask
% s.maskname = 'nucedge_';
% 
% 
% %%% Quantification parameters
% s.bgcmoscorrection = 1;
% s.bias = [1 1 1];
% s.signal_foreground = [0 0 0];
% s.bgsubmethod = {'global nuclear','global nuclear','global nuclear'}; 
%                                  % Options:'global nuclear','global cyto','tophat','semi-local nuclear', 'none'
% s.compression = 4;               % For bg subtraction
% s.bgprctile = [25 25 25];        % Set bg as perctile for each channel
% s.sigblur = [3 3 3];
% s.localbgmeasure = [0 0 0];      % Measure local background in box around cell
% s.ringcalc = [0 1 ];
% s.ringthresh = [0 0];            % Threshold for foreground ring pixels
% s.punctacalc = 0;                % Analyze puncta in nucleus
% s.punctaThresh = [125 150 200];
% s.varThresh = [75 100 125];
% 
% %%% Segmentation parameters
% s.firstsegmethod = 'concavity';  % Options: 'concavity', 'log', 'single adapt', 'double adapt', 'single', 'double'
% s.nucr = 24;                     % 10x bin1 = 12 and 20x bin2 = 12
% s.blurradius  =  3;              % 10x: 3
% s.soliditythresh = 0.8;          % Minimum solidity of nucleus allowed
% s.debrisarea = 400;              % 10x = 100 % 150 if no mitosis
% s.boulderarea = 6000;            % 10x = 1500
% s.blobthreshold = -0.02;         % For blobdector 'log' segmentation


%% ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timetotal = tic;
shot = [num2str(row),'_',num2str(col),'_',num2str(site)];
if ~exist(fullfile(s.savepath,['tracedata_',shot,'.mat']),'file')
    
    %% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Initialize path variables
    switch s.microscope
        case 'IXM'
            rawdir = fullfile(s.imagepath,shot,[shot '_']);
        case 'nikon'
            rawdir = fullfile(s.imagepath,shot);
    end
    maskdir = fullfile(s.maskpath,shot);
    nucname = s.signals{1};
    parameternum = 4+2*length(s.signals)+sum(s.ringcalc)+sum(s.localbgmeasure);
    if s.punctacalc
        parameternum = parameternum + 4 + 2*length(s.punctaThresh) + 2*length(s.varThresh);
    end
    if ~exist(s.savepath,'dir')
        mkdir(s.savepath);
    end
    if ~exist(maskdir,'dir') && s.maskwrite
        mkdir(maskdir);
    end
    
    %%% Load bgcmos/s.bias
    if ~isempty(s.bgcmospath)
        load(s.bgcmospath,'cmosoffset');
        bgcmos  =  cmosoffset;
    elseif strcmp(s.microscope,'nikon')
        bgcmos = 100;
    end
    
    for i = 1:length(s.signals)
        if s.bias(i)
            load(fullfile(s.biaspath,[s.signals{i},num2str(site),'.mat']));
            bias_cell{i} = bias;
        end
    end
    
    %% Imaging processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% Load images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for j = 1:length(s.signals)
        switch s.microscope
            case 'IXM'
                rawdir = fullfile(s.imagepath,shot);
                fileName = [shot '_' s.signals{j},num2str(s.frame),'.tif'];
                raw{j} = single(imread(fullfile(rawdir,fileName)));
            case 'nikon'
                rawdir = fullfile(s.imagepath,shot);
                shotSearch = sprintf(s.formatCode,s.frame-1,rowColumnTowellName(row,col),site-1);
                fileName = findFile(rawdir,shotSearch);
                raw{j} = single(imread(fullfile(rawdir,fileName),j));
        end
        [height,width]=size(raw{j});
        
        if s.bgcmoscorrection
            raw{j} = (raw{j}-bgcmos);
        end
        if s.bias(j)
            raw{j} = raw{j}./bias_cell{j};
        end
        if s.postbin> 0
            raw{j} = imresize(raw{j},settings.postbin);
            [height,width]=size(raw{j});
        end
        if s.signal_foreground(j)
            normnuc = (raw{1}-prctile(raw{1}(:),1));
            normnuc(normnuc < 1) = 1;
            normnuc = normnuc/prctile(normnuc(:),99);
            normraw = (raw{j}-prctile(raw{j}(:),1));
            normraw(normraw < 1) = 1;
            normraw = normraw/prctile(normraw(:),99);
            sumraw = normnuc + normraw;
            %foreground{j} = threshmask_adapt(raw{j},s.blurradius);
            foreground{j} = threshmask(sumraw,s.blurradius);
            if sum(foreground{j}(:))/length(foreground{j}(:)) > .9 % account for large foregrounds
                foreground{j} = zeros(height, width);
            end
        end
    end
    
    %%% Clear unnecessary variables
    clear normnuc normraw sumraw
    
    %% Segment nuclei %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    segparams = struct('nucr', s.nucr,'blurradius', s.blurradius,'debrisarea', s.debrisarea,...
        'boulderarea',s.boulderarea, 'blobthreshold', s.blobthreshold);
    nuc_mask = segmentImage(raw{1},s.firstsegmethod, segparams);
    
    whole_mask = nuc_mask;
    nuc_mask = bwareaopen(nuc_mask,s.debrisarea);
    
    %%% adaptive bridging settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mask_props = regionprops(nuc_mask,'Area','Solidity');
    sizes = [mask_props(:).Area];
    solid = [mask_props(:).Solidity];
    medSize = prctile(sizes,50);
    gate = sizes >= medSize | solid < s.soliditythresh;
    
    %         %%% adaptive bridging settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %         mask_props = regionprops(nuc_mask,raw{1},'Area','Circularity','PixelIdxList');
    %         inten = [];
    %         varInt = [];
    %         varIm = stdfilt(raw{1},true(5));
    %         nuc_mask_erode = imerode(nuc_mask,strel('disk',5));
    %         varIm(~nuc_mask_erode) = NaN;
    %         for n = 1:length(mask_props)
    %             pix = mask_props(n).PixelIdxList;
    %             inten(n) = mean(raw{1}(pix));
    %             varInt(n) = nanmean(varIm(pix))/inten(n);
    %         end
    %         sizes = [mask_props(:).Area];
    %         circ = [mask_props(:).Circularity];
    %         medSize = median(sizes);
    %         medInten = median(inten);
    %         varThresh = .3;
    %         isMitotic = circ<.4 & sizes < medSize*1.5 & inten > medInten;
    %         isLarge = sizes >= medSize;
    %         isDim = inten < 200;
    %         gate = isLarge & ~isMitotic & ~isDim;
    
    %%% Clean up segmentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nuc_mask = segmentdeflections_bwboundaries_adaptive(nuc_mask,s.nucr,gate);
    %nuc_mask = segmentdeflections_bwboundaries(nuc_mask,s.nucr,s.debrisarea);
    nuc_mask = bwareaopen(nuc_mask,s.debrisarea);
    %
    %         mask_props = regionprops(nuc_mask,raw{1},'Area','PixelList');
    %         stillLarge = find([mask_props.Area] > medSize*4);
    %
    %         for n = 1:length(stillLarge)
    %             pixels = mask_props(stillLarge(n)).PixelList;
    %             minX = min(pixels(:,1));
    %             maxX = max(pixels(:,1));
    %             minY = min(pixels(:,2));
    %             maxY = max(pixels(:,2));
    %             image = raw{1}(minX:maxX,minY:maxY);
    %         end
    %
    
    nuc_mask = excludelargeandwarped_3(nuc_mask,s.boulderarea,s.soliditythresh);
    
    %%% remove border objects %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    antiborder_mask=imclearborder(nuc_mask);
    border_mask=nuc_mask-antiborder_mask;
    nuc_mask=logical(nuc_mask-border_mask);
    
    %%% Check segmentation
    if debug_mode
        keyboard;
    end
    %{
        %%% debugging: view images %%%%%%%%%%
        %Check nuclear segmentation
        extractmask = bwmorph(nuc_mask,'remove');
        tempframe = zeros(height,width,3);
        tempframe(:,:,1) = imadjust(mat2gray((raw{1})));
        tempframe(:,:,2) = extractmask;;
        %tempframe(:,:,3) = marker_mask;
        figure,imshow(tempframe);

        %Check foreground segmentation
        extractmask = bwmorph(foreground{2},'remove');
        tempframe = zeros(height,width,3);
        tempframe(:,:,1) = imadjust(mat2gray((raw{2})));
        tempframe(:,:,2) = extractmask;;
        %tempframe(:,:,3) = marker_mask;
        figure,imshow(tempframe);

        nuc_info = struct2cell(regionprops(nuc_mask,'Area')');
        nuc_area = squeeze(cell2mat(nuc_info(1,1,:)));
        hist(nuc_area,100);

        anti_mask = bwareaopen(nuc_mask,s.debrisarea);
        temp_mask = nuc_mask-anti_mask;
        extractmask = bwmorph(temp_mask,'remove');

        anti_mask = bwareaopen(nuc_mask,3500);
        extractmask = bwmorph(anti_mask,'remove');
    %}
    
    %% background subtract: masked %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for j = 1:length(s.signals)
        if s.signal_foreground(j)
            mask_total  =  whole_mask | foreground{j};
            nanmask{j} = imdilate(mask_total,strel('disk',s.nucr));
            nanmaskcyto{j} = imdilate(mask_total,strel('disk',s.nucr*2));
        else
            nanmask{j} = imdilate(whole_mask,strel('disk',s.nucr));
            nanmaskcyto{j} = imdilate(whole_mask,strel('disk',s.nucr*2));
        end
        
        %%% Blur signal
        if s.sigblur > 0
            blurIm = imfilter(raw{j},fspecial('disk',s.sigblur(j)),'symmetric');
        else
            blurIm  =  raw{j};
        end
        
        %%% Background subtract signals
        switch s.bgsubmethod{j}
            case 'global nuclear'
                real{j} = bgsubmasked_global_NR(blurIm,nanmask{j},1,s.compression,s.bgprctile(j));
                unblurred{j} = bgsubmasked_global_NR(raw{j},nanmask{j},1,s.compression,s.bgprctile(j));
            case 'global cyto'
                real{j} = bgsubmasked_global_2(blurIm,nanmaskcyto{j},1,s.compression,s.bgprctile(j));
                unblurred{j} = bgsubmasked_global_2(raw{j},nanmaskcyto{j},1,s.compression,s.bgprctile(j));
            case 'tophat'
                real{j} = imtophat(blurIm,strel('disk',s.nucr,0));
                unblurred{j} = imtophat(raw{j},strel('disk',s.nucr,0));
            case 'semi-local nuclear'
                real{j} = bgsubmasked_global_2(blurIm,nanmask{j},11,s.compression,s.bgprctile(j));
                unblurred{j} = bgsubmasked_global_2(raw{j},nanmask{j},11,s.compression,s.bgprctile(j));
            case 'none'
                real{j} = blurIm;
                unblurred{j} = raw{j};
        end
    end
    %%% Clear unnecessary variables
    clear raw mask_total foreground mask_total blurIm
    
    %% Extract features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [nuc_label,numcells] = bwlabel(nuc_mask);
    nuc_info = regionprops(nuc_mask,real{1},'Area','Centroid','MeanIntensity','PixelIdxList');
    nuc_area = [nuc_info.Area]';
    nuc_center = cell2mat({nuc_info.Centroid}');
    nuc_density = [nuc_info.MeanIntensity]';
    
    %%% calculate masses
    nuc_mass = nuc_density.*nuc_area;
    nanvec = ones(numcells,1)*NaN;
    
    %%% Clear unnecessary variables
    clear nuc_mask
    
    %%% Make and save mask
    extractmask = bwmorph(nuc_label,'remove'); % used for next time point
    if s.maskwrite
        imwrite(uint16(extractmask),fullfile(maskdir,[shot,'_',s.maskname,num2str(s.frame),'.tif']));
    end
    
    maxtrackedID=max(nuc_label(:));
    border_label=bwlabel(border_mask);
    border_label=border_label+maxtrackedID;
    border_label(border_label==maxtrackedID)=0;
    all_label=nuc_label+border_label;
    
    %% Measure features%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% Initialize measurement variables%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sigringmedian{1} = [];
    localbg{1} = [];
    for j = 1:length(s.signals)
        sigmean{j} = nanvec;
        sigmedian{j} = nanvec;
        if s.ringcalc(j)
            sigringmedian{j} = nanvec;
        end
        if s.localbgmeasure(j)
            localbg{j} = nanvec;
        end
    end
    
    %%% Cyto ring calculation
    if any(s.ringcalc)
        if (s.magnification == 10 && s.binsize == 1) || (s.magnification == 20 && s.binsize == 2)
            innerrad = 1; outerrad = 5;
        else
            innerrad = 2; outerrad = 10;
        end
        ring_label = getcytoring_thicken(all_label,innerrad,outerrad,real{2});
        ring_info = regionprops(ring_label,'PixelIdxList');
        
        %%% Check cytoring
        if debug_mode
            keyboard;
        end
        %{
            %%% debugging: view images %%%%%%%%%%
            temp_ring  =  ring_label > 0;
            extractmask = bwmorph(temp_ring,'remove');
            signal  =  real{1};
            intensities  =  signal.*temp_ring;
            tempframe = zeros(height,width,3);
            tempframe(:,:,1) = imadjust(mat2gray(real{1}));
            tempframe(:,:,2) = extractmask;;
            %tempframe(:,:,3) = marker_mask;
            figure,imshow(tempframe);
        %}
    end
    
    %% Measure for each signal
    for j = 1:length(s.signals)
        real_temp = real{j};
        real_masked = real_temp;
        real_masked(nanmask{j}) = NaN;
        
        %%% Puncta analysis
        if j == s.punctacalc
            p = punctaModule(unblurred{j}, nuc_label, nanvec, s.punctaThresh, s.varThresh,debug_mode, numlivecells, nuc_info,cellid);
        elseif s.punctacalc == 0
            p = createPunctaStruct();
        end
        
        %%% Measurements for each cell
        for cc = 1:numcells
            sigmean{j}(cc) = mean(real_temp(nuc_info(cc).PixelIdxList));
            sigmedian{j}(cc) = median(real_temp(nuc_info(cc).PixelIdxList));
            if s.localbgmeasure(j)
                localbg{j}(cc) = getlocalbg(real_masked, nuc_info(cc).Centroid,100, 25, .15);
            end
            if s.ringcalc(j)
                if ~(cc > numel(ring_info))
                    ringall = real_temp(ring_info(cc).PixelIdxList);
                    ringall(ringall>prctile(ringall,98)) = [];
                    ringforeground = ringall(ringall>s.ringthresh(j));
                    if numel(ringforeground)<50
                        ringforeground = ringall;
                    end
                    if numel(ringall)>= 50
                        sigringmedian{j}(cc) = nanmedian(ringforeground);
                    end
                end
            end
        end
    end
    
    %% Save data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tracedata(:,1,:) = ...
        [nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass,...
        cell2mat(sigmean),cell2mat(sigmedian),cell2mat(localbg),...
        cell2mat(sigringmedian), ...
        p.PCNAMean, p.filtMean, p.varMean, p.varStd, ...
        cell2mat(p.filtMaskedMean), cell2mat(p.filtMaskedArea), ...
        cell2mat(p.varMaskedMean), cell2mat(p.varMaskedArea)];
    jitters = [0 0];
    
    
    %%% save data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    save(fullfile(s.savepath,['tracedata_',shot,'.mat']),'tracedata','jitters');
    names = output_names_puncta(s.signals,s.localbgmeasure, s.ringcalc, s.punctacalc, s.punctaThresh, s.varThresh);
    save(fullfile(s.savepath,'settings_live.mat'),'s','names'); %saves all the settings to the Data folder. will overwrite with the most recent
    
    toc(timetotal);
end
end


function p = punctaModule(punctaSig, nuc_label, nanvec, punctaThresh, varThresh, debug_mode, numlivecells, nuc_info,cellid)
p = createPunctaStruct();

%% Initialize PCNA storage
p.PCNAMean = nanvec;
p.filtMean = nanvec;
p.varMean = nanvec;
p.varStd = nanvec;
for t = 1:length(punctaThresh)
    p.filtMaskedMean{t} = nanvec;
    p.filtMaskedArea{t} = nanvec;
end
for t = 1:length(varThresh)
    p.varMaskedMean{t} = nanvec;
    p.varMaskedArea{t} = nanvec;
end

%% Preprocess PCNA image
punctaBlur = imgaussfilt(punctaSig, 8); %blurred to remove nucleoli
%punctaSemiBlur = imgaussfilt(punctaSig, 1);
punctaSemiBlur = punctaSig;

%Masking to remove nucleoli and border
nuclear_mask = nuc_label > 0;
punctaDiff = punctaSig - punctaBlur;
erode_nuc = imerode(nuclear_mask, strel('disk',7,0));
nucleolus_mask = (punctaDiff < -10);
nucleolus_mask = imdilate(nucleolus_mask, strel('disk',1,0));
analysis_mask = erode_nuc & ~nucleolus_mask;

%Image transformations quantification
punctaFilt = imtophat(punctaSig, strel('disk',2,0));
punctaVar = stdfilt(punctaSemiBlur,ones(3));

%Mask transformed images
nanVar = punctaVar;
nanVar(~analysis_mask) = NaN;

%% Check masking/puncta
if debug_mode
    keyboard;
end
%{
      extractmask = bwmorph(analysis_mask,'remove');
      tempframe=imadjust(mat2gray(punctaVar));
      tempframe(:,:,2)=(extractmask);
      tempframe(:,:,3)=0;
      figure,imshow(tempframe)

      imagesc(punctaSig)
      imtool(punctaFilt)
      imtool(punctaVar.*analysis_mask)
%}

%% Measurements for each cell
for n = 1:numlivecells
    cc = cellid(n);
    p.PCNAMean(n) = mean(punctaSemiBlur(nuc_info(cc).PixelIdxList));
    p.filtMean(n) = mean(punctaFilt(nuc_info(cc).PixelIdxList));
    p.varMean(n) = nanmean(nanVar(nuc_info(cc).PixelIdxList));
    p.varStd(n) = nanstd(nanVar(nuc_info(cc).PixelIdxList));
end

%% Measure threshold transformed images
for t = 1:length(punctaThresh)
    filtMask = punctaFilt > punctaThresh(t) & nuclear_mask;
    punctaFiltMasked = punctaFilt * filtMask{t};
    for n = 1:numlivecells
        p.filtMaskedMean{t}(n) = nanmean(punctaFiltMasked(nuc_info(cc).PixelIdxList));
        p.filtMaskedArea{t}(n) = sum(filtMask(nuc_info(cc).PixelIdxList));
    end
end
for t = 1:length(varThresh)
    varMask = punctaVar > varThresh(t) & analysis_mask;
    punctaVarMasked = punctaVar * varMask{t};
    for n = 1:numlivecells
        p.varMaskedMean{t}(n) = nanmean(punctaVarMasked(nuc_info(cc).PixelIdxList));
        p.varMaskedArea{t}(n) = sum(varMask(nuc_info(cc).PixelIdxList));
    end
end
end


function p = createPunctaStruct()

% PCNA variables
p.PCNAMean = [];
p.filtMean = [];
p.varMean = [];
p.varStd = [];
p.filtMaskedMean{1} = [];
p.filtMaskedArea{1} = [];
p.varMaskedMean{1} = [];
p.varMaskedArea{1} = [];

end


function mask = segmentImage(image, segMethod, s)
switch segMethod
    case 'concavity'
        %imblur = imfilter(image,fspecial('gaussian',s.s.blurradius),'symmetric');
        imblur = imfilter(image,fspecial('gaussian',5),'symmetric');
        mask = logical(ThreshImage(imblur));
        mask = (imfill(mask,'holes'));
        mask_noholes = (imfill(mask,'holes'));
        holes = mask_noholes & ~mask;
        holes = bwareaopen(holes,round(s.nucr^2*pi));
        mask = mask_noholes & ~holes;
        
    case 'log'
        mask = blobdetector_4(log(image),s.nucr,s.blobthreshold,s.debrisarea);
        %         mask = activecontour(image,mask,10);
    case 'single adapt'
        mask = threshmask_adapt(image,s.blurradius);
        mask = markershed_filter(mask,round(s.nucr*2/3),6);
    case 'double adapt'
        mask = threshmask_adapt(image,s.blurradius);
        mask = markershed_filter(mask,round(s.nucr*2/3),6);
        mask = secondthresh(image,s.blurradius,mask,s.boulderarea*2);
    case 'single'
        mask = threshmask(image,s.blurradius);
        mask = markershed_filter(mask,round(s.nucr*2/3),6);
    case 'double'
        mask = threshmask(image,s.blurradius);
        mask = markershed_filter(mask,round(s.nucr*2/3),6);
        mask = secondthresh(image,s.blurradius,mask,s.boulderarea*2);
end
end