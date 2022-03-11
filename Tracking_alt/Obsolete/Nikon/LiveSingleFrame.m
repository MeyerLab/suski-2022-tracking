function LiveSingleFrame(row,col,site,debug_mode)
row = 3;
col = 2;
site = 1;
debug_mode = 0;
timetotal = tic;

%% SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Paths
experiment_name='C183-live1frame';
image_drive = 'I:\4TB8\Data\';
savepath=['F:\Data\C-Cdt1\',experiment_name,'\Data\'];
imagepath=[image_drive,experiment_name,'\Raw\'];
biaspath=[image_drive,experiment_name,'\Bias\'];
maskpath=[image_drive,experiment_name,'\Mask\'];
bgcmospath='F:\MATLAB\BGimages\IX_XL\cmosoffset_bin1.mat';

%%% General parameters
startFrame=1;   %Frame to start analyzing from
endFrame=1;  %Frame to stop analyzing
magnification=10; %10=10x or 20=20x
binsize=1; %1=bin 1 or 2=bin 2
signals={'CFP_','YFP_'};
maskwrite = 1; %1 = save an image with the nucmask; 0 = dont save an image with the nucmask
maskname = 'nucedge_';

%%% Quantification parameters
bgcmoscorrection = 1;
bias = [1 1];
signal_foreground = [0 0];
bgsubmethod = {'global nuclear','constant YFP'}; %Options:'global nuclear','global cyto','tophat','semi-local nuclear', 'none'
compression = 4;
sigblur = [3 3];
localbgmeasure = [0 0];
ringcalc = [0 1];
ringthresh = [0 5 ];
punctacalc = 0;
punctaThresh = [175 200 225];
varThresh = [75 100 125];

%%% Segmentation parameters
firstsegmethod = 'concavity'; %Options: 'concavity', 'log', 'single adapt', 'double adapt', 'single', 'double'
nucr = 12; %10x bin1 = 12 and 20x bin2 = 12
blurradius  =  3; %10x: 3
soliditythresh = 0.8;
debrisarea = 100; % 10x = 100
boulderarea = 1500; % 10x = 1500
blobthreshold = -0.02;

%%% Tracking parameters
maxjump = nucr*3;
masschangethreshold = 0.30;
areachangethreshold = 0.60;
daughtervariance = 0.10;

%% TRACKING ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
shot = [num2str(row),'_',num2str(col),'_',num2str(site)];
if ~exist([savepath,'tracedata_',shot,'.mat'],'file')

    %% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Initialize path variables
    rawdir = [imagepath,shot,'\',shot,'_'];
    maskdir = [maskpath,shot,'\',shot,'_'];
    nucname = signals{1};
    parameternum = 4+2*length(signals)+sum(ringcalc)+sum(localbgmeasure);
    if punctacalc
        parameternum = parameternum + 4 + 2*length(punctaThresh) + 2*length(varThresh);
    end
    
    if ~exist(savepath,'dir')
        mkdir(savepath);
    end
    if ~exist(maskdir,'dir') && maskwrite
        mkdir(maskdir);
    end
    
    %%% Load bgcmos/bias
    load(bgcmospath,'cmosoffset');
    bgcmos  =  cmosoffset;
        [height,width]=size(bgcmos);

    for i = 1:length(signals)
        if bias(i)
            load([biaspath,signals{i},num2str(site),'.mat']);
            bias_cell{i} = bias;
        end
    end
    
    %% Imaging processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%% Load images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1:length(signals)
        raw{i} = double(imread([rawdir,signals{i},num2str(1),'.tif']));
        if bgcmoscorrection
            raw{i} = (raw{i}-bgcmos);
        end
        if bias(i)
            raw{i} = raw{i}./bias_cell{i};
        end
        if signal_foreground(i)
            %foreground{j} = threshmask_adapt(raw{j},blurradius);
            foreground{i} = threshmask(raw{i},blurradius);
            if sum(foreground{i}(:))/length(foreground{i}(:)) > .9 % account for large foregrounds
                foreground{i} = zeros(height, width);
            end
        end
    end
    
    %%% Segment nuclei %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    segparams = struct('nucr', nucr,'blurradius', blurradius,'debrisarea', debrisarea,...
        'boulderarea',boulderarea, 'blobthreshold', blobthreshold);
    nuc_mask = segmentImage(raw{1},firstsegmethod, segparams);
    whole_mask = nuc_mask;
    
    %%% Clean up segmentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nuc_mask = bwareaopen(nuc_mask,debrisarea);
    nuc_mask = segmentdeflections_bwboundaries(nuc_mask,nucr,debrisarea);
    nuc_mask = excludelargeandwarped_3(nuc_mask,boulderarea,soliditythresh);
    nuc_mask = imclearborder(nuc_mask);
    
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
        tempframe(:,:,1) = imadjust(mat2gray((raw{1})));
        tempframe(:,:,2) = extractmask;;
        %tempframe(:,:,3) = marker_mask;
        figure,imshow(tempframe);

        nuc_info = struct2cell(regionprops(nuc_mask,'Area')');
        nuc_area = squeeze(cell2mat(nuc_info(1,1,:)));
        hist(nuc_area,100);

        anti_mask = bwareaopen(nuc_mask,debrisarea);
        temp_mask = nuc_mask-anti_mask;
        extractmask = bwmorph(temp_mask,'remove');

        anti_mask = bwareaopen(nuc_mask,3500);
        extractmask = bwmorph(anti_mask,'remove');
    %}
    
    %%% background subtract: masked %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1:length(signals)
        if signal_foreground(i)
            mask_total  =  whole_mask | foreground{i};
            nanmask{i} = imdilate(mask_total,strel('disk',nucr));
            nanmaskcyto{i} = imdilate(mask_total,strel('disk',nucr*2));
        else
            nanmask{i} = imdilate(whole_mask,strel('disk',nucr));
            nanmaskcyto{i} = imdilate(whole_mask,strel('disk',nucr*2));
        end
        
        %%% Blur signal
        if sigblur > 0
            blurIm = imfilter(raw{i},fspecial('disk',sigblur(i)),'symmetric');
        else
            blurIm  =  raw{i};
        end
        
        %%% Background subtract signals
        switch bgsubmethod{i}
            case 'global nuclear'
                real{i} = bgsubmasked_global_NR(blurIm,nanmask{i},1,compression,25);
                unblurred{i} = bgsubmasked_global_NR(raw{i},nanmask{i},1,compression,25);
            case 'global cyto'
                real{i} = bgsubmasked_global_2(blurIm,nanmaskcyto{i},1,compression,25);
                unblurred{i} = bgsubmasked_global_2(raw{i},nanmaskcyto{i},1,compression,25);
            case 'tophat'
                real{i} = imtophat(blurIm,strel('disk',nucr,0));
                unblurred{i} = imtophat(raw{i},strel('disk',nucr,0));
            case 'semi-local nuclear'
                real{i} = bgsubmasked_global_2(blurIm,nanmask{i},11,compression,25);
                unblurred{i} = bgsubmasked_global_2(raw{i},nanmask{i},11,compression,25);
            case 'none'
                real{i} = blurIm;
                unblurred{i} = raw{i};
            case 'constant YFP'
                real{i} = blurIm - 200;
                unblurred{i} = raw{i} - 200;
        end
    end
    
    %%% Extract features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [nuc_label,numcells] = bwlabel(nuc_mask);
    nuc_info = struct2cell(regionprops(nuc_mask,real{1},'Area','Centroid','MeanIntensity'));
    nuc_area = (cell2mat(nuc_info(1,:)))';
    nuc_center = (cell2mat(nuc_info(2,:)'));
    nuc_density = (cell2mat(nuc_info(3,:)))';
    
    %%% calculate masses
    nuc_mass = nuc_density.*nuc_area;
   
    %%% Save mask
    extractmask = bwmorph(nuc_label,'remove');
    if maskwrite
        imwrite(uint16(extractmask),[maskdir,maskname,num2str(1),'.tif']);
    end
    
    nuc_info = regionprops(nuc_label,'PixelIdxList','Centroid');
    nanvec = ones(numcells,1)*NaN;
    
    %% Measure features%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% Initialize measurement variables%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sigringmedian{1} = [];
    localbg{1} = [];
    
    % PCNA variables
    PCNAMean = [];
    filtMean = [];
    varMean = [];
    varStd = [];
    filtMaskedMean{1} = [];
    filtMaskedArea{1} = [];
    varMaskedMean{1} = [];
    varMaskedArea{1} = [];
    
    for i = 1:length(signals)
        sigmean{i} = nanvec;
        sigmedian{i} = nanvec;
        if ringcalc(i)
            sigringmedian{i} = nanvec;
        end
        if localbgmeasure(i)
            localbg{i} = nanvec;
        end
    end
    
    if any(ringcalc)
        if (magnification == 10 && binsize == 1) || (magnification == 20 && binsize == 2)
            innerrad = 1; outerrad = 5;
        else
            innerrad = 2; outerrad = 10;
        end
        ring_label = getcytoring_thicken(nuc_label,innerrad,outerrad,real{2});
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
            tempframe(:,:,1) = imadjust(mat2gray(real{2}));
            tempframe(:,:,2) = extractmask;;
            %tempframe(:,:,3) = marker_mask;
            figure,imshow(tempframe);
        %}
    end
    
    %%% Measure for each signal
    for i = 1:length(signals)
        real_temp = real{i};
        real_masked = real_temp;
        real_masked(nanmask{i}) = NaN;
        
        %% Puncta analysis
        if i == punctacalc
            %Preprocess PCNA image
            punctaSig = unblurred{punctacalc};
            punctaBlur = imgaussfilt(punctaSig, 8); %blurred to remove nucleoli
            %punctaSemiBlur = imgaussfilt(punctaSig, 1);
            punctaSemiBlur = punctaSig;
            
            %Masking to remove nucleoli and border
            nuclear_mask = nuc_label > 0;
            punctaDiff = real_temp - punctaBlur;
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
            
            %Threshold transformed images
            for t = 1:length(punctaThresh)
                filtMask{t} = punctaFilt > punctaThresh(t) & nuclear_mask;
                punctaFiltMasked{t} = punctaFilt * filtMask{t};
            end
            for t = 1:length(varThresh)
                varMask{t} = punctaVar > varThresh(t) & analysis_mask;
                punctaVarMasked{t} = punctaVar * varMask{t};
            end
            
            %Initialize PCNA storage
            PCNAMean = nanvec;
            filtMean = nanvec;
            varMean = nanvec;
            varStd = nanvec;
            for t = 1:length(punctaThresh)
                filtMaskedMean{t} = nanvec;
                filtMaskedArea{t} = nanvec;
            end
            for t = 1:length(varThresh)
                varMaskedMean{t} = nanvec;
                varMaskedArea{t} = nanvec;
            end
            
            % Check masking/puncta
            if debug_mode
                keyboard;
            end
            %{
               extractmask = bwmorph(analysis_mask,'remove');
               tempframe=imadjust(mat2gray(punctaVar));
               tempframe(:,:,2)=(extractmask);
               tempframe(:,:,3)=0;
               figure,imshow(tempframe)
            %}
        end
        
        %% Measurements for each cell
        for cc = 1:numcells
            sigmean{i}(cc) = mean(real_temp(nuc_info(cc).PixelIdxList));
            sigmedian{i}(cc) = median(real_temp(nuc_info(cc).PixelIdxList));
            if localbgmeasure(i)
                localbg{i}(cc) = getlocalbg(real_masked, nuc_info(cc).Centroid,100, 25, .15);
            end
            if ringcalc(i)
                if ~(cc > numel(ring_info))
                    ringall = real_temp(ring_info(cc).PixelIdxList);
                    ringall(ringall>prctile(ringall,98)) = [];
                    ringforeground = ringall(ringall>ringthresh(i));
                    if numel(ringforeground)<50
                        ringforeground = ringall;
                    end
                    if numel(ringall)>= 50
                        sigringmedian{i}(cc) = nanmedian(ringforeground);
                    end
                end
            end
            if i == punctacalc
                PCNAMean(cc) = mean(punctaSemiBlur(nuc_info(cc).PixelIdxList));
                filtMean(cc) = mean(punctaFilt(nuc_info(cc).PixelIdxList));
                varMean(cc) = nanmean(nanVar(nuc_info(cc).PixelIdxList));
                varStd(cc) = nanstd(nanVar(nuc_info(cc).PixelIdxList));
                for t = 1:length(punctaThresh)
                    filtMaskedMean{t}(cc) = nanmean(punctaFiltMasked{t}(nuc_info(cc).PixelIdxList));
                    filtMaskedArea{t}(cc) = sum(filtMask{t}(nuc_info(cc).PixelIdxList));
                end
                for t = 1:length(varThresh)
                    varMaskedMean{t}(cc) = nanmean(punctaVarMasked{t}(nuc_info(cc).PixelIdxList));
                    varMaskedArea{t}(cc) = sum(varMask{t}(nuc_info(cc).PixelIdxList));
                end
            end
        end
    end
    
    %% Save data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tracedata(:,1,:) = ...
        [nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass,...
        cell2mat(sigmean),cell2mat(sigmedian),cell2mat(localbg),...
        cell2mat(sigringmedian), PCNAMean, ...
        filtMean, varMean, varStd, cell2mat(filtMaskedMean), cell2mat(filtMaskedArea), ...
        cell2mat(varMaskedMean), cell2mat(varMaskedArea)];
   jitters = [0 0];
        
    %%% save data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    save([savepath,'tracedata_',shot,'.mat'],'tracedata','jitters');
    names = output_names_puncta(signals,localbgmeasure, ringcalc, punctacalc, punctaThresh, varThresh);
    save([savepath,'settings_live.mat'],'names'); %saves all the settings to the Data folder. will overwrite with the most recent
    
    toc(timetotal);
    clear all;
end
end

function mask = segmentImage(image, segMethod, s)
switch segMethod
    case 'concavity'
        %imblur = imfilter(image,fspecial('gaussian',s.blurradius),'symmetric');
        imblur = imfilter(image,fspecial('gaussian',5),'symmetric');
        mask = ThreshImage(imblur);
        mask = logical(imfill(mask,'holes'));
    case 'log'
        mask = blobdetector_4(log(image),s.nucr,s.blobthreshold,s.debrisarea);
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