function Timelapse(row,col,site,debug_mode)
row = 4;
col = 2;
site = 1;
debug_mode = 0;
timetotal = tic;

%% SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Paths
s.experiment_name='C179-live';
s.image_drive = 'I:\4TB8\Data\';
s.savepath=['F:\Data\C-Cdt1\',s.experiment_name,'\Data\'];
s.imagepath=[s.image_drive,s.experiment_name,'\Raw\'];
s.biaspath=[s.image_drive,s.experiment_name,'\Bias\'];
s.maskpath=[s.image_drive,s.experiment_name,'\Mask\'];
s.bgcmospath='F:\MATLAB\BGimages\IX_XL\cmosoffset_bin1.mat';

%%% General parameters
s.startFrame=1;   %Frame to start analyzing from
s.endFrame=93;  %Frame to stop analyzing
s.magnification=10; %10=10x or 20=20x
s.binsize=1; %1=bin 1 or 2=bin 2
s.signals={'CFP_','YFP_','RFP_'};
s.maskwrite = 1; %1 = save an image with the nucmask; 0 = dont save an image with the nucmask
s.maskname = 'nucedge_';

%%% Quantification parameters
s.bgcmoscorrection = 1;
s.bias = [1 1 1];
s.signal_foreground = [0 0 0];
s.bgsubmethod = {'global nuclear','global nuclear','global nuclear'}; %Options:'global nuclear','global cyto','tophat','semi-local nuclear', 'none'
s.compression = 4;
s.bgprctile = [25 25 25];
s.sigblur = [3 3 3];
s.localbgmeasure = [0 0 0];
s.ringcalc = [0 1 0];
s.ringthresh = [0 50 50];
s.punctacalc = 0;
s.punctaThresh = [125 150 200];
s.varThresh = [75 100 125];

%%% Segmentation parameters
s.firstsegmethod = 'log'; %Options: 'concavity', 'log', 'single adapt', 'double adapt', 'single', 'double'
s.secondsegmethod = 'log'; %Options: 'concavity', 'log', 'single adapt', 'double adapt', 'single', 'double'
s.nucr = 12; %10x bin1 = 12 and 20x bin2 = 12
s.blurradius  =  3; %10x: 3
s.soliditythresh = 0.8;
s.debrisarea = 150; % 10x = 100 % 150 if no mitosis
s.boulderarea = 1500; % 10x = 1500
blobthreshold = -0.02;

%%% Tracking parameters
s.maxjump = s.nucr*3;
s.masschangethreshold = 0.30;
s.areachangethreshold = 0.60;
s.daughtervariance = 0.10;

%% TRACKING ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
shot = [num2str(row),'_',num2str(col),'_',num2str(site)];
if ~exist([s.savepath,'tracedata_',shot,'.mat'],'file')
    
    %% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Initialize path variables
    rawdir = [s.imagepath,shot,'\',shot,'_'];
    maskdir = [s.maskpath,shot,'\',shot,'_'];
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
    
    %%% Initialize tracking variables
    frames = s.startFrame:s.endFrame;
    totalframes = numel(frames);
    badframes = ones(s.endFrame,1)*NaN;
    if s.startFrame>1
        badframes(1:s.startFrame-1) = 0;
    end
    jitters = zeros(s.endFrame,2);
    blocksize = 10000;
    maxcellnum = blocksize;
    tracedata = ones(maxcellnum,s.endFrame,parameternum)*NaN;
    tracking = ones(maxcellnum,5)*NaN;
    [firstgoodindex,blurthreshhigh,blurthreshlow,numthresh,badframes,height,width] = ...
        timelapsesetup_4(rawdir,nucname,frames,s.nucr,blobthreshold,s.debrisarea,badframes,s.maskwrite);
    regheight = 1:0.5*height; regwidth = 1:0.5*width;
    segparams = struct('nucr', s.nucr,'blurradius', s.blurradius,'debrisarea', s.debrisarea,...
        'boulderarea',s.boulderarea, 'blobthreshold', blobthreshold);
    trackparams = {s.nucr,s.maxjump,s.debrisarea,s.masschangethreshold,s.areachangethreshold,s.daughtervariance};
    
    %%% Load bgcmos/s.bias
    load(s.bgcmospath,'cmosoffset');
    bgcmos  =  cmosoffset;
    for i = 1:length(s.signals)
        if s.bias(i)
            load([s.biaspath,s.signals{i},num2str(site),'.mat']);
            bias_cell{i} = bias;
        end
    end
    
    %% Imaging processing for each frame %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = firstgoodindex:totalframes
        f = frames(i);
        fprintf('frame %0.0f\n',f);
        timeframe = tic;
        
        %%% Load images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for j = 1:length(s.signals)
            raw{j} = double(imread([rawdir,s.signals{j},num2str(f),'.tif']));
            if s.bgcmoscorrection
                raw{j} = (raw{j}-bgcmos);
            end
            if s.bias(j)
                raw{j} = raw{j}./bias_cell{j};
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
        if i == firstgoodindex
            nuc_mask = segmentImage(raw{1},s.firstsegmethod, segparams);
        else
            if isequal(s.secondsegmethod, 'apriori')
                nuc_mask = threshmask(raw{1});
            else
                nuc_mask = segmentImage(raw{1},s.secondsegmethod, segparams);
            end
            %%% Calculate jitter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            lastgoodframe = find(badframes == 0,1,'last');
            [reljitx,reljity] = registerimages(imfill(extractmask(regheight,regwidth),'holes'),nuc_mask(regheight,regwidth));
            jitters(f,:) = jitters(lastgoodframe,:)+[reljitx,reljity];
            if isequal(s.secondsegmethod, 'apriori')
                [nuc_mask,marker_mask] = apriori_markermask(nuc_mask,nuc_center,jitters(f,:));
            end
        end
        %nuc_mask = imopen(nuc_mask,strel('disk',3));
        whole_mask = nuc_mask;
        
        %%% Clean up segmentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        nuc_mask = bwareaopen(nuc_mask,s.debrisarea);
        nuc_mask = segmentdeflections_bwboundaries(nuc_mask,s.nucr,s.debrisarea);
        nuc_mask = bwareaopen(nuc_mask,s.debrisarea);
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
        nuc_info = struct2cell(regionprops(nuc_mask,real{1},'Area','Centroid','MeanIntensity')');
        nuc_area = squeeze(cell2mat(nuc_info(1,1,:)));
        
        %%% detect bad frame
        mednuc = median(nuc_area);
        if i>firstgoodindex+1 && (numcells<numthresh || mednuc>blurthreshhigh || mednuc<blurthreshlow)
            fprintf('badframe: frame %0.0f\n',f);
            badframes(f) = 1;
            extractmask = bwmorph(nuc_mask,'remove');
            if s.maskwrite
                imwrite(uint16(extractmask),[maskdir,s.maskname,num2str(f),'.tif']);
            end
            continue;
        end
        blurthreshhigh = 1.2*mednuc;
        blurthreshlow = 0.8*mednuc;
        numthresh = 0.5*numcells;
        nuc_center = squeeze(cell2mat(nuc_info(2,1,:)))';
        nuc_density = squeeze(cell2mat(nuc_info(3,1,:)));
        
        %%% calculate masses
        nuc_mass = nuc_density.*nuc_area;
        curdata = [nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass];
        
        %%% Clear unnecessary variables
        clear nuc_mask
        %% Track cells %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if i>firstgoodindex
            nuc_center(:,1) = nuc_center(:,1)+jitters(f,1);
            nuc_center(:,2) = nuc_center(:,2)+jitters(f,2);
            %%% temporarily store values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            curdata = [nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass];
            debugpackage = {extractmask,jitters(lastgoodframe,:),[reljitx,reljity]};
            %%% track & correct merges (update centers, masses and labels) %%%%
            [tracedata,curdata,tracking,nuc_label] = ...
                adaptivetrack_9(f,lastgoodframe,f,tracedata,curdata,tracking,real{1},nuc_label,jitters(f,:),trackparams,debugpackage);
            badframes(f) = 0;
            
            %%% Check tracking
            if debug_mode
                keyboard;
            end
            %{
            matchinglabel=prev_label.*nuc_label;
            tempframe=zeros(height,width,3);
            tempframe(:,:,1)=prev_label>0;
            tempframe(:,:,2)=nuc_label>0;
            tempframe(:,:,3)=matchinglabel;
            figure,imshow(tempframe);
            %}
        end
        
        %%% Make and save mask
        extractmask = bwmorph(nuc_label,'remove'); % used for next time point
        if s.maskwrite
            imwrite(uint16(extractmask),[maskdir,s.maskname,num2str(f),'.tif']);
        end
        prev_label = nuc_label;
        
        %%% Reformat features
        cellid = find(~isnan(curdata(:,1)));
        numlivecells = numel(cellid);
        curdata = curdata(cellid,:);
        nuc_center = curdata(:,[1 2]);
        nuc_area = curdata(:,3);
        nuc_mass = curdata(:,4);
        nuc_info = regionprops(nuc_label,'PixelIdxList','Centroid', 'BoundingBox');
        nanvec = ones(numlivecells,1)*NaN;
        
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
            for n = 1:numlivecells
                cc = cellid(n);
                sigmean{j}(n) = mean(real_temp(nuc_info(cc).PixelIdxList));
                sigmedian{j}(n) = median(real_temp(nuc_info(cc).PixelIdxList));
                if s.localbgmeasure(j)
                    localbg{j}(n) = getlocalbg(real_masked, nuc_info(cc).Centroid,100, 25, .15);
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
                            sigringmedian{j}(n) = nanmedian(ringforeground);
                        end
                    end
                end
            end
        end
        
        %% Save data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tracedata(cellid,f,:) = ...
            [nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass,...
            cell2mat(sigmean),cell2mat(sigmedian),cell2mat(localbg),...
            cell2mat(sigringmedian), ...
            p.PCNAMean, p.filtMean, p.varMean, p.varStd, ...
            cell2mat(p.filtMaskedMean), cell2mat(p.filtMaskedArea), ...
            cell2mat(p.varMaskedMean), cell2mat(p.varMaskedArea)];
        if maxcellnum-max(cellid)<blocksize
            tempdata = ones(blocksize,s.endFrame,parameternum)*NaN;
            temptrack = ones(blocksize,5)*NaN;
            tracedata = [tracedata;tempdata];
            tracking = [tracking;temptrack];
            maxcellnum = maxcellnum+blocksize;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        toc(timeframe);
    end
    [tracedata,genealogy,jitters] = postprocessing_nolinking(tracedata,cellid,jitters,badframes,tracking,maxcellnum,s.nucr);
    %%% save data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    save([s.savepath,'tracedata_',shot,'.mat'],'tracedata','genealogy','jitters');
    names = output_names_puncta(s.signals,s.localbgmeasure, s.ringcalc, s.punctacalc, s.punctaThresh, s.varThresh);
    save([s.savepath,'settings_live.mat'],'s','names'); %saves all the settings to the Data folder. will overwrite with the most recent
    
    toc(timetotal);
    clear all;
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

%Threshold transformed images
for t = 1:length(punctaThresh)
    filtMask{t} = punctaFilt > punctaThresh(t) & nuclear_mask;
    punctaFiltMasked{t} = punctaFilt * filtMask{t};
end
for t = 1:length(varThresh)
    varMask{t} = punctaVar > varThresh(t) & analysis_mask;
    punctaVarMasked{t} = punctaVar * varMask{t};
end

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
    for t = 1:length(punctaThresh)
        p.filtMaskedMean{t}(n) = nanmean(punctaFiltMasked{t}(nuc_info(cc).PixelIdxList));
        p.filtMaskedArea{t}(n) = sum(filtMask{t}(nuc_info(cc).PixelIdxList));
    end
    for t = 1:length(varThresh)
        p.varMaskedMean{t}(n) = nanmean(punctaVarMasked{t}(nuc_info(cc).PixelIdxList));
        p.varMaskedArea{t}(n) = sum(varMask{t}(nuc_info(cc).PixelIdxList));
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