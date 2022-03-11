function Timelapse(s,row,col,site,debug_mode,varargin)
%Timelapse timelapse tracking analysis. Inputs parameter struct s, row,
%          col, site. debug_mode for setting up analysis, varargin is
%          ignored. Compatible with IX micro and nikon naming schemes.
%%% REQUIRED SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% Paths
% s.image_drive = '*';    % Location of images base folder
% s.savepath=fullfile('*\Data'); % Output location of mat files
% s.imagepath=fullfile(s.image_drive,'Raw');                       % All sites in same folder for nikon
% s.biaspath=fullfile(s.image_drive,'Bias');                       % Location of bias files
% s.maskpath=fullfile(s.image_drive,'Mask');                       % Future location of mask saving
% s.bgcmospath='';
% % Formatting code for nikon files, 1$ is timepoint, 2$ is well, 3$ is point
% s.formatCode = 'Time%1$05d_Well%2$s_*.tiff'; %'Time%1$05d_Well%2$s_Point%2$s_%3$04d_*.tiff' if multipoint
% s.microscope = 'nikon';
% s.nucChannel = 4;

% %%% General parameters
% s.startFrame = 1; 
% s.endFrame = getMaxFrame(fullfile(s.imagepath,'2_2_1'),'Time(\d+)_');
% s.magnification=20;               % 10=10x or 20=20x
% s.binsize=2;                      % 1=bin 1 or 2=bin 2
% s.postbin = 0;             % 0 for no software bin, number for scaling factor (e.g. 0.5 for 2x2 bin);
% s.signals={'CFP_','YFP_','RFP_','IRFP_'};
% s.maskwrite = 1;                  % 1 = save an image with the nucmask; 0 = dont save an image with the nucmask
% s.maskname = 'nucedge_';
% s.register = 0;                   % Register images between timepoints for all
% s.register_exception = [1:5];        % If don't register images, manually put in sites to register
% 
% %%% Quantification parameters
% s.bgcmoscorrection = 1;
% s.bias = [0 0 0 0];
% s.signal_foreground = [0 0 0 0];
% s.bgsubmethod = {'global nuclear','global cyto','global cyto','global nuclear'}; 
%                                   % Options:'global nuclear','global cyto','tophat','semi-local nuclear', 'none'
% s.compression = 4;                % For bg subtraction
% s.bgprctile = [25 25 25 25];         % Set bg as perctile for each channel
% s.sigblur = [3 3 0 3]; 
% s.localbgmeasure = [0 0 0 0];       % Measure local background in box around cell
% s.ringcalc = [0 1 1 0];
% s.ringthresh = [0 0 0 0];         % Threshold for foreground ring pixels
% s.punctacalc = 0;                 % Analyze puncta in nucleus
% s.punctaThresh = [125 150 200];
% s.varThresh = [75 100 125];
% 
% %%% Segmentation parameters
% s.firstsegmethod = 'log contour';  % Options: 'concavity', 'log', 'single adapt', 'double adapt', 'single', 'double', 'log contour'
% s.secondsegmethod = 'log contour'; % Options: 'concavity', 'log', 'single adapt', 'double adapt', 'single', 'double','log contour'
% s.nucr = 12;                     % 10x bin1 = 12 and 20x bin2 = 12
% s.blurradius  =  3;              % 10x: 3
% s.soliditythresh = 0.8;          % Minimum solidity of nucleus allowed
% s.debrisarea = 50;              % 10x = 100 % 150 if no mitosis
% s.boulderarea = 1500;            % 10x = 1500
% s.blobthreshold = -0.02;         % For blobdector 'log' segmentation
% s.split_mult = 1;
% %%% Tracking parameters
% s.maxjump = s.nucr*3;
% s.masschangethreshold = 0.30;
% s.areachangethreshold = 0.60;
% s.daughtervariance = 0.10;
% 

%% TRACKING ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    nucname = s.signals{s.nucChannel};

    parameternum = 4+2*length(s.signals)+sum(s.ringcalc)+sum(s.localbgmeasure);
    if s.punctacalc
        parameternum = parameternum + 4 + 3*length(s.punctaThresh) + 3*length(s.varThresh);
    end
    
    if ~exist(s.savepath,'dir')
        mkdir(s.savepath);
    end
    if ~exist(maskdir,'dir') && s.maskwrite
        mkdir(maskdir);
    end
    
    %%% Initialize tracking variables
    frames = s.startFrame:s.endFrame(end);
    totalframes = numel(frames);
    badframes = ones(s.endFrame,1)*NaN;
    if s.startFrame>1
        badframes(1:s.startFrame-1) = 0;
    end
    jitters = zeros(s.endFrame(end),2);
    blocksize = 10000;
    maxcellnum = blocksize;
    tracedata = ones(maxcellnum,s.endFrame,parameternum)*NaN;
    tracking = ones(maxcellnum,5)*NaN;
    switch s.microscope
        case 'nikon'
            [firstgoodindex,blurthreshhigh,blurthreshlow,numthresh,badframes,height,width] = ...
                timelapsesetup_nikon_channel(rawdir,row,col,site, s.formatCode,frames,s.nucr,s.blobthreshold,s.debrisarea,badframes,s.maskwrite,s.nucChannel);
        case 'IXM'
            [firstgoodindex,blurthreshhigh,blurthreshlow,numthresh,badframes,height,width] = ...
                timelapsesetup_4(rawdir,nucname,frames,s.nucr,s.blobthreshold,s.debrisarea,badframes,s.maskwrite);
    end
    regheight = 1:0.5*height; regwidth = 1:0.5*width;
    segparams = struct('nucr', s.nucr,'blurradius', s.blurradius,'debrisarea', s.debrisarea,...
        'boulderarea',s.boulderarea, 'blobthreshold', s.blobthreshold);
    trackparams = {s.nucr,s.maxjump,s.debrisarea,s.masschangethreshold,s.areachangethreshold,s.daughtervariance};
    
    %%% Load bgcmos/s.bias
    if ~isempty(s.bgcmospath)
        load(s.bgcmospath,'cmosoffset');
        bgcmos  =  cmosoffset;
    else
        bgcmos = 100 * ones(height,width);
    end
    
    for i = 1:length(s.signals)
        if s.bias(i)
            load(fullfile(s.biaspath,[s.signals{i},'_',num2str(site),'.mat']));
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
            switch s.microscope
                case 'IXM'
                    rawdir = fullfile(s.imagepath,shot);
                    fileName = [shot '_' s.signals{j},'_',num2str(f),'.tif'];
                    raw{j} = single(imread(fullfile(rawdir,fileName)));
                case 'nikon'
                    rawdir = fullfile(s.imagepath,shot);
                    shotSearch = sprintf(s.formatCode,f-1,rowColumnTowellName(row,col),site-1);
                    fileName = findFile(rawdir,shotSearch);
                    raw{j} = single(imread(fullfile(rawdir,fileName),j));
            end
        
            if s.bgcmoscorrection
                raw{j} = (raw{j}-bgcmos);
                raw{j}(raw{j} < 0.001) = 0.001;
            end
            if s.bias(j)
                raw{j} = raw{j}./bias_cell{j};
            end
            if s.postbin> 0
                raw{j} = imresize(raw{j},s.postbin);
                [height,width]=size(raw{j});
            end
            if s.signal_foreground(j)
                normnuc = (raw{s.nucChannel}-prctile(raw{s.nucChannel}(:),1));
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
            nuc_mask = segmentImage(raw{s.nucChannel},s.firstsegmethod, segparams); 
        else
            if isequal(s.secondsegmethod, 'apriori')
                nuc_mask = threshmask(raw{s.nucChannel});
            else
                nuc_mask = segmentImage(raw{s.nucChannel},s.secondsegmethod, segparams);
            end
            %%% Calculate jitter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            lastgoodframe = find(badframes == 0,1,'last');
            if s.register | any(f == s.register_exception)  %register images if set
                [reljitx,reljity] = registerimages(imfill(extractmask(regheight,regwidth),'holes'),nuc_mask(regheight,regwidth));
            else
                reljitx = 0;
                reljity = 0;
            end
            jitters(f,:) = jitters(lastgoodframe,:)+[reljitx,reljity];
            if isequal(s.secondsegmethod, 'apriori')
                [nuc_mask,marker_mask] = apriori_markermask(nuc_mask,nuc_center,jitters(f,:));
            end

        end
        %nuc_mask = imopen(nuc_mask,strel('disk',3));
        whole_mask = nuc_mask;
        nuc_mask = bwareaopen(nuc_mask,s.debrisarea);
        
        %%% adaptive bridging settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [bounds,label_mask]=bwboundaries(nuc_mask,'noholes');
        mask_props = regionprops(label_mask,raw{s.nucChannel},'Area');
        sizes = [mask_props(:).Area];
%         solid = [mask_props(:).Solidity];
        medSize = prctile(sizes,50);
        gate = sizes >= medSize*s.split_mult;% | solid < s.soliditythresh;
        

        %%% Mitotic gate
%         mask_props = regionprops(nuc_mask,raw{s.nucChannel},'Area','Circularity','PixelIdxList');
%         inten = [];
%         varInt = [];
%         varIm = stdfilt(raw{s.nucChannel},true(5));
%         nuc_mask_erode = imerode(nuc_mask,strel('disk',5));
%         varIm(~nuc_mask_erode) = NaN;
%         for n = 1:length(mask_props)
%             pix = mask_props(n).PixelIdxList;
%             inten(n) = mean(raw{s.nucChannel}(pix));
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
        nuc_mask = segmentdeflections_gated(nuc_mask,s.nucr, bounds,gate);
        nuc_mask = bwareaopen(nuc_mask,s.debrisarea);
        %nuc_mask = excludelargeandwarped_3(nuc_mask,s.boulderarea,s.soliditythresh);
        
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
        tempframe(:,:,1) = imadjust(mat2gray((raw{s.nucChannel})));
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
        clear raw mask_total foreground mask_total blurIm label_mask bounds
        
        %% Extract features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [nuc_label,numcells] = bwlabel(nuc_mask);
        nuc_info = regionprops(nuc_mask,real{s.nucChannel},'Area','Centroid','MeanIntensity');
        nuc_area = vertcat(nuc_info.Area);
        
        %%% detect bad frame
        mednuc = median(nuc_area);
        if i>firstgoodindex+1 && (numcells<numthresh || mednuc>blurthreshhigh || mednuc<blurthreshlow)
            fprintf('badframe: frame %0.0f\n',f);
            badframes(f) = 1;
            extractmask = bwmorph(nuc_mask,'remove');
            if s.maskwrite
                imwrite(uint16(extractmask),fullfile(maskdir,[shot,'_',s.maskname,num2str(f),'.tif']));
            end
            continue;
        end
        blurthreshhigh = 1.2*mednuc;
        blurthreshlow = 0.8*mednuc;
        numthresh = 0.5*numcells;
        nuc_center = vertcat(nuc_info.Centroid);
        nuc_density = vertcat(nuc_info.MeanIntensity);
        
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
                adaptivetrack_NR(f,lastgoodframe,f,tracedata,curdata,tracking,real{s.nucChannel},nuc_label,jitters(f,:),trackparams,debugpackage);
            badframes(f) = 0;
            
            %%% Check tracking
            if debug_mode
                keyboard;
            end
            %{
            matchinglabel=prev_label==nuc_label & nuc_label;
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
                imwrite(uint16(extractmask),fullfile(maskdir,[shot,'_',s.maskname,num2str(f),'.tif']));
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
                p = punctaModule(unblurred{j}, nuc_label, nanvec, s.punctaThresh, s.punctaFiltSize, s.varThresh,debug_mode, numlivecells, nuc_info,cellid);
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
            cell2mat(p.filtMaskedMean), cell2mat(p.filtMaskedArea),cell2mat(p.filtMaskedNum), ...
            cell2mat(p.varMaskedMean), cell2mat(p.varMaskedArea), cell2mat(p.varMaskedNum)];
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
    save(fullfile(s.savepath,['tracedata_',shot,'.mat']),'tracedata','genealogy','jitters','-v7.3');
    names = output_names_puncta_number(s.signals,s.localbgmeasure, s.ringcalc, s.punctacalc, s.punctaThresh, s.varThresh);
    save(fullfile(s.savepath,'settings_live.mat'),'s','names'); %saves all the settings to the Data folder. will overwrite with the most recent
    
    toc(timetotal);
end
end


function p = punctaModule(punctaSig, nuc_label, nanvec, punctaThresh, filtSize, varThresh, debug_mode, numlivecells, nuc_info,cellid)
p = createPunctaStruct();

%% Initialize PCNA storage
p.PCNAMean = nanvec;
p.filtMean = nanvec;
p.varMean = nanvec;
p.varStd = nanvec;
for t = 1:length(punctaThresh)
    p.filtMaskedMean{t} = nanvec;
    p.filtMaskedArea{t} = nanvec;
    p.filtMaskedNum{t} = nanvec;
end
for t = 1:length(varThresh)
    p.varMaskedMean{t} = nanvec;
    p.varMaskedArea{t} = nanvec;
    p.varMaskedNum{t} = nanvec;

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
punctaFilt = imtophat(punctaSig, strel('disk',filtSize,0));
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
    punctaFiltMasked = punctaFilt .* filtMask;
    filt_label = bwlabel(filtMask);

    for n = 1:numlivecells
        cc = cellid(n);
        p.filtMaskedMean{t}(n) = nanmean(punctaFiltMasked(nuc_info(cc).PixelIdxList));
        p.filtMaskedArea{t}(n) = sum(filtMask(nuc_info(cc).PixelIdxList));
        numPuncta = length(unique(filt_label(nuc_info(cc).PixelIdxList))) - 1;
        p.filtMaskedNum{t}(n)  = numPuncta;
    end
end
for t = 1:length(varThresh)
    varMask = punctaVar > varThresh(t) & analysis_mask;
    punctaVarMasked = punctaVar .* varMask;
    var_label = bwlabel(varMask);

    for n = 1:numlivecells
        cc = cellid(n);
        p.varMaskedMean{t}(n) = nanmean(punctaVarMasked(nuc_info(cc).PixelIdxList));
        p.varMaskedArea{t}(n) = sum(varMask(nuc_info(cc).PixelIdxList));
        numPuncta = length(unique(var_label(nuc_info(cc).PixelIdxList))) - 1;
        p.varMaskedNum{t}(n)  = numPuncta;
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
p.filtMaskedNum{1} = [];
p.varMaskedMean{1} = [];
p.varMaskedArea{1} = [];
p.varMaskedNum{1} = [];

end


function mask = segmentImage(image, segMethod, s)
switch segMethod
    case 'concavity'
        %imblur = imfilter(image,fspecial('gaussian',s.blurradius),'symmetric');
        imblur = imfilter(image,fspecial('gaussian',s.blurradius),'symmetric');
        mask = logical(ThreshImage(imblur));
        mask = imopen(mask,strel('disk',3)); % make mask less jagged
        mask = (imfill(mask,'holes'));
        mask_noholes = (imfill(mask,'holes'));
        holes = mask_noholes & ~mask;
        holes = bwareaopen(holes,round(s.nucr^2*pi));
        mask = mask_noholes & ~holes;

    case 'log'
        mask = blobdetector_4(log(image),s.nucr,s.blobthreshold,s.debrisarea);
    case 'log contour'
        mask = blobdetector_4(log(image),s.nucr,s.blobthreshold,s.debrisarea);
        imblur = imfilter(image,fspecial('gaussian',5),'symmetric');
        mask = activecontour(log(imblur),mask,10);
        mask = imopen(mask,strel('disk',3)); % make mask less jagged
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