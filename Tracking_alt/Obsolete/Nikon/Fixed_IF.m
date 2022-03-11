function Fixed_IF(row,col,site,debug_mode)
% row = 2;
% col = 2;
% site = 2;
% debug_mode = 0;
settings.debug_mode = debug_mode;

%% SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Paths
settings.data_path = 'F:\Data\C-Cdt1\C163-live\Data\';
settings.IF_imagesessions = {'H:\4TB6\Data\C163-IF1\','H:\4TB6\Data\C163-IF2\'};
settings.bgcmospath = 'F:\MATLAB\BGimages\IX_XL\cmosoffset_bin1.mat';
settings.maskpath = 'H:\4TB6\Data\C163-IFcrop\';
settings.crop_save = 'H:\4TB6\Data\C163-IFcrop\';

%%% File Options
settings.maskwrite_option = 1; %1 = save an image with the nucmask; 0 = dont save an image with the nucmask
maskname = 'nucedge_';

%%% General parameters
settings.magnification = 20; %10 = 10x or 20 = 20x
settings.binsize = 1; %1 = bin 1 or 2 = bin 2
settings.postbin = .5;
settings.signals = {{'DAPI_','YFP_','FarRed_'}, {'DAPI_','YFP_'}};

%%% Quantification parameters
settings.bias = {[1 1 1],[1 1]};
settings.biasall = {[0 0 0], [0 0]};
settings.sigblur = {[0 0 0], [0 0]};
settings.signal_foreground = {[0 1 1],[0 1]};
settings.ringcalc = {[1 1 1],[0 1]};
settings.ringthresh = {[100 100 100 ], [100 100 ]};
settings.punctacalc = {[0 0 0 ],[0 0]};
settings.punctaThresh = {{[],[],[]},{[],[]}};
settings.punctatopsize = 2;
settings.localbg = {[0 0 0],[0 0]};
settings.minringsize = 100;
settings.bleedthrough = {[0 0 0],[0 0]};
settings.bleedthroughpth = {};
settings.bgcmoscorrection = 1; %1 = correct for shading; 0 = dont correct for shading;
settings.bgsubmethod = {{'global nuclear','global cyto', 'global nuclear'}, ...
     {'global nuclear', 'global cyto'}}; %'global nuclear','global cyto','tophat','semi-local nuclear'
settings.bgperctile = {[25 25 25],[25 25]};

%%% Segmentation parameters
settings.segmethod = {'concavity', 'thresh'}; %'log' or 'single' or 'double'
settings.nucr = 12; %10x bin1 = 12 and 20x bin2 = 12
settings.debrisarea = 100; %100
settings.boulderarea = 1500; %1500 % less restrictive for C152 for >4N cells
settings.blobthreshold = -0.03;
settings.blurradius = 3;
settings.soliditythresh = 0.80;
settings.compression = 4;
settings.exclude = {};

%% PROCESS IMAGES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timetotal = tic;
%%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set paths
datadir = settings.data_path;
shot = [num2str(row),'_',num2str(col),'_',num2str(site)];
% break if excluded or already exisits
if exist([datadir,'IFdata_',shot,'.mat']) | any(ismember(settings.exclude,shot))
    return
end
if ~exist(datadir,'dir')
    mkdir(datadir);
end
maskdir = [settings.maskpath,shot];
if ~exist(maskdir,'dir') && settings.maskwrite_option
    mkdir(maskdir);
end

%%% Set channnels
names = settings.signals;
nucname = names{1}{1};

%%% Set segmentation parameters
nucr = settings.nucr;
debrisarea = settings.debrisarea;
boulderarea = settings.boulderarea;
blobthreshold = settings.blobthreshold;
solidity = settings.soliditythresh;

%% Load auxilliary files %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Load bgcmos
load([settings.bgcmospath],'cmosoffset');
bgcmos = cmosoffset;
[height,width]=size(bgcmos);

%%% Load bias
for i = 1:length(names)
    for j = 1:length(names{i})
        if settings.biasall{i}(j)
            biasdir = [settings.IF_imagesessions{i} 'Bias_all\'];
        else
            biasdir = [settings.IF_imagesessions{i} 'Bias\'];
        end
        if settings.bias{i}(j)
            load([biasdir,names{i}{j},num2str(site),'.mat']);
            bias_cell{i}{j} = bias;
        end
    end
end

%% Load images and segment %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(names)
    rawdir = [settings.IF_imagesessions{i},'Raw\',shot,'\',shot,'_'];
    for j = 1:length(names{i})
        raw{i}{j} = single(imread([rawdir,names{i}{j},num2str(1),'.tif']));
        if settings.bgcmoscorrection
            raw{i}{j} = (raw{i}{j}-bgcmos);
            raw{i}{j}(raw{i}{j}<0) = 0;
        end
        if settings.bias{i}(j)
            raw{i}{j} = raw{i}{j}./bias_cell{i}{j};
        end
        if settings.postbin
            raw{i}{j} = imresize(raw{i}{j},settings.postbin);
            [height,width]=size(raw{i}{j});
        end
    end
    
    %% Segment nuclei %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    segmethod = settings.segmethod{i};
    blurradius = settings.blurradius;
    switch segmethod
        case 'log'
            nuc_mask{i} = blobdetector_4(log(raw{i}{1}),nucr,blobthreshold,debrisarea);
        case 'thresh'
            nuc_mask{i} = threshmask(raw{i}{1},blurradius);
        case 'single'
            nuc_mask{i} = threshmask(raw{i}{1},blurradius);
            nuc_mask{i} = markershed_filter(nuc_mask{i},round(nucr*2/3),4);
        case 'double'
            nuc_mask{i} = threshmask(raw{i}{1},blurradius);
            nuc_mask{i} = markershed_filter(nuc_mask{i},round(nucr*2/3),4);
            nuc_mask{i} = secondthresh(raw{i}{1},blurradius,nuc_mask{i},boulderarea*2);
        case 'double marker'
            nuc_mask{i} = threshmask(raw{i}{1},blurradius);
            nuc_mask{i} = secondthresh_all(raw{i}{1},blurradius,nuc_mask{i});
            nuc_mask{i} = markershed_filter(nuc_mask{i},round(nucr*.4),4);
        case 'multithresh'
            nuc_mask{i} = threshmask_multi(raw{i}{1},blurradius,3, 1e-02);
            nuc_mask{i} = markershed_filter(nuc_mask{i},round(nucr*.4),4);
            %nuc_mask{i} = secondthresh(raw{i}{1},blurradius,nuc_mask{i},boulderarea);
        case 'concavity'
            blur=imfilter(raw{i}{1},fspecial('gaussian',blurradius),'symmetric');
            nuc_mask{i}=ThreshImage(blur);
            nuc_mask{i}=imfill(nuc_mask{i},'holes');
    end
    
    %%% Split and filter out nuclei only for primary mask
    if i==1
        whole_mask = logical(nuc_mask{i});
        nuc_mask{i} = bwareaopen(nuc_mask{i},debrisarea);
        nuc_mask{i} = segmentdeflections_bwboundaries(nuc_mask{i},nucr,debrisarea);
        nuc_mask{i} = excludelargeandwarped_3(nuc_mask{i},boulderarea,solidity);
    end
    nuc_mask{i} = bwareaopen(nuc_mask{i},debrisarea);
    nuc_mask{i} = imclearborder(nuc_mask{i});
end

%%% Clear unnecessary variables
clear bias bias_cell cmosoffset imblur

%% Check  segmentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if settings.debug_mode
    keyboard;
end
%% debugging: view images %%%%%%%%%%
%{
        %check mask with first channel
        extractmask = bwmorph(nuc_mask{1},'remove');
        tempframe = imadjust(mat2gray(raw{1}{1}));
        tempframe(:,:,2) = extractmask;
        tempframe(:,:,3) = 0;
        figure,imshow(tempframe);
       
        nuc_info = struct2cell(regionprops(nuc_mask{1},'Area')');
        nuc_area = squeeze(cell2mat(nuc_info(1,1,:)));
        figure, histogram(nuc_area,100);
    %}
    
%% Align images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate jitter between IF rounds
regheight = 1:0.5*height;
regwidth = 1:0.5*width;
jitmatx = [];
jitmaty = [];
if length(names) > 1
    for i = 1:length(names)
        [reljitx,reljity] = registerimages(nuc_mask{1}(regheight,regwidth),nuc_mask{i}(regheight,regwidth));
        jitmatx = [jitmatx; reljitx];
        jitmaty = [jitmaty; reljity];
    end
else
    jitmatx = [jitmatx 0];
    jitmaty = [jitmaty 0];
end

%%% Align and crop IF rounds
cropcoors = getcropcoors([height width], jitmatx, jitmaty);
for i = 1:length(names)
    for j = 1:length(names{i})
        raw{i}{j} = raw{i}{j}(cropcoors(i+1,1):cropcoors(i+1,2),cropcoors(i+1,3):cropcoors(i+1,4));
        if settings.signal_foreground{i}(j) %calculate foreground based on sum of signal
            normnuc = (raw{1}{1}-prctile(raw{1}{1}(:),1));
            normnuc(normnuc < 1) = 1;
            normnuc = normnuc/prctile(normnuc(:),99);
            normraw = (raw{i}{j}-prctile(raw{i}{j}(:),1));
            normraw(normraw < 1) = 1;
            normraw = normraw/prctile(normraw(:),99);
            sumraw = normnuc + normraw;
            foreground_aligned{i}{j} = threshmask(sumraw,blurradius);
            if sum(foreground_aligned{i}{j}(:))/length(foreground_aligned{i}{j}(:)) > .9 % account for large foregrounds
                foreground_aligned{i}{j} = zeros(height, width);
            end
        end
    end
end
nuc_mask_aligned = nuc_mask{1}(cropcoors(1,1):cropcoors(1,2),cropcoors(1,3):cropcoors(1,4));
whole_mask_aligned = whole_mask(cropcoors(1,1):cropcoors(1,2),cropcoors(1,3):cropcoors(1,4));

%%% Clear unnecessary variables
clear nuc_mask normnuc normraw sumraw

%% Subtract background %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
compression = settings.compression;
for i=1:length(names)
    for j = 1:length(names{i})
        %%% Set foreground mask
        if settings.signal_foreground{i}(j)
            mask_total  =  whole_mask_aligned | foreground_aligned{i}{j};
            nanmask_aligned{i}{j}  = imdilate(mask_total,strel('disk',nucr));
            nanmaskcyto_aligned{i}{j} = imdilate(mask_total,strel('disk',nucr*2));
        else
            nanmask_aligned{i}{j}  = imdilate(whole_mask_aligned,strel('disk',nucr));
            nanmaskcyto_aligned{i}{j} = imdilate(whole_mask_aligned,strel('disk',nucr*2));
        end
        
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
    %%% correct IF for bleedthrough
    for j = 1:length(names{i})
        if settings.bleedthrough{i}(j)
            load([settings.bleedthroughpth{i}{j}],'bleedthroughrate');
            realbleedthrough = real{i}{j+1}*bleedthroughrate(2);
            real{i}{j} = real{i}{j}-realbleedthrough;
        end
    end
end

%% Check alignment segmentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if settings.debug_mode & length(names) > 1
    keyboard;
end
%{
%%% Check nuclear mask with primary nuclear channel
extractmask = bwmorph(nuc_mask_aligned,'remove');
tempframe = imadjust(mat2gray(raw{1}{1}));
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
    tempframe = imadjust(mat2gray(raw{session}{1}));
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
tempframe=imadjust(mat2gray(real{1}{1}));
tempframe(:,:,2)=imadjust(mat2gray(real{2}{1}));
tempframe(:,:,3)=0;
figure,imshow(tempframe)
%}

%%% Clear unnecessary variables
clear nuc_mask blur mask_total raw whole_mask foreground_aligned 

%% Puncta processing
for i = 1:length(names)
    for j = 1:length(names{i})
        if settings.punctacalc{i}(j)
            punctaFilt = imtophat(unblurred{i}{j}, strel('disk',settings.punctatopsize,0));
            erodeMask = imerode(nuc_mask_aligned,strel('disk',3,0));
            if settings.debug_mode
                keyboard;
            end
            for t = 1:length(settings.punctaThresh{i}{j})
                punctaMask{i}{j}{t} = punctaFilt > settings.punctaThresh{i}{j}(t) & erodeMask;
                punctaIntensity{i}{j}{t} = punctaFilt .* punctaMask{i}{j}{t};
                puncta_label{i}{j}{t} = bwlabel(punctaMask{i}{j}{t});
            end
        end
    end
end

%% Save cropped/corrected images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(settings.crop_save) & ~settings.debug_mode
    for i = 1:length(names)
        for j = 1:length(names{i})
            save_dir = [settings.crop_save,shot,'\'];
            if ~exist(save_dir)
                mkdir(save_dir)
            end
            save_name = [save_dir shot,'_' names{i}{j} '_round' num2str(i) '.tif'];
            imwrite(uint16(real{i}{j}/4), save_name);
        end
    end
    %%% Save mask
    extractmask = bwmorph(nuc_mask_aligned,'remove');
    save_name = [save_dir shot,'_' 'mask' '.tif'];
    imwrite(uint16(extractmask),save_name);
end

%%% Save mask
if settings.maskwrite_option
    extractmask = bwmorph(nuc_mask_aligned,'remove');
    imwrite(uint16(extractmask),[maskdir '\' maskname shot '.tif']);
end
%% MEASSURE FEATURES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nuc_info = struct2cell(regionprops(nuc_mask_aligned,'Area','Centroid')');
nuc_area = squeeze(cell2mat(nuc_info(1,1,:)));
nuc_center = squeeze(cell2mat(nuc_info(2,1,:)))';
nuc_label = bwlabel(nuc_mask_aligned);
numcells = numel(nuc_area);
nuc_info = regionprops(nuc_label,'PixelIdxList','Centroid');

%% Initialize storage variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Ring variables
if any(cell2mat(settings.ringcalc))
    if (settings.magnification == 10 & settings.binsize == 1) | ...
            (settings.magnification == 20 & (settings.binsize == 2 | settings.postbin))
        innerrad = 1; outerrad = 5;
    else
        innerrad = 2; outerrad = 10;
    end
    ring_label = getcytoring_thicken(nuc_label,innerrad,outerrad,real{1}{1});
    ring_info = regionprops(ring_label,'PixelIdxList');
    %%% Check cytoring
    if debug_mode
        keyboard;
    end
    %{
        %%% debugging: view images %%%%%%%%%%
        temp_ring  =  ring_label > 0;
        extractmask = bwmorph(temp_ring,'remove');
        signal  =  real{1}{2};
        intensities  =  signal.*temp_ring;
        [height,width] = size(real{1}{1});
        tempframe = zeros(height,width,3);
        tempframe(:,:,1) = imadjust(mat2gray(real{1}{2}));
        tempframe(:,:,2) = extractmask;;
        %tempframe(:,:,3) = marker_mask;
        figure,imshow(tempframe);
    %}
end

%%% Storage variables
nanvec = ones(numcells,1)*NaN;
sigringmedian{1} = [];
localbg{1} = [];
punctaArea{1}{1} = [];
punctaIntensity_store{1}{1} = [];
punctaNumber{1}{1} = [];

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
            punctaIntensity_store{i}{j} = repmat(nanvec,1,length(settings.punctaThresh{i}{j}));
            punctaNumber{i}{j} = repmat(nanvec,1,length(settings.punctaThresh{i}{j}));
        end
    end
end

%% Measure cell features%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(names)
    for j = 1:length(names{i})
        %Measure features%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        real_temp = real{i}{j};
        real_temp_masked = real_temp;
        real_temp_masked(nanmask_aligned{i}{j}) = NaN;
        for cc = 1:numcells
            sigmean{i}{j}(cc) = mean(real_temp(nuc_info(cc).PixelIdxList));
            sigmedian{i}{j}(cc) = median(real_temp(nuc_info(cc).PixelIdxList));
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
                    if numel(ringall)>100
                        sigringmedian{i}{j}(cc) = nanmedian(ringforeground);
                    end
                end
                if settings.punctacalc{i}(j)
                    for t = 1:length(settings.punctaThresh{i}{j})
                        punctaArea{i}{j}(cc,t) = sum(punctaMask{i}{j}{t}(nuc_info(cc).PixelIdxList));
                        punctaIntensity_store{i}{j}(cc,t) = sum(punctaIntensity{i}{j}{t}(nuc_info(cc).PixelIdxList));
                        numPuncta = length(unique(puncta_label{i}{j}{t}(nuc_info(cc).PixelIdxList))) - 1;
                        punctaNumber{i}{j}(cc,t) = numPuncta;
                    end
                end 
            end
        end
    end
end


%% store data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IFdata = [nuc_center(:,1),nuc_center(:,2),nuc_area,...
    cell2mat([sigmean{:}]),cell2mat([sigmedian{:}]),cell2mat([localbg{:}]),...
    cell2mat([sigringmedian{:}]),...
    cell2mat([punctaArea{:}]), cell2mat([punctaIntensity_store{:}]), cell2mat([punctaNumber{:}])];
save([datadir,'IFdata_',shot,'.mat'],'IFdata');
header = output_names_multiIF_puncta_fixed(settings.signals,settings.localbg, settings.ringcalc,settings.punctacalc,settings.punctaThresh);
save([datadir,'settings_IF.mat'],'settings','header'); %saves all the settings to the Data folder. will overwrite with the most recent
toc(timetotal);

close all;