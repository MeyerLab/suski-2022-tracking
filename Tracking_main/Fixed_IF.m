function Fixed_IF(settings, row,col,site,debug_mode,varargin)
%Fixed_IF  IF segmetation analysis. Inputs parameter struct s, row,
%          col, site. debug_mode for setting up analysis, varargin not used
%            Compatible with IX micro and nikon naming schemes. 
%% REQUIRED SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% Paths
% settings.data_path = 'F:\Data\C-Cdt1\C198-live\Data\';      % Output location of mat files
% settings.IF_imagesessions = {'H:\8TB4\Data\C198-IF1\','H:\8TB4\Data\C198-IF2\'};    
%                                                             % Cell array containing all imaging sessions to be matched
% settings.bgcmospath = '';
% settings.crop_save = 'H:\8TB4\Data\C198-IFcrop\';           % Output cropped images, leave empty if no save
% settings.mask_save = 'H:\8TB4\Data\C198-IFcrop\';           % Output mask, leave empty if no save
% settings.maskname = 'mask';
% settings.primaryMaskRound = 1;
% settings.maskIndex = [1 1];

% %%% General parameters
% settings.microscope = 'nikon';
% % Formatting code for nikon files, 1$ is timepoint, 2$ is well, 3$ is point
% settings.formatCodeLive = 'Time%1$05d_Well%2$s_Point%2$s_%3$04d_*.tiff';
% settings.formatCodeFixed = 'Time%1$05d_Well%2$s_Point%2$s_%3$04d_*.tiff';
% settings.magnification = 20;     % 10 = 10x or 20 = 20x
% settings.binsize = 1;            % 1 = bin 1 or 2 = bin 2
% settings.postbin = [.5 .5];           % 0 for no software bin, number for scaling factor (e.g. 0.5 for 2x2 bin);
% settings.signals = {{'DAPI','YFP','FarRed'}, {'DAPI','YFP','FarRed'}};  % Each imaging session in a cell inside cell array
%
% %%% Segmentation parameters
% settings.segmethod = {'concavity','thresh'};  % 'log' or 'single' or 'double'
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
% settings.bias = {[1 1 1], [1 1 1]};
% settings.biasall = {[0 0 0], [0 0 0]};
% settings.sigblur = {[0 0 0], [0 0 0]};
% settings.signal_foreground = {[0 0 0], [0 0 0]};
% settings.bleedthrough = {[0 0 0], [0 0 0]};   % Calculate bleedthrough for channel
% settings.bleedthroughslope = {};              % Cell array of cell arrays
% settings.bleedthroughoff = {};                % Cell array of cell arrays
% settings.ringcalc = {[1 1 1], [0 1 1]};
% settings.ringthresh = {[0 0 0 ], [0 0 0]};    % Threshold for foreground ring pixels
% settings.punctacalc = {[0 0 0], [0 0 0]};     % Analyze puncta in nucleus
% settings.punctaThresh = {{[],[],[]},{}};      % Threshold for tophat filter for puncta
% settings.punctatopsize = 2;                   % Top hat filter size
% settings.cytopuncta = {[0 1 0]};
% settings.thickenradius = 2*settings.nucr; 

% settings.localbg = {[0 0 0], [0 0 0]};        % Measure local background in box around cell
% settings.minringsize = 100;
% settings.bgcmoscorrection = 1;                % 1 = correct for shading; 0 = dont correct for shading;
% settings.compression = 4;
% settings.bgsubmethod = {{'global nuclear','global nuclear', 'global nuclear'}, ...
%     {'global nuclear', 'global nuclear', 'global nuclear'}}; % 'global nuclear','global cyto','tophat','semi-local nuclear'
% settings.bgperctile = {[25 25 25], [25 25 25]};  % Set bg as perctile for each channel
% settings.frameIF=1;
%  


%% PROCESS IMAGES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timetotal = tic;
%%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set paths
datadir = settings.data_path;
shot = [num2str(row),'_',num2str(col),'_',num2str(site)];
% break if excluded or already exisits
if exist(fullfile(datadir,['IFdata_',shot,'.mat']))
    return
end
if ~exist(datadir,'dir')
    mkdir(datadir);
end
maskdir = fullfile(settings.mask_save,shot);
if ~exist(maskdir,'dir') && ~isempty(settings.mask_save)
    mkdir(maskdir);
end

%%% Set channnels
names = settings.signals;

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
            load(fullfile(biasdir,[names{i}{j},'_',num2str(site),'.mat']));
            bias_cell{i}{j} = bias;
        end
    end
end

%% Load images and segment %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(names)
    for j = 1:length(names{i})
        switch settings.microscope
            case 'IXM'
                rawdir = fullfile(settings.IF_imagesessions{i},'Raw',shot);
                fileName = [shot '_' names{i}{j},'_',num2str(1),'.tif'];
                raw{i}{j} = single(imread(fullfile(rawdir,fileName)));
            case 'nikon'
                rawdir = fullfile(settings.IF_imagesessions{i},'Raw',shot);
                shotSearch = sprintf(settings.formatCodeFixed,0,rowColumnTowellName(row,col),site-1);
                fileName = findFile(rawdir,shotSearch);
                raw{i}{j} = double(imread(fullfile(rawdir,fileName),j));
        end
        [height,width]=size(raw{i}{j});
        if settings.bgcmoscorrection
            raw{i}{j} = (raw{i}{j}-cmosoffset);
            raw{i}{j}(raw{i}{j}<0) = 0;
        end
        if settings.bias{i}(j)
            raw{i}{j} = raw{i}{j}./bias_cell{i}{j};
        end
        if settings.postbin(i)
            raw{i}{j} = imresize(raw{i}{j},settings.postbin(i));
            [height,width]=size(raw{i}{j});
        end
    end
    
    %% Segment nuclei %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    segmethod = settings.segmethod{i};
    blurradius = settings.blurradius;
    switch segmethod
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
            blur=imfilter(raw{i}{1},fspecial('gaussian',blurradius),'symmetric');
            nuc_mask{i}=ThreshImage(blur);
            nuc_mask{i}=imfill(nuc_mask{i},'holes');
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
%% debugging: view images %%%%%%%%%%
%{
        %check mask with first channel
          for session = 1:length(names)
                extractmask = bwmorph(nuc_mask{session},'remove');
                tempframe = imadjust(mat2gray(raw{session}{settings.maskIndex(session)}));
                tempframe(:,:,2) = extractmask;
                tempframe(:,:,3) = 0;
                figure,imshow(tempframe);
          end
        nuc_info = struct2cell(regionprops(nuc_mask{1},'Area')');
        nuc_area = squeeze(cell2mat(nuc_info(1,1,:)));
        figure, histogram(nuc_area,100);
    %}
    
%% Align images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate jitter between IF rounds
jitmatx = [];
jitmaty = [];

if length(names) > 1
    regheight = 1:0.5*height;
    regwidth = 1:0.5*width;
    for i = 1:length(names)
        [reljitx,reljity] = registerimages(nuc_mask{1}(regheight,regwidth),nuc_mask{i}(regheight,regwidth));
        jitmatx = [jitmatx; reljitx];
        jitmaty = [jitmaty; reljity];
    end
    
    %%% Error if large offset
    if max(jitmatx) > settings.maxjit | max(jitmaty) > settings.maxjit
        disp(['Max jitter too high ', num2str(max(jitmatx)),', ', num2str(max(jitmaty))]);
        regheight = 1:height;
        regwidth = 1:width;
        jitmatx = [];
        jitmaty = [];
        for i = 1:length(names)
            [reljitx,reljity] = registerimages(nuc_mask{1}(regheight,regwidth),nuc_mask{i}(regheight,regwidth));
            jitmatx = [jitmatx; reljitx];
            jitmaty = [jitmaty; reljity];
        end        
        disp(['Fixed to ', num2str(max(jitmatx)),', ', num2str(max(jitmaty))]);
        
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

%%% Clear unnecessary variables
clear nuc_mask normnuc normraw sumraw exclude_mask

%% Subtract background %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
compression = settings.compression;
[height,width]=size(nuc_mask_aligned);
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
    clear nanmask_aligned nanmaskcyto_aligned whole_mask_aligned foreground_aligned
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
            if settings.punctacalc{i}(j)                
                unblurred{i}{j} = unblurred{i}{j}-realbleedthrough;
            end
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
extractmask = boundarymask(nuc_mask_aligned);
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
    tempframe = imadjust(mat2gray(raw{session}{settings.maskIndex(session)}));
    tempframe(:,:,2) = extractmask;
    tempframe(:,:,3) = 0;
    figure,imshow(tempframe);
end

%%% Check foreground
round = 2;
channel = 2;
extractmask = boundarymask(foreground_aligned{round}{channel});
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
    for i = 1:length(names)
        for j = 1:length(names{i})
            save_dir = fullfile(settings.crop_save,shot);
            if ~exist(save_dir)
                mkdir(save_dir)
            end
            save_name = fullfile(save_dir,[ shot,'_' names{i}{j} '_round' num2str(i) '.tif']);
            imwrite(uint16(real{i}{j}/4), save_name);
        end
    end
    %%% Save mask
    extractmask = bwmorph(nuc_mask_aligned,'remove');
    save_name = fullfile(save_dir,[ shot,'_' 'mask' '.tif']);
    imwrite(uint16(extractmask),save_name);
end

%%% Save mask
if ~isempty(settings.mask_save)
    extractmask = bwmorph(nuc_mask_aligned,'remove');
    imwrite(uint16(extractmask),fullfile(maskdir, [ shot '_' settings.maskname '.tif']));
end

%% Label masks %%%%%^%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nuc_label = bwlabel(nuc_mask_aligned);
nuc_info = regionprops(nuc_label,'PixelIdxList','Centroid','Area');
nuc_area = vertcat(nuc_info.Area);
nuc_center = vertcat(nuc_info.Centroid);
numcells = numel(nuc_area);

filtered_label = bwlabel(exclude_mask_aligned);
filtered_label(filtered_label > 0) = filtered_label(filtered_label > 0) + numcells;
all_obj = nuc_label + filtered_label;


%% Create secondary masks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Cyto ring
if any(cell2mat(settings.ringcalc))
    if (settings.magnification == 10 & settings.binsize == 1) | ...
            (settings.magnification == 20 & (settings.binsize == 2 | settings.postbin))
        innerrad = 1; outerrad = 5;
    else
        innerrad = 2; outerrad = 10;
    end
    ring_label = getcytoring_thicken(all_obj,innerrad,outerrad,[]);
    ring_info = regionprops(ring_label,'PixelIdxList');
    
    %%% Check cytoring
    if debug_mode
        keyboard;
    end
    %{
            %%% debugging: view images %%%%%%%%%%
            temp_ring  =  ring_label > 0;
            extractmask = bwmorph(temp_ring,'remove');
            signal  =  real{1}{1};
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
            signal  =  real{1}{1};
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
        real_temp = real{i}{j};
        if settings.localbg{i}(j)
            real_temp_masked = real_temp;
            real_temp_masked(nanmask_aligned{i}{j}) = NaN;
        end       
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
                punctaMask = punctaFilt{i}{j} > settings.punctaThresh{i}{j}(t);
                punctaIntensity_Masked = punctaFilt{i}{j} .* punctaMask;
                puncta_label = bwlabel(punctaMask);
                
                %%% Insert code to output puncta masks here %%%
                for cc = 1:numcells
                    punctaArea{i}{j}(cc,t) = sum(punctaMask(mask_info(cc).PixelIdxList));
                    punctaIntensity{i}{j}(cc,t) = sum(punctaIntensity_Masked(mask_info(cc).PixelIdxList));
                    numPuncta = length(unique(puncta_label(mask_info(cc).PixelIdxList))) - 1;
                    punctaNumber{i}{j}(cc,t) = numPuncta;
                    
                end
            end
        end
    end
end

%% store data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IFdata = [nuc_center(:,1),nuc_center(:,2),nuc_area,...
    cell2mat([sigmean{:}]),cell2mat([sigmedian{:}]),cell2mat([localbg{:}]),...
    cell2mat([sigringmedian{:}]),...
    cell2mat([punctaArea{:}]), cell2mat([punctaIntensity{:}]), cell2mat([punctaNumber{:}])];
save(fullfile(datadir,['IFdata_',shot,'.mat']),'IFdata');
header = output_names_multiIF_puncta_fixed(settings.signals,settings.localbg, settings.ringcalc,settings.punctacalc,settings.punctaThresh);
save(fullfile(datadir,'settings_IF.mat'),'settings','header'); %saves all the settings to the Data folder. will overwrite with the most recent
toc(timetotal);

close all;