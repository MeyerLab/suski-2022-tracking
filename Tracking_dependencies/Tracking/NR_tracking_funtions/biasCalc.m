function biasCalc(p,debug)
%biasCalc Calculate bias based for either nikon or IX micro renamed files.
%   Takes in parameter struct s, and debug (shows sample segmented image). At end of
%   analysis, outputs visualiation of bias calculation per frame.
%%%% Required input parameters
%
%     Paramter struct s:
%     %%% Directories
%     p.biaspath='H:\8TB4\Data\C198-live\Bias';  % Output path for bias
%     p.imagepath='H:\8TB4\Data\C198-live\Raw';  % Image directory. If nikon, single folder for all sites
%     p.shadingpath='F:\MATLAB\BGimages\IX_Liu\cmosoffset_bin1.mat'; % cmos offset file for IX micro
%     p.names={
%         'CFP';
%         'YFP';
%         'RFP';
%         };
%         p.row_mat = [2:6];
%     p.col_mat = [2:11];
%     p.site_mat = [1:4];
%     p.frame_mat=[1];
%     p.biasAll = 0; % Save averaged bias across sites
%                    % Formatting code for nikon files, 1$ is timepoint, 2$ is well, 3$ is point
%     p.formatCode = 'Time%1$05d_Well%2$s_Point%2$s_%3$04d_*.tiff';
%     p.microscope = 'IXM';
%
%     %%% Settings
%     p.maskforeground = 1;   % 1 to mask nuclei
%     p.nucradius=12;%12;
%     p.blur_radius = 5;
%     p.adaptive = [0 0 0];   % use add channels together to segment foreground
%     p.dilate = {.75,1,1};   % multiples of nucradius, default .5, to expand mask
%     p.foreground_calc = 0;  % 1 for foreground calc (average foreground signal)
%     p.method = 'block';     % or pixel
%     p.blocknum= 15;
%     p.prctilethresh = 50;
%     p.compress = .25;
%     p.prctile_thresh=[0 100]; % pixel remove outliers
%     p.sigma = 25;             % pixel smooth window
%

%%% Initialize variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(p.shadingpath)
    load(p.shadingpath,'cmosoffset');
    bgcmos = cmosoffset;
    [height,width] = size(bgcmos);
else %no cmosoffset file, assume 100
    site = p.site_mat(1);
    col = p.col_mat(1);
    row = p.row_mat(1);
    shot = [num2str(row),'_',num2str(col),'_',num2str(site)];
    switch p.microscope
        %Nikon naming
        case 'nikon'
            shotSearch = sprintf(p.formatCode,0,rowColumnTowellName(row,col),site-1);
            fileName = findFile(fullfile(p.imagepath,shot),shotSearch);
            tempim = double(imread(fullfile(p.imagepath,shot,fileName),1));
            
        case 'IXM'
            fileName =  [shot '_' p.names{1},'_',num2str(1),'.tif'];
            tempim = double(imread(fullfile(p.imagepath,shot,fileName),1));
    end
    
    [height, width] = size(tempim);
    bgcmos = 100 * ones(height,width);
end

numsites = numel(p.site_mat);
if ~exist(p.biaspath,'dir')
    mkdir(p.biaspath);
end
if p.biasAll
    if ~exist([p.biaspath,'_all'],'dir')
        mkdir([p.biaspath,'_all']);
    end
end
bias_save = {};
block_save = {};
if strcmp(p.method,'block')
    for i = 1:length(p.names)
        block_save{i} = [];
    end
end


%% Calculate bias
for s = 1: numsites
    %%% Initialize intermediate variables
    switch p.method
        case 'pixel'
            accumulate_back = {};
            back_count = {};
            for i = 1:length(p.names)
                accumulate_back{i} = zeros(height,width);
                back_count{i} = zeros(height,width);
            end
        case 'block'
            biasstack = {};
            for i = 1:length(p.names)
                biasstack{i} = [];
            end
    end
    
    site = p.site_mat(s);
    for c = 1: numel(p.col_mat)
        col = p.col_mat(c);
        for r = 1:numel(p.row_mat)
            row = p.row_mat(r);
            for f = 1:numel(p.frame_mat)
                frame = p.frame_mat(f);
                shot = [num2str(row),'_',num2str(col),'_',num2str(site)];
                fprintf([shot,'\n']);
                
                %Find file name
                if strcmp(p.microscope, 'nikon')
                    shotSearch = sprintf(p.formatCode,frame-1,rowColumnTowellName(row,col),site-1);
                    fileName = findFile(fullfile(p.imagepath,shot),shotSearch);
                end
                
                %%% Read in all channels
                for i = 1:length(p.names)
                    try
                        switch p.microscope
                            %Nikon naming
                            case 'nikon'
                                raw{i} = double(imread(fullfile(p.imagepath,shot,fileName),i));
                            case 'IXM'
                                raw{i} = double(imread(fullfile(p.imagepath,shot,[shot '_' p.names{i},'_',num2str(frame),'.tif'])));
                        end
                    catch
                        warning([shot ' not found']);
                        continue;
                    end
                    raw{i} = raw{i} - bgcmos;
                    raw{i}(raw{i} < 1) = 1;                  
                end
                
                %%% Mask
                if p.maskforeground
                    % Calc and save nuclear mask if needed
                    if any(~p.adaptive)
                        mask = threshmask(raw{p.nuc_channel},p.blur_radius);
                        if p.foreground_calc
                            mask = ~mask;
                        end
                    end
                    % Save normalized nuclear signal if needed
                    if any(p.adaptive)
                        rawcomp = imresize(raw{p.nuc_channel}, p.compress);
                        norm_nuc = raw{p.nuc_channel}/prctile(rawcomp(:),90);
                    end
                    
                else
                    mask = false(size(raw{i}));
                end
                
                %%% Loop through channels for calc
                for i = 1:length(p.names)                    
                    %%% Mask Image
                    if p.maskforeground
                        if p.adaptive(i)
                            rawcomp = imresize(raw{i}, p.compress);
                            norm_raw = raw{i}/prctile(rawcomp(:),90);
                            overlay = norm_raw+norm_nuc;
                            overlay(overlay<1) = 1;
                            mask_ind = threshmask(log(overlay),p.blur_radius);
                            if p.foreground_calc
                                mask_ind = ~mask_ind;
                            end
                        else
                            mask_ind = mask;
                        end
                        mask_ind = imdilate(mask, strel('disk',ceil(p.nucradius*p.dilate{i})));
                    else
                        mask_ind = mask;
                    end
                    

                    blur = imfilter(raw{i},fspecial('disk',p.blur_radius),'symmetric');
                    blurnan = blur;
                    blurnan(mask_ind) = NaN;
                    
                    
                    if debug
                        extractmask=bwmorph(mask_ind,'remove');
                        [height,width]=size(mask_ind);
                        RGB=zeros(height,width,3);
                        RGB(:,:,1)=imadjust(mat2gray(raw{i}));
                        RGB(:,:,2)=extractmask;
                        figure,imshow(RGB);
                        keyboard;
                    end
                    
                    %%% Compile background
                    switch p.method
                        case 'pixel'
                            blurnancomp = imresize(blurnan, p.compress);
                            outliers = blurnan < prctile(blurnancomp(:),p.prctile_thresh(1)) | blurnan >prctile(blurnancomp(:),p.prctile_thresh(2));
                            mask_ind = mask_ind | outliers;
                            blurnan(mask_ind) = NaN;
                            blurnancomp = imresize(blurnan,p.compress);
                            blur(mask_ind)=0;
                            %                     imagesc(blur);
                            %                     title(p.names{i});
                            %                     pause;
                            blur = blur/nanmedian(blurnancomp(:));
                            accumulate_back{i} = accumulate_back{i} + blur;
                            back_count{i} = back_count{i} + ~mask_ind;
                        case 'block'
                            bgblock=blockpercentile_blockimage(blurnan,p.blocknum,p.prctilethresh);
                            midrc=ceil(p.blocknum/2);
                            refval=bgblock(midrc,midrc);
                            if ~isnan(refval)
                                bgblocknorm=bgblock/refval;
                                biasstack{i}=cat(3,biasstack{i},bgblocknorm);
                            end
                    end
                    
                end
            end
        end
    end
    
    for i = 1:length(p.names)
        switch p.method
            case 'pixel'
                rawback = accumulate_back{i}./back_count{i};
                rawback(isnan(rawback)) = 1;
                blurback = imgaussfilt(rawback,sigma);
                %blurback = imfilter(rawback,fspecial('disk',blursize),'symmetric');
                bias = blurback/median(blurback(:));
                
            case 'block'
                blockbias=nanmedian(biasstack{i},3);
                bias=imresize(blockbias,[height width],'bicubic');
                block_save{i} = cat(3,block_save{i}, biasstack{i});
        end
        save(fullfile(p.biaspath,[p.names{i},'_',num2str(site),'.mat']),'bias');
        bias_save{i,s} = bias;
    end
end

if p.biasAll
    for i = 1:length(p.names)
        switch p.method
            case 'pixel'
                bias = zeros(height,width);
                for j = 1:size(bias_save,2)
                    bias = bias + bias_save{i,j};
                end
                bias = bias/size(bias_save,2);
            case 'block'
                blockbias = nanmedian(block_save{i},3);
                bias=imresize(blockbias,[height width],'bicubic');
                bias_all{i} = bias;
        end
        for s = 1:numsites
            save(fullfile([p.biaspath,'_all'],[p.names{i},'_',num2str(p.site_mat(s)),'.mat']),'bias');
        end
    end
end

%%% Output bias
figure('pos',[10 10  numsites*100 length(p.names)*100])
for i = 1:length(p.names)
    for s = 1:numsites
        subplot(length(p.names),numsites,sub2ind([ numsites length(p.names)], s,i))
        outputsmall = imresize(bias_save{i,s},.25);
        imagesc(outputsmall);
        title([p.names{i} ' Site: ' num2str(p.site_mat(s))]);
    end
end
end
