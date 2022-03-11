function BiasCalc_combine_all_options
%%% file paths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imagepath='E:\4TB11\Data\';
shadingpath='F:\MATLAB\BGimages\IX_XL\';
experimentpath='C193-IF2\';
Biasdir=[imagepath,experimentpath,'Bias'];
%%% setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nucradius=12;%12;
names={
    'DAPI';
    'FarRed';
    };
row_mat = [2:5];
col_mat = [1];
site_mat = [1:37];
frame_mat=[1];
%%% Settings
blur_radius = 5;
method = 'block'; % or pixel
maskforeground = 0; %1 to mask nuclei
adaptive = [0 0 0 0];
dilate = {.5,.5,.5,.5}; %nucradius/2
foreground_calc = 0; % 1 for foreground calc
blocknum= 15;
prctile_thresh=[0 100];
sigma = 25;
compress = .25; 
prctilethresh = 50;



%%% Initialize variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numrows = numel(row_mat); 
numcols = numel(col_mat); 
numsites = numel(site_mat); 
numframes = numel(frame_mat);

load([shadingpath,'cmosoffset_bin1.mat'],'cmosoffset');
bgcmos = cmosoffset;
[height,width] = size(bgcmos);

if ~exist(Biasdir,'dir')
    mkdir(Biasdir);
end
if ~exist([Biasdir,'_all'],'dir')
    mkdir([Biasdir,'_all']);
end
bias_save = {};
block_save = {};
if strcmp(method,'block')
    for i = 1:length(names)
        block_save{i} = [];
    end
end


%% Calculate bias
for s = 1:numsites
    accumulate_back={};
    back_count = {};
    for i = 1:length(names)
        switch method
            case 'pixel'
                accumulate_back{i} = zeros(height,width);
                back_count{i} = zeros(height,width);
            case 'block'
                biasstack{i} = [];
        end
    end
    site = site_mat(s);
    for c = 1:numcols
        col = col_mat(c);
        for r = 1:numrows
            row = row_mat(r);
            for f = 1:numframes
                frame = frame_mat(f);
                shot = [num2str(row),'_',num2str(col),'_',num2str(site)];
                fprintf([shot,'\n']);
                rawdir = [imagepath,experimentpath,'Raw\',shot,'\',shot,'_'];
               
                %%% Loop through channels
                for i = 1:length(names)
                    try
                        raw{i} = double(imread([rawdir,names{i},'_',num2str(frame),'.tif']));
                    catch
                        warning([shot ' not found']);
                        continue;
                    end
                    raw{i} = raw{i} - bgcmos;
                    
                    %%% Mask Image
                    switch maskforeground
                        case 0
                            mask = false(size(raw{i}));
                        case 1
                            rawcomp = imresize(raw{i}, compress);
                            if adaptive(i)
                                normraw{i} = raw{i}/median(rawcomp(:));
                                overlay{i} = normraw{i}+normraw{1};
                                overlay{i}(overlay{i}<1) = 1;
                            else
                                overlay{i} = raw{1};
                            end
                            mask = threshmask(overlay{i},blur_radius);
                            if ~foreground_calc
                                mask = imdilate(mask, strel('disk', ceil(nucradius*dilate{i})));
                            else
                                mask = ~mask;
                                mask = imdilate(mask, strel('disk',ceil(nucradius*dilate{i})));
                            end
                    end
                    blur = imfilter(raw{i},fspecial('disk',blur_radius),'symmetric');
                    blurnan = blur;
                    blurnan(mask) = NaN;
                    %{
%%% debugging: view images %%%%%%%%%%
%%
 extractmask=bwmorph(mask,'remove');
[height,width]=size(mask);
RGB=zeros(height,width,3);
RGB(:,:,1)=imadjust(mat2gray(raw{i}));
RGB(:,:,2)=extractmask;
figure,imshow(RGB);
%%
                    %}
                    %%% Compile background
                    switch method
                        case 'pixel'
                            blurnancomp = imresize(blurnan, compress);
                            outliers = blurnan < prctile(blurnancomp(:),prctile_thresh(1)) | blurnan >prctile(blurnancomp(:),prctile_thresh(2));
                            mask = mask | outliers;
                            blurnan(mask) = NaN;
                            blurnancomp = imresize(blurnan,compress);
                            blur(mask)=0;
                            %                     imagesc(blur);
                            %                     title(names{i});
                            %                     pause;
                            blur = blur/nanmedian(blurnancomp(:));
                            accumulate_back{i} = accumulate_back{i} + blur;
                            back_count{i} = back_count{i} + ~mask;
                        case 'block'
                            bgblock=blockpercentile_blockimage(blurnan,blocknum,prctilethresh);
                            midrc=ceil(blocknum/2);
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
    
    for i = 1:length(names)
        switch method
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
        save([Biasdir,'\',names{i},'_',num2str(site),'.mat'],'bias');
        bias_save{i,s} = bias;
    end
end


for i = 1:length(names)
    switch method
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
        save([Biasdir,'_all\',names{i},'_',num2str(site_mat(s)),'.mat'],'bias');
    end
end

%%% Output bias
figure('pos',[10 10  numsites*100 length(names)*100])
for i = 1:length(names)
    for s = 1:numsites
        subplot(length(names),numsites,sub2ind([ numsites length(names)], s,i))
        outputsmall = imresize(bias_save{i,s},.25);
        imagesc(outputsmall);
        title([names{i} ' Site: ' num2str(site_mat(s))]);
    end
end
end
