% nikon_rename('J:\D112-IF\20200707_150240_897','J:\D112-IF\Raw')

%% Calculate Bias IF1 (empty well)
clear p
%%%Directories
p.biaspath='K:\4TB10\Data\A48-IF1\Bias';                       % Output path for bias
p.imagepath='K:\4TB10\Data\A48-IF1\Raw';                       % Image directory. If nikon, single folder for all sites
p.shadingpath=''; % cmos offset file for IX micro
p.names={
    'DAPI';
    'FarRed';

    };
p.row_mat = [2:4];
p.col_mat = [1];
p.site_mat = [1:37];
p.frame_mat=[1];
p.biasAll = 1;  % Save averaged bias across sites
                % Formatting code for nikon files, 1$ is timepoint, 2$ is well, 3$ is point
p.formatCode = 'Time%1$05d_Well%2$s_Point%2$s_%3$04d_*.tiff'; 
p.microscope = 'IXM';

%%% Settings
p.maskforeground = 0;  % 1 to mask nuclei
p.nucradius=12;%12;
p.blur_radius = 5;
p.adaptive = [0 0 0];  % use add channels together to segment foreground
p.dilate = {.75,1,1};  % multiples of nucradius, default .5, to expand mask
p.foreground_calc = 0; % 1 for foreground calc (average foreground signal)
p.method = 'block';    % or pixel
p.blocknum= 15;
p.prctilethresh = 50;
p.compress = .25;
p.prctile_thresh=[0 100]; % pixel remove outliers
p.sigma = 25;             % pixel smooth window

% biasCalc(p,0);


 %% Prep processing Fixed
 clear settings
%%% Paths
settings.data_path = 'F:\Data\A-Optimization\A49-IF\';
settings.IF_imagesessions = {'K:\4TB10\Data\A48-IF1\'};
settings.bgcmospath = 'C:\Users\nalin\Documents\GitHub\cell-cycle-tracking\BGimages\IX_Liu\cmosoffset_bin1.mat';
settings.crop_save = '';
settings.mask_save = 'K:\4TB10\Data\A48-IF1\Raw';           % Output mask, leave empty if no save
settings.maskname = 'mask';
settings.primaryMaskRound = 1;
settings.maskIndex = [1];

%%% General parameters
settings.microscope = 'IXM';
settings.formatCodeLive = 'Time%1$05d_Well%2$s_Point%2$s_%3$04d_*.tiff';
settings.formatCodeFixed = 'Time%1$05d_Well%2$s_Point%2$s_%3$04d_*.tiff';
settings.magnification = 20; %10 = 10x or 20 = 20x
settings.binsize = 1; %1 = bin 1 or 2 = bin 2
settings.postbin = .5; %0 for no bin, number for scaling factor;
settings.signals = {{'DAPI','FarRed'}};
settings.maxjit = 500;

%%% Quantification parameters
settings.bias = {[1 1 ]};
settings.biasall = {[1 1]};
settings.sigblur = {[0 0]};
settings.signal_foreground = {[0 0 ]};
settings.bleedthrough = {[0 0 ]};
settings.bleedthroughslope = {};
settings.bleedthroughoff = {};
settings.ringcalc = {[1 1 ]};
settings.ringthresh = {[0 0]};
settings.punctacalc = {[0 0 ]};
settings.punctaThresh = {{[],[],[]},{}};
settings.punctatopsize = 2;
settings.localbg = {[1 1 ]};
settings.minringsize = 100;
settings.bgcmoscorrection = 1; %1 = correct for shading; 0 = dont correct for shading;
settings.bgsubmethod = {{'global nuclear','none'}}; %'global nuclear','global cyto','tophat','semi-local nuclear'
settings.bgperctile = {[25 25]};
settings.frameIF=1;

%%% Segmentation parameters
settings.segmethod = {'concavity'}; %'log' or 'single' or 'double'
settings.nucr = 12; %10x bin1 = 12 and 20x bin2 = 12
settings.debrisarea = 100; %100
settings.boulderarea = 4500; %1500 % less restrictive for C152 for >4N cells
settings.blobthreshold = -0.03;
settings.blurradius = 3;
settings.soliditythresh = 0.50;
settings.compression = 4;
settings.badFrameCheck = .25; %0 for no check
settings.split_mult = 1;

Fixed_IF(settings,2,2,1,1)


%% Run add IF
rows = [2:7];
cols = [2:11];
sites = [1:37];
Parallelwell({@Fixed_IF},settings,rows,cols,sites)
