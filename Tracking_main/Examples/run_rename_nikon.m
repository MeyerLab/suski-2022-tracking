%% Calculate Bias live
% clear s
% s.biaspath='J:\C196_two_stage\Bias';
% s.imagepath='J:\C196_two_stage\Raw1';
% s.shadingpath='';
% s.nucradius=24;%12;
% s.names={
%     'CFP';
%     'RFP';
%     };
% s.formatCode = 'Time%1$05d_Well%2$s_Point%2$s_%3$04d_*.tiff';
% s.microscope = 'nikon';
% s.row_mat = [6];
% s.col_mat = [2:11];
% s.site_mat = [1:16];
% s.frame_mat=[10];
% %%% Settings
% s.blur_radius = 5;
% s.method = 'block'; % or pixel
% s.maskforeground = 1; %1 to mask nuclei
% s.adaptive = [1 1 0];
% s.dilate = {.75,1.25,.5}; %nucradius/2
% s.foreground_calc = 0; % 1 for foreground calc
% s.blocknum= 9;
% s.prctilethresh = 50;
% s.compress = .25;
% s.prctile_thresh=[0 100]; %pixel
% s.sigma = 25; % pixel
% 
% biasCalc(s,0);

%% Calculate Bias IF1
% clear s
% s.biaspath='J:\Nikon\C196_two_stage_IF1\Bias';
% s.imagepath='J:\Nikon\C196_two_stage_IF1\Raw';
% s.shadingpath='';
% s.nucradius=24;%12;
% s.names={
%     'CFP';
%     'YFP';
%     'RFP';
%     };
% s.formatCode = 'Time%1$05d_Well%2$s_Point%2$s_%3$04d_*.tiff';
% s.microscope = 'nikon';
% s.row_mat = [2:5];
% s.col_mat = [1];
% s.site_mat = [1:16];
% s.frame_mat=[1];
% %%% Settings
% s.blur_radius = 5;
% s.method = 'block'; % or pixel
% s.maskforeground = 0; %1 to mask nuclei
% s.adaptive = [1 1 0];
% s.dilate = {.75,1.25,.5}; %nucradius/2
% s.foreground_calc = 0; % 1 for foreground calc
% s.blocknum= 15;
% s.prctilethresh = 50;
% s.compress = .25;
% s.prctile_thresh=[0 100]; %pixel
% s.sigma = 25; % pixel
% 
% biasCalc(s,0);

%% Calculate Bias IF1
% clear s
% s.biaspath='J:\Nikon\C196_two_stage_IF2\Bias';
% s.imagepath='J:\Nikon\C196_two_stage_IF2\Raw';
% s.shadingpath='';
% s.nucradius=24;%12;
% s.names={
%     'DAPI';
%     'YFP';
%     'FarRed';
%     };
% s.formatCode = 'Time%1$05d_Well%2$s_Point%2$s_%3$04d_*.tiff';
% s.microscope = 'nikon';
% s.row_mat = [2:5];
% s.col_mat = [1];
% s.site_mat = [1:16];
% s.frame_mat=[1];
% %%% Settings
% s.blur_radius = 5;
% s.method = 'block'; % or pixel
% s.maskforeground = 0; %1 to mask nuclei
% s.adaptive = [1 1 0];
% s.dilate = {.75,1.25,.5}; %nucradius/2
% s.foreground_calc = 0; % 1 for foreground calc
% s.blocknum= 15;
% s.prctilethresh = 50;
% s.compress = .25;
% s.prctile_thresh=[0 100]; %pixel
% s.sigma = 25; % pixel
% 
% biasCalc(s,0);

%% Process
%Parallelwell



%% Prep processing
%%% Paths
settings.data_path = 'F:\Data\C-Cdt1\C196-live\Data_old\';
settings.IF_imagesessions = {'J:\Nikon\C196_two_stage_IF2\', 'J:\Nikon\C196_two_stage_IF1\'};
settings.live_imagepath = 'J:\Nikon\C196_two_stage\Raw2\';
settings.bgcmospath = '';
settings.crop_save = 'J:\Nikon\C196_two_stage-IFcrop\';


%%% General parameters
settings.microscope = 'nikon';
settings.formatCodeLive = 'Time%1$05d_Well%2$s_Point%2$s_%3$04d_*.tiff';
settings.formatCodeFixed = 'Time%1$05d_Well%2$s_Point%2$s_%3$04d_*.tiff';
settings.magnification = 20; %10 = 10x or 20 = 20x
settings.binsize = 1; %1 = bin 1 or 2 = bin 2
settings.postbin = 0; %0 for no bin, number for scaling factor;
settings.match10x20x = 0;
settings.numsites = [];
settings.scaleLive = 1;
settings.signals = {{'DAPI_','YFP_','FarRed_'}, {'CFP_','YFP_','RFP_'}};
settings.nucLive = 'CFP_';
settings.manualLiveFrame = [84];

%%% Quantification parameters
settings.bias = {[1 1 1], [1 1 1]};
settings.biasall = {[0 0 0], [0 0 0]};
settings.sigblur = {[0 0 0], [0 0 0]};
settings.signal_foreground = {[0 0 0], [0 0 0]};
settings.bleedthrough = {[0 0 0], [0 0 0]};
settings.bleedthroughslope = {};
settings.bleedthroughoff = {};
settings.ringcalc = {[1 1 1], [0 0 1]};
settings.ringthresh = {[0 0 0 ], [0 0 0]};
settings.punctacalc = {[0 0 0], [0 0 0]};
settings.punctaThresh = {{[],[],[]},{}};
settings.punctatopsize = 2;
settings.localbg = {[1 1 1], [0 0 1]};
settings.minringsize = 100;
settings.bgcmoscorrection = 1; %1 = correct for shading; 0 = dont correct for shading;
settings.bgsubmethod = {{'global nuclear','global nuclear', 'global nuclear'}, ...
    {'global nuclear', 'global nuclear', 'global nuclear'}}; %'global nuclear','global cyto','tophat','semi-local nuclear'
settings.bgperctile = {[25 25 25], [25 25 25]};
settings.frameIF=1;

%%% Segmentation parameters
settings.segmethod = {'concavity','thresh'}; %'log' or 'single' or 'double'
settings.nucr = 24; %10x bin1 = 12 and 20x bin2 = 12
settings.debrisarea = 400; %100
settings.boulderarea = 4500; %1500 % less restrictive for C152 for >4N cells
settings.blobthreshold = -0.03;
settings.blurradius = 3;
settings.soliditythresh = 0.80;
settings.compression = 4;
settings.badFrameCheck = .25; %0 for no check

%%% Tracking parameters
settings.distthresh = 3*settings.nucr;
settings.arealowthresh = -.4;
settings.areahighthresh = .4;

Timelapse_addIF(settings,6,9,4,0)


%% Run timelapse
% rows = [5:6];
% cols = [2:11];
% sites = [1:16];
% Parallelwell({@Timelapse_addIF},settings,rows,cols,sites)
