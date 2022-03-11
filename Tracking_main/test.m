 %% Prep processing
%%% Paths
settings.data_path = 'F:\Data\C-Cdt1\C219-live\Data\';
settings.IF_imagesessions = {'H:\C219-IF2\'};
settings.live_imagepath = 'H:\C219-live\Raw';
settings.bgcmospath = '';
settings.crop_save = '';
settings.mask_save = '';           % Output mask, leave empty if no save


%%% General parameters
settings.microscope = 'nikon';
settings.formatCodeLive = 'Time%1$05d_Well%2$s_Point%2$s_%3$04d_*.tiff';
settings.formatCodeFixed = 'Time%1$05d_Well%2$s_Point%2$s_%3$04d_*.tiff';
settings.magnification = 20; %10 = 10x or 20 = 20x
settings.binsize = 1; %1 = bin 1 or 2 = bin 2
settings.postbin = .5; %0 for no bin, number for scaling factor;
settings.match10x20x = 1;
settings.numsites = 9;
settings.scaleLive = 1;
settings.maxjit = 500;
settings.signals = {{'DAPI_','YFP_','FarRed_'}};
settings.nucLive = 'CFP_';
settings.manualLiveFrame = [];

%%% Quantification parameters
settings.bias = {[1 1 1]};
settings.biasall = {[0 0 0]};
settings.sigblur = {[0 0 0]};
settings.signal_foreground = {[0 0 0]};
settings.bleedthrough = {[0 0 0]};
settings.bleedthroughslope = {};
settings.bleedthroughoff = {};
settings.ringcalc = {[1 1 1]};
settings.ringthresh = {[0 0 0 ]};
settings.punctacalc = {[0 1 0]};
settings.punctaThresh = {{[],[0:10:100],[]},{}};
settings.punctatopsize = 2;
settings.localbg = {[1 1 1]};
settings.minringsize = 100;
settings.bgcmoscorrection = 1; %1 = correct for shading; 0 = dont correct for shading;
settings.bgsubmethod = {{'global nuclear','global nuclear', 'global nuclear'}}; %'global nuclear','global cyto','tophat','semi-local nuclear'
settings.bgperctile = {[25 25 25]};
settings.frameIF=1;

%%% Segmentation parameters
settings.segmethod = {'concavity'}; %'log' or 'single' or 'double'
settings.nucr = 12; %10x bin1 = 12 and 20x bin2 = 12
settings.debrisarea = 100; %100
settings.boulderarea = 1500; %1500 % less restrictive for C152 for >4N cells
settings.blobthreshold = -0.03;
settings.blurradius = 3;
settings.soliditythresh = 0.50;
settings.compression = 4;
settings.badFrameCheck = .25; %0 for no check

%%% Tracking parameters
settings.distthresh = 3*settings.nucr;
settings.arealowthresh = -.5;
settings.areahighthresh = .6;

settings.cytopuncta = {[0 1 0]};
settings.thickenradius = 2*settings.nucr;

Timelapse_addIF_beta(settings,2,7,5,1)