%% Calculate live image bias us old bias for this run since hard to do DHB



%% Prep live single frame processing 
%%% Paths
s.experiment_name='C201-live';                     % Experiment name
s.image_drive = 'E:\Nikon\C201_live\';   % Location of images base folder
s.savepath=fullfile('F:\Data\C-Cdt1\',s.experiment_name,'Data'); % Output location of mat files
s.imagepath=fullfile(s.image_drive,'Raw');                       % All sites in same folder for nikon
s.biaspath=fullfile(s.image_drive,'Bias');                       % Location of bias files
s.maskpath=fullfile(s.image_drive,'Mask');                       % Future location of mask saving
s.bgcmospath='';
% Formatting code for nikon files, 1$ is timepoint, 2$ is well, 3$ is point
s.formatCode = 'Time%1$05d_Well%2$s_Point%2$s_%3$04d_*.tiff';
s.microscope = 'nikon';

%%% General parameters
s.frame = 1; 
s.magnification=10;              % 10=10x or 20=20x
s.binsize=1;                     % 1=bin 1 or 2=bin 2
s.postbin=0;                     % 0 for no software bin, number for scaling factor (e.g. 0.5 for 2x2 bin);
s.signals={'CFP_','YFP_','RFP_'};
s.maskwrite = 1;                 % 1 = save an image with the nucmask; 0 = dont save an image with the nucmask
s.maskname = 'nucedge_';


%%% Quantification parameters
s.bgcmoscorrection = 1;
s.bias = [1 1 1];
s.signal_foreground = [0 1 0];
s.bgsubmethod = {'global nuclear','global nuclear','global nuclear'}; 
                                 % Options:'global nuclear','global cyto','tophat','semi-local nuclear', 'none'
s.compression = 4;               % For bg subtraction
s.bgprctile = [25 25 25];        % Set bg as perctile for each channel
s.sigblur = [3 3 3];
s.localbgmeasure = [0 0 0];      % Measure local background in box around cell
s.ringcalc = [0 1 0];
s.ringthresh = [0 0 0];            % Threshold for foreground ring pixels
s.punctacalc = 0;                % Analyze puncta in nucleus
s.punctaThresh = [125 150 200];
s.varThresh = [75 100 125];

%%% Segmentation parameters
s.firstsegmethod = 'concavity';  % Options: 'concavity', 'log', 'single adapt', 'double adapt', 'single', 'double'
s.nucr = 12;                     % 10x bin1 = 12 and 20x bin2 = 12
s.blurradius  =  3;              % 10x: 3
s.soliditythresh = 0.8;          % Minimum solidity of nucleus allowed
s.debrisarea = 100;              % 10x = 100 % 150 if no mitosis
s.boulderarea = 1500;            % 10x = 1500
s.blobthreshold = -0.02;         % For blobdector 'log' segmentation


%% Run timelapse on sample well
% Timelapse_singleFrame(s,6,3,3,1)

%%% Run timelapse parallelized
% rows = [2:7];
% cols = [2:9];
% sites = [1:9];
% Parallelwell({@Timelapse_singleFrame},s,rows,cols,sites)

%% Calculate Bias IF (empty well)
clear p
%%%Directories
p.biaspath='E:\Nikon\C201-IF1\Bias';                       % Output path for bias
p.imagepath='E:\Nikon\C201-IF1\Raw';                       % Image directory. If nikon, single folder for all sites
p.shadingpath=''; % cmos offset file for IX micro
p.names={
    'DAPI';
    'FarRed';
    };
p.row_mat = [2:5];
p.col_mat = [1];
p.site_mat = [1:9];
p.frame_mat=[1];
p.biasAll = 0;  % Save averaged bias across sites
                % Formatting code for nikon files, 1$ is timepoint, 2$ is well, 3$ is point
p.formatCode = 'Time%1$05d_Well%2$s_Point%2$s_%3$04d_*.tiff'; 
p.microscope = 'nikon';

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

biasCalc(p,0);


%% Prep fixed IF to live processing
%%% Paths
sIF.data_path = 'F:\Data\C-Cdt1\C201-live\Data\';      % Output location of mat files
sIF.IF_imagesessions = {'E:\Nikon\C201-IF1\', 'E:\Nikon\C201-IF2\'};    
                                                            % Cell array containing all imaging sessions to be matched
sIF.live_imagepath = 'E:\Nikon\C201_live\Raw\';    % Live imaging raw image path
sIF.bgcmospath = '';
sIF.crop_save = 'E:\Nikon\C201-IFcrop\';           % Output cropped images, leave empty if no save
sIF.mask_save = 'E:\Nikon\C201-IFcrop\';           % Output mask, leave empty if no save

%%% General parameters
sIF.microscope = 'nikon';
% Formatting code for nikon files, 1$ is timepoint, 2$ is well, 3$ is point
sIF.formatCodeLive = 'Time%1$05d_Well%2$s_Point%2$s_%3$04d_*.tiff';
sIF.formatCodeFixed = 'Time%1$05d_Well%2$s_Point%2$s_%3$04d_*.tiff';
sIF.magnification = 10;     % 10 = 10x or 20 = 20x
sIF.binsize = 1;            % 1 = bin 1 or 2 = bin 2
sIF.postbin = 0;           % 0 for no software bin, number for scaling factor (e.g. 0.5 for 2x2 bin);
sIF.match10x20x = 0;
sIF.numsites = 4;         % Number of sites imaged for 20x to 10x (double check order is correct)
sIF.scaleLive = 1;          % Scale live to have same pixel size as final IF image
sIF.maxjit = 500;      % maximum jitter between IF images before throwing warning
sIF.signals = {{'DAPI_','FarRed_'},{'DAPI_'}};  % Each imaging session in a cell inside cell array
sIF.nucLive = 'CFP_';       % Label for IX micro live nucleus (assume first tiff in stack for nikon)
sIF.manualLiveFrame = [];   % Manually set last frame for live imaging, otherwise use tracedata if empty

%%% Quantification parameters
sIF.bias = {[1 1 ], [1]};
sIF.biasall = {[0 0 ] , [0]};
sIF.sigblur = {[0 0 ], [0]};
sIF.signal_foreground = {[0 0 ], [0]};
sIF.bleedthrough = {[0 0 ], [0]};   % Calculate bleedthrough for channel
sIF.bleedthroughslope = {};              % Cell array of cell arrays
sIF.bleedthroughoff = {};                % Cell array of cell arrays
sIF.ringcalc = {[1 1 ], [1]};
sIF.ringthresh = {[0 0 ], [0]};    % Threshold for foreground ring pixels
sIF.punctacalc = {[0 0 ], [0]};     % Analyze puncta in nucleus
sIF.punctaThresh = {{[],[]},{}};      % Threshold for tophat filter for puncta
sIF.punctatopsize = 2;                   % Top hat filter size
sIF.localbg = {[0 0 ], [0]};        % Measure local background in box around cell
sIF.minringsize = 100;
sIF.bgcmoscorrection = 1;                % 1 = correct for shading; 0 = dont correct for shading;
sIF.compression = 4;
sIF.bgsubmethod = {{'global nuclear','global nuclear'},{'global nuclear'}}; % 'global nuclear','global cyto','tophat','semi-local nuclear'
sIF.bgperctile = {[25 25], [25]};  % Set bg as perctile for each channel
sIF.frameIF=1;
 
%%% Segmentation parameters
sIF.segmethod = {'concavity','thresh'};  % 'log' or 'single' or 'double'
sIF.nucr = 12;                           % 10x bin1 = 12 and 20x bin2 = 12
sIF.debrisarea = 100;                    % 100
sIF.boulderarea = 1500;                  % 1500 
sIF.blobthreshold = -0.03;               % For blobdector 'log' segmentation
sIF.blurradius = 3;                      % Blur for blobdetector
sIF.soliditythresh = 0.80;               % Minimum solidity of nucleus allowed
sIF.badFrameCheck = .25;                 % Min fraction of cells untracked before logging error, 0 for no check

%%% Tracking parameters
sIF.distthresh = 3*sIF.nucr;
sIF.arealowthresh = -.4;
sIF.areahighthresh = .5;

%% Run add IF on sample well
Timelapse_addIF(sIF,7,2,9,1)


%%% Run timelapse parallelized
% rows = [2:7];
% cols = [2:9];
% sites = [1:9];
% Parallelwell({@Timelapse_addIF},sIF,rows,cols,sites)
