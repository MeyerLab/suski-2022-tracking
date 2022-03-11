% addpath(genpath('/media/meyer/Nalin 8TB7/cell-cycle-tracking/Tracking_dependencies'))

%% Rename images
% nikon_rename({'/media/meyer/Nalin 8TB7/D114-live/20200713_103202_564',...
%     '/media/meyer/Nalin 8TB7/D114-live/20200713_143635_287'},'/media/meyer/Nalin 8TB7/D114-live/Raw/','COPY',false)


%% Calculate live image bias
clear p
%%%Directories
p.biaspath='/media/meyer/Nalin 8TB7/D114-live/Bias';     % Output path for bias
p.imagepath='/media/meyer/Nalin 8TB7/D114-live/Raw';     % Image directory. If nikon, single folder for all sites
p.shadingpath=''; %cmos offset file for IX micro
p.names={
    'CFP';
    'YFP';
    'RFP';
    };
p.row_mat = [2:7];
p.col_mat = [2:11];
p.site_mat = [1:9];
p.frame_mat=[1 50];
p.biasAll = 0;    % Save averaged bias across sites
                  % Formatting code for nikon files, 1$ is timepoint, 2$ is well, 3$ is point
p.formatCode = 'Time%1$05d_Well%2$s_Point%2$s_%3$04d_*.tiff'; 
p.microscope = 'nikon';

%%% Settings
p.maskforeground = 1;  % 1 to mask nuclei
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


% biasCalc(p,0)

%% Prep processing
clear s
%%% Paths
s.experiment_name='D114-live';                % Experiment name
s.image_drive = 'F:\D114-live';    % Location of images base folder
s.savepath=fullfile('E:\Data',s.experiment_name,'Data'); % Output location of mat files
s.imagepath=fullfile(s.image_drive,'Raw');                       % All sites in same folder for nikon
s.biaspath=fullfile(s.image_drive,'Bias');                       % Location of bias files
s.maskpath=fullfile(s.image_drive,'Mask');                       % Future location of mask saving
s.bgcmospath='';
% Formatting code for nikon files, 1$ is timepoint, 2$ is well, 3$ is point
s.formatCode = 'Time%1$05d_Well%2$s_Point%2$s_%3$04d_*.tiff'; 
s.microscope = 'nikon';

%%% General parameters
s.startFrame = 1; 
% s.endFrame = 114;
s.endFrame = getMaxFrame(fullfile(s.imagepath,'2_2_1'),'Time(\d+)_');

s.magnification=10;               % 10=10x or 20=20x
s.binsize=1;                      % 1=bin 1 or 2=bin 2
s.postbin = 0;             % 0 for no software bin, number for scaling factor (e.g. 0.5 for 2x2 bin);
s.signals={'CFP_','YFP_','RFP_'};
s.maskwrite = 1;                  % 1 = save an image with the nucmask; 0 = dont save an image with the nucmask
s.maskname = 'nucedge_';
s.register = 0;                   % Register images between timepoints for all
s.register_exception = [19 20 21 22];        % If don't register images, manually put in sites to register

%%% Quantification parameters
s.bgcmoscorrection = 1;
s.bias = [1 1 1];
s.signal_foreground = [0 0 0];
s.bgsubmethod = {'global nuclear','global nuclear','global nuclear'}; 
                                  % Options:'global nuclear','global cyto','tophat','semi-local nuclear', 'none'
s.compression = 4;                % For bg subtraction
s.bgprctile = [25  25 25];         % Set bg as perctile for each channel
s.sigblur = [3  3 3]; 
s.localbgmeasure = [0  0 0];       % Measure local background in box around cell
s.ringcalc = [0  1 1];
s.ringthresh = [0 0 0];         % Threshold for foreground ring pixels
s.punctacalc = 0;                 % Analyze puncta in nucleus
s.punctaThresh = [125 150 200];   % Threshold for tophat filter for puncta
s.varThresh = [75 100 125];       % Threshold for variance filter puncta 

%%% Segmentation parameters
s.firstsegmethod = 'log contour';  % Options: 'concavity', 'log', 'single adapt', 'double adapt', 'single', 'double'
s.secondsegmethod = 'log contour'; % Options: 'concavity', 'log', 'single adapt', 'double adapt', 'single', 'double'
s.nucr = 12;                     % 10x bin1 = 12 and 20x bin2 = 12
s.blurradius  =  3;              % 10x: 3
s.soliditythresh = 0.7;          % Minimum solidity of nucleus allowed
s.split_mult = 2;
s.debrisarea = 100;              % 10x = 100 % 150 if no mitosis
s.boulderarea = 1500;            % 10x = 1500
s.blobthreshold = -0.01;         % For blobdector 'log' segmentation

%%% Tracking parameters
s.maxjump = s.nucr*3;
s.masschangethreshold = 0.30;
s.areachangethreshold = 0.60;
s.daughtervariance = 0.10;

%%% Run timelapse on sample site
%Timelapse(s,2,2,4,0)

%% Run timelapse parallelized
rows = [2:7];
cols = [2:11];
sites = [1:9];
Parallelwell({@Timelapse},s,rows,cols,sites)

