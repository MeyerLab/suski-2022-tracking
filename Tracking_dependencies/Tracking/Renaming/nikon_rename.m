function nikon_rename(ORIGINAL_DIR, DESTINATION_DIR, varargin)
% Rename multiple part imaging sequence taken on nikon, put into site
% folders
% args:
% ORIGINAL_DIR - either single path or cell array of paths of files which
% are to be renamed. If cell array, will continue timepoint counter
% DESTINATION_DIR - folder to copy to
% optional args:
% COPY - if true, keep original file
% START_FRAME - manual start frame
% SINGLE_SITE - For stitched images


%% set/fix parameters 
if ~iscell(ORIGINAL_DIR)
    ORIGINAL_DIR = {ORIGINAL_DIR};
end

%detect parameters
p = inputParser;
addParameter(p,'COPY',true,@islogical);
addParameter(p,'START_FRAME',0,@isnumeric);
addParameter(p,'SINGLE_SITE',false,@islogical);

parse(p,varargin{:});
COPY = p.Results.COPY;
PREV_FRAME = p.Results.START_FRAME - 1;
SINGLE_SITE = p.Results.SINGLE_SITE;


time1 = tic;  %Starts Timer
if ~exist(DESTINATION_DIR,'dir')
    mkdir(DESTINATION_DIR);
end
%% Rename files for timelapse imaging
for d = 1:length(ORIGINAL_DIR)
    time2 = tic;  % Timepoint timer
    curr_dir  = ORIGINAL_DIR{d};
    files = getFilenames(curr_dir); % Finds timepoint folders inside main
    if ~isempty(files)
        currFrames = getTokens(files, 'Time(\d+)_'); % Finds timepoint numbers
        currFrames = cellfun(@str2num, currFrames);
        currWell = getTokens(files,'Well([A-H]\d{2})_');
        if ~SINGLE_SITE
            currSite = cellfun(@str2num,getTokens(files,'Point[A-H]\d{2}_(\d+)'));
        else
            currSite = zeros(size(currFrames));
        end
        maxFrameCurr = max(currFrames);
        renameFrame = currFrames + PREV_FRAME+1;
        for f = 1:length(files)
            oldFile = files{f};
            newTime = sprintf('Time%05d_',renameFrame(f));
            newFile = regexprep(oldFile,'Time(\d+)_',newTime);
            newFile = regexprep(newFile,'.tiff','_renamed.tiff');
            if strcmp(oldFile,newFile)
                disp(['Error with file renaming: ' oldFile]);
                continue;
            end
            %%% Rename File
            fprintf('Dir %d: Renaming frame %d to frame %d\n',d,currFrames(f),renameFrame(f));
            oldPath = fullfile(curr_dir,oldFile);
            new_dir = fullfile(DESTINATION_DIR,[wellNameToRowColumn(currWell{f}) '_' num2str(currSite(f)+1)]);
            if ~exist(new_dir,'dir')
                mkdir(new_dir);
            end
            destPath = fullfile(new_dir,newFile);
            if ~exist(destPath)
                if COPY
                    copyfile(oldPath,destPath,'f');
                else
                    movefile(oldPath,destPath,'f');
                end
            else
                disp([destPath ' already exists']);
            end
        end
    else
        disp(['Cannot find files in ' curr_dir]);
    end
    PREV_FRAME = max(renameFrame);
    disp(['Directory took ' num2str(toc(time1)) ' sec']); % Displays time it takes to rename all the files in each Timepoint
end

disp(['Total elapsed time: ' num2str(toc(time2)) ' sec']); % Displays total time to rename all files
end
