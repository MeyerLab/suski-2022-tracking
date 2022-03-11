function nikon_rename_time(ORIGINAL_DIR, DESTINATION_DIR, varargin)
% Rename multiple part imaging sequence taken on nikon
% args:
% ORIGINAL_DIR - either single path or cell array of paths of files which
% are to be renamed. If cell array, will continue timepoint counter
% DESTINATION_DIR - folder to copy to, if has images, renamed timepoint
% will start at the next timepoint
% optional args:
% COPY - if true, keep original file
% PREV_FRAME - manual previous frame (renamed files will start at
% PREV_FRAME+1)


%% set/fix parameters 
if ~iscell(ORIGINAL_DIR)
    ORIGINAL_DIR = {ORIGINAL_DIR};
end

% Detect max frame of destination dir
filesDest = getFilenames(DESTINATION_DIR);
filesDest = filesDest(~contains(filesDest,'renamed.tiff'));
maxFrameDest = getTokens(filesDest, 'Time(\d+)_');
if ~isempty(maxFrameDest)
    maxFrameDest = max(cellfun(@str2num,maxFrameDest));
else
    maxFrameDest = -1;
end
%detect parameters
p = inputParser;
addParameter(p,'COPY',true,@islogical);
addParameter(p,'PREV_FRAME',maxFrameDest,@isnumeric);

parse(p,varargin{:});
COPY = p.Results.COPY;
PREV_FRAME = p.Results.PREV_FRAME;


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
        maxFrameCurr = max(currFrames);
        renameFrame = currFrames + PREV_FRAME+1;
        for f = 1:length(files)
            oldFile = files{f};
            newTime = sprintf('Time%05d_',renameFrame(f));
            newFile = regexprep(oldFile,'Time(\d+)_',newTime);
            newFile = regexprep(newFile,'.tiff','_renamed.tiff');
            if strcmp(oldFile,newFile)
                disp(['Error with file: ' oldFile]);
                continue;
            end
            %%% Rename File
            fprintf('Dir %d: Renaming frame %d to frame %d\n',d,currFrames(f),renameFrame(f));
            oldPath = fullfile(curr_dir,oldFile);
            destPath = fullfile(DESTINATION_DIR,newFile);
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
