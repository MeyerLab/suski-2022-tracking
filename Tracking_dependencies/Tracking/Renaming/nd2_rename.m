function nd2_rename(ORIGINAL_DIR,DESTINATION_DIR,CHANNELS,manualwells, workers)
% Define paths and settings
% manualwells = {};
% ORIGINAL_DIR = {'J:\C185_PCNA_45_fast15\20190607_123645_700\'};
% DESTINATION_DIR = 'J:\C185-live-nikon\Raw';
% CHANNELS = {'CFP','YFP','RFP'};
% workers = 2;
% delete(gcp('nocreate'))
% parpool(workers);
% nikon_rename(ORIGINAL_DIR,DESTINATION_DIR,CHANNELS,manualwells, workers)
% delete(gcp('nocreate'))

time1 = tic;  %Starts Timer
%% Rename files for timelapse imaging
totalFrames = 1; % Running count of time point
for d = 1:length(ORIGINAL_DIR)
    curr_dir  = ORIGINAL_DIR{d};
    files = dir([curr_dir '*.nd2']);
    filenames = extractfield(files,'name');    
    %for w = 1:length(filenames)
    parfor (w = 1:length(filenames),workers)
        well_name = getTokens(filenames{w}, 'Well([A-H]\d{2})_'); % Finds timepoint folder numbers
        row_column = wellNameToRowColumn(well_name);
        if isempty(manualwells) | any(ismember(manualwells,row_column))
            disp(row_column);
            path = [curr_dir filenames{w}];
%             reader = loci.formats.Memoizer(bfGetReader(),0);
%             reader.setId(path);
            reader = bfGetReader(path);
            numsites = reader.getSeriesCount();
            numchannels = reader.getSizeC();
            numtime = reader.getSizeT();
            
            STmat = nd2Map(numtime,numsites);
            iZ = 1;
            for s = 1:numsites
                tempFrames = totalFrames;
                row_col_site = [row_column '_' num2str(s)];
                if ~exist(fullfile(DESTINATION_DIR,row_col_site),'dir')
                    mkdir(fullfile(DESTINATION_DIR,row_col_site));
                end
                for t = 1:numtime
                    disp([row_col_site '_' num2str(tempFrames)]);
                    iT = STmat{s,t}(2);
                    iSeries = STmat{s,t}(1);
                    reader.setSeries(iSeries -1);
                    for c = 1:numchannels
                        iC = c;
                        dest_file = fullfile(DESTINATION_DIR, row_col_site, [row_col_site, '_', CHANNELS{c}, '_', num2str(tempFrames), '.tif']);
                        if ~exist(dest_file,'file')
                            iPlane = reader.getIndex(iZ - 1, iC -1, iT - 1) + 1;
                            raw = bfGetPlane(reader,iPlane);
                            imwrite(raw,dest_file);
                        end
                    end
                    tempFrames = tempFrames + 1;
                end
            end
        end
    end
    path = [curr_dir filenames{1}];
    reader = bfGetReader(path);
    numtime = reader.getSizeT();
    totalFrames = totalFrames + numtime;
end

toc(time1)  %Displays total time to rename all files
end
