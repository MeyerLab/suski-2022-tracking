function IX_Micro_rename_single_channel
% Define paths and settings
ORIGINAL_DIR = {'G:\Nalin\20190214-C147-live\2019-02-15\215', 'G:\Nalin\20190214-C147-live\2019-02-15\216'};
DESTINATION_DIR = 'G:\4TB5\Data\C147-live\Raw';
CHANNELS = {'CFP','YFP','RFP'};
LIVE = true;
COPY = false; % true for copy file, false for move file
SINGLESITE = false;
time1 = tic;  %Starts Timer

if LIVE
    %% Rename files for timelapse imaging
    frame = 1; % Running count of time point
    for d = 1:length(ORIGINAL_DIR)
        curr_dir  = ORIGINAL_DIR{d};
        sub_dirs = getSubdirectories(curr_dir); % Finds timepoint folders inside main
        if ~isempty(sub_dirs)
            times = getTokens(sub_dirs, 'TimePoint_(\d+)$'); % Finds timepoint folder numbers
            times = cellfun(@str2num, times);
            max_time = max(times);
            for t = 1:max_time
                disp(frame)
                time2 = tic;  % Timepoint timer
                file_names = getFilenames([curr_dir,'\TimePoint_',num2str(t)]); % List of files inside timepoint folder
                file_names = file_names(boolRegExp(file_names, '\.tif$') & ~boolRegExp(file_names, 'thumb')); % remove unneeded files
                if isempty(file_names) % Empty if folder is empty or missing
                    disp(['Time point ' num2str(frame) ' missing']);
                else
                    % Rename each file
                    for f = 1:length(file_names)
                        [row_col_site, channel] = getWellSite(file_names(f), SINGLESITE); % parse filename to get well and channel
                        if ~exist([DESTINATION_DIR,'\',row_col_site],'dir')
                            mkdir([DESTINATION_DIR,'\',row_col_site]);
                        end
                        old_file = [curr_dir,'\TimePoint_',num2str(t),'\',file_names{f}];
                        dest_file = [DESTINATION_DIR, '\', row_col_site, '\', row_col_site, '_', CHANNELS{channel}, '_', num2str(frame), '.tif'];
                        if ~exist(dest_file)
                            if COPY
                                copyfile(old_file,dest_file,'f');
                            else
                                movefile(old_file,dest_file,'f');
                            end
                        end
                    end
                end
                frame = frame + 1;
                toc(time2) % Displays time it takes to rename all the files in each Timepoint
            end
        end
    end
    
else
    %% Rename files for fixed imaging
    file_names = getFilenames(ORIGINAL_DIR{1}); % List of files inside  folder
    file_names = file_names(boolRegExp(file_names, '\.tif$') & ~boolRegExp(file_names, 'thumb')); % remove unneeded files4
    % Rename each file
    for f = 1:length(file_names)
        [row_col_site, channel] = getWellSite(file_names(f), SINGLESITE); % parse filename to get well and channel
        disp([row_col_site ' ' CHANNELS{channel}]);
        old_file = [ORIGINAL_DIR{1}, '\', file_names{f}];
        new_file = [DESTINATION_DIR, '\', row_col_site, '\', row_col_site, '_', CHANNELS{channel}, '_1', '.tif'];
        if ~exist(new_file)
            if ~exist([DESTINATION_DIR, '\', row_col_site], 'dir')
                mkdir([DESTINATION_DIR, '\', row_col_site]);
            end
            if COPY
                copyfile(old_file, new_file, 'f');
            else
                movefile(old_file, new_file, 'f');
            end
            
        end
    end
    
end
toc(time1)  %Displays total time to rename all files
end


function [row_col_site, channel] = getWellSite(name, SINGLESITE)
%Returns "row_col_site" and channel number
well_name = char(getTokens(name, '^[^_]+_([A-H]\d{2})_'));  %Finds well name in format wellColumn (eg B04)
row_column = wellNameToRowColumn(well_name);                            %Converts well name to format row_column (eg 2_4)

if ~SINGLESITE
    site = char(getTokens(name, '^[^_]+_[A-H]\d{2}_s(\d+)\w{8}-'));
    channel = getTokens(name, '^[^_]+_[A-H]\d{2}_s\d+_w(\d)'); % Channel number
else
    site = '1';
    channel = getTokens(name, '^[^_]+_[A-H]\d{2}+_w(\d)'); % Channel number
end
if ~isempty(channel{1})
    channel = str2num(channel{1});
else
    channel = 1;
end

row_col_site = [row_column, '_', site];
end