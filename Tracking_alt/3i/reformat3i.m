function reformat3i
paths = {'F:\Data\3i Data\PCNA 45\C176-3i-ind\site 1\', ...
    'F:\Data\3i Data\PCNA 45\C176-3i-ind\site 2\', ...
    'F:\Data\3i Data\PCNA 45\C176-3i-ind\site 3\'...
    'F:\Data\3i Data\PCNA 45\C176-3i-ind\site 4\'};
rows = [0:3];
cols = [0:3];
channels = {'445_50um_','515_50um_','561_50um_'};

outchannel = {'CFP_','YFP_','RFP_'};
outpath = 'F:\Data\3i Data\PCNA 45\C176-3i-ind\Raw\';
copy = true;
            

if ~exist(outpath)
    mkdir(outpath)
end

for p = 1:length(paths)
    path = paths{p};
    files = dir([path '*.tiff']);
    filenames = extractfield(files,'name');
    for r = 1:length(rows)
        row = rows(r);
        for c = 1:length(cols)
            col = cols(c);
            prefix = sprintf('%03u_%03u_%s',row,col);
            [~, numbers] = regexp(filenames, [prefix '\w*_(\d+).tiff'],'match','tokens');
            tokens = {};
            for i = 1:length(numbers)
                if ~isempty(numbers{i})
                    tokens = [tokens numbers{i}{1}];
                end
            end
            timepoints = string(unique(tokens));
            timepoints = sort(timepoints);
            
            dest_folder = [outpath num2str(p) '_' num2str(row+1) '_' num2str(col+1) '\'];
            if ~exist(dest_folder)
                mkdir(dest_folder)
            end
            
            for chan = 1:length(channels)
                channel = channels{chan};
                for t = 1:length(timepoints)
                    old_file = strcat(path,prefix, channel, timepoints(t), '.tiff');
                    new_file = strcat(dest_folder,num2str(p),'_',num2str(row+1),'_',num2str(col+1),'_', outchannel{chan}, num2str(t), '.tif');
                    if copy
                        copyfile(old_file,new_file,'f');
                    else
                        movefile(old_file,new_file,'f');
                    end
                end
            end
        end
    end
end
end
