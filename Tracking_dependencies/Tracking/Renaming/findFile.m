function filenames=findFile(path,search)
% searches for a single file based on name
d0=dir(fullfile(path,search));
filenames={d0([d0.isdir]==0).name};
filenames=filenames(:);

if length(filenames) > 1
    fprintf(['Error: cannot find single file based on ' search '\n']);
    filenames = [];
elseif length(filenames) == 0
    fprintf(['Error: cannot find any files based on ' search '\n']);
    filenames = [];
else
    filenames = filenames{1};
end

