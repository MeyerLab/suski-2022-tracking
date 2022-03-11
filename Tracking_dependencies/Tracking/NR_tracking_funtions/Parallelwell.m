function Parallelwell(funcIn, settings, rows, cols, sites, varargin)
% Run function handle funcIn in parallel for the given settings, rows,
% cols, sites. if manualwells selected (cell array of row_col_site), put as
% varargin with 'manualwells'. Also selects path for error logs in current
% directory (usually where running processing functions). If enter
% parameter 'logPath', saves log files there.

%%% Manual well control
p = inputParser;
addParameter(p,'manualwells',{});
addParameter(p,'logPath',cd);
parse(p,varargin{:});
manualwells = p.Results.manualwells;
logPath = p.Results.logPath;



%%% Initialize variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numrows=length(rows);
numcols=length(cols); 
numsites=length(sites);
shots=numrows*numcols*numsites;
if ~isempty(manualwells)
    shots=size(manualwells,1);
end

%%% log folders
logFolder = {};
for i = 1:length(funcIn)
    logFolder{i} = fullfile(logPath,'error_logs',[char(funcIn{i}) '_' char(datetime('now','Format','HHmmss-MMddyyyy')), '_logs\']);
end
        
%% Run parallel processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time1=tic;
parfor shot=1:shots
    %%% Calculate row/col/site for parallel worker
    if ~isempty(manualwells)
        row=manualwells(shot,1);
        col=manualwells(shot,2);
        site=manualwells(shot,3);
    else
        siteidx=mod(shot,numsites);
        if siteidx==0
            siteidx=numsites;
        end
        site=sites(siteidx);
        colidx=mod(ceil(shot/numsites),numcols);
        if colidx==0
            colidx=numcols;
        end
        col=cols(colidx);
        rowidx=ceil(shot/(numcols*numsites));
        row=rows(rowidx);
    end
    fprintf('Shot %02d_%02d_%02d started\n',row,col,site);
    time2 = tic;
    
    %%% Get worker ID
    t = getCurrentTask();
    if isempty(t)
        ID = 0;
    else
        ID = t.ID;
    end
    filePrefix = ['Worker' num2str(ID)];
    
    %% Run code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    try
        for i = 1:length(funcIn)
            funcIn{i}(settings,row,col,site,0,logFolder{i},filePrefix);
        end
    catch ME
        disp(['Error: ',num2str(row),'_',num2str(col),'_',num2str(site),' Worker ' num2str(ID)]);
        fprintf('%s \n', ME.message);
    end
    elapsedTime2 = toc(time2);
    fprintf('Shot %02d_%02d_%02d finished after %07.2f min\n',row,col,site,elapsedTime2/60);
end

elapsedTime1 = toc(time1);
fprintf('Total elapsed time: %.1f min \n', elapsedTime1/60);
fclose all;

%% Combine log files %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(funcIn)
    errorFiles = dir([logFolder{i},'Worker*error.txt']);
    logFiles = dir([logFolder{i},'Worker*log.txt']);
    errorAllID = fopen([logFolder{i},'all_error.txt'],'a');
    logAllID = fopen([logFolder{i},'all_log.txt'],'a');
    
    for j=1:length(errorFiles)
        fid = fopen([logFolder{i},errorFiles(j).name],'r');
        fwrite(errorAllID, fread(fid,'*char'),'*char');
        fclose(fid);
    end
    
    
    for j=1:length(logFiles)
        fid = fopen([logFolder{i},logFiles(j).name],'r');
        fwrite(logAllID, fread(fid,'*char'),'*char');
        fclose(fid);
    end
end