clear all;
close all;
clc;

%% Designate wells to analyze %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rows =[2];
cols = [2:11];
sites= [1:16];

manualwells = [
    1 10 1;
    ];
manualcontrol = 0;

%%% Initialize variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numrows=length(rows);
numcols=length(cols); 
numsites=length(sites);
shots=numrows*numcols*numsites;
if manualcontrol==1
    shots=size(manualwells,1);
end

%% Run parallel processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time1=tic;
data_path = 'F:\Data\C-Cdt1\C163-live\Data\';
logFolder = [data_path char(datetime('now','Format','HHmmss-MMddyyyy')), '-logs\'];
%logFolder = 'F:\Data\C-Cdt1\C118-live\Data\111918-11012018-logs\';
parfor shot=1:shots
    %%% Calculate row/col/site for parallel worker
    if manualcontrol==1
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
    ID = t.ID;
    filePrefix = ['Worker' num2str(ID)];
    
    %% Run code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    try
      %Timelapse(row,col,site,0 )
      %Timelapse_addIF(row,col,site,0,logFolder,filePrefix)
      Fixed_IF(row,col,site,0 )
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
errorFiles = dir([logFolder,'Worker*error.txt']);
logFiles = dir([logFolder,'Worker*log.txt']);
errorAllID = fopen([logFolder,'all_error.txt'],'a');
logAllID = fopen([logFolder,'all_log.txt'],'a');

for i=1:length(errorFiles)
    fid = fopen([logFolder,errorFiles(i).name],'r');
    fwrite(errorAllID, fread(fid,'*char'),'*char');
    fclose(fid);
end


for i=1:length(logFiles)
    fid = fopen([logFolder,logFiles(i).name],'r');
    fwrite(logAllID, fread(fid,'*char'),'*char');
    fclose(fid);
end