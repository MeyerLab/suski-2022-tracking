function [S]=gathertracedata_NR_small(shot, settings)
%%% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([settings.dataDir,'tracedata_',shot,'.mat'],'tracedata','genealogy','jitters');
if IFoption
    load([datadir,IFlabel,shot,'.mat'],'IFdata');
else
    IFdata=ones(size(genealogy))*NaN;
end
%%% gate by genealogy %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gatemother=zeros(size(genealogy));
mothers=genealogy(~isnan(genealogy));
ismother = false(size(genealogy));
ismother(unique(mothers)) = true;
cellid=(1:length(genealogy))';
nonmothers=cellid(~ismember(cellid,mothers));
if motheroption==0 %no gating
    gatemother(:)=1;
elseif motheroption==1 %only traces that end in mitosis
    gatemother(mothers)=1;
elseif motheroption==2 %only traces that don't end in mitosis
    gatemother(nonmothers)=1;
end
if daughteroption==0
    gatedaughter=ones(size(genealogy)); %no gating
elseif daughteroption==1
    gatedaughter=~isnan(genealogy); %only traces that start with mitosis
elseif daughteroption==2
    gatedaughter=isnan(genealogy); %only traces that don't start with mitosis
end
gategenealogy=gatemother & gatedaughter;
%%% gate by presence of IF data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if IFoption
    gateIF=~isnan(IFdata(:,1));
else
    gateIF=ones(size(genealogy));
end
%%% combine gating %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
samplecells=gategenealogy & gateIF;
samplecellsID=find(samplecells);
%%% get full durations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
tracestats=getstats(tracedata,genealogy); %start,end,length,genealogy
catch
    keyboard
end
% tracedata=linkancestry(tracedata,tracestats,samplecellsID);
tracedata=linkancestry_NR(tracedata,tracestats,samplecellsID);
% tracedata=linkancestry_4generations_NR(tracedata,tracestats,samplecellsID);
%%% record mother stats %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if daughteroption==1
    motherstats=getmotherstatsonly(tracedata,tracestats,samplecellsID);
else
    motherstats=ones(size(genealogy))*NaN;
end
%%% remove gated out cells %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tracedata=tracedata(samplecells,:,:);
genealogy=genealogy(samplecells,:,:);
tracestats=tracestats(samplecells,:);
motherstats=motherstats(samplecells,:);
ismother = logical(ismother(samplecells));
allIFdata = IFdata;
unmatched = isnan(allIFdata(:,1));
allIFdata(unmatched,:) = [];

IFdata=IFdata(samplecells,:);
samplecellsID;
end