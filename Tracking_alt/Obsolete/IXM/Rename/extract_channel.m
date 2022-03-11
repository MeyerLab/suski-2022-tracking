%% Name_change_IX_Micro_1site_3color_timelapse.m
% This script converts the filenames of image files generated using ix
% microXL into a more useful format
% The script assumes that image files are stored within subfolders named
% 'Timepoint_XX' of the folder 'dir'. The files are renamed and moved into
% 'dir'. The output format is 'Row_Column_Site_Channel_Timepoint.tif'
% Written by Steve Cappell
% 120919

%% Init
clear
close all
clc

%% Load Directory with images
original_dir='H:\Data\E1015-live\Raw';     %'/Users/scappell/Documents/Meyer_Lab/Data/2012-07-04/32';
destination_dir='H:\Data\E1015-live\Channels';    %'/Users/scappell/Documents/Meyer_Lab/Data/2012-07-04/32';

%% Designate Channels
channels={'CFP'};  %Change this to match the correct channels used in the experiment  %e.g. {'Channel1','Channel2','Channel3'}  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alldir=getSubdirectories(original_dir); %Creates list of all the subfolders

for i=1:length(alldir);
    disp(alldir(i))                                                                               %Displays what timepoint you are on
    filenames=getFilenames([original_dir,'/',char(alldir(i))]);                           %Returns list of all file names inside the folder
    
    for j=1:length(channels);
        ins=strfind(filenames,channels(j));
        list = find(not(cellfun('isempty', ins)));
        for k=1:length(list);
            newFileDir=[destination_dir,'/',char(channels(j)),'/',char(filenames(list(k)))];
            oldFileDir=[original_dir,'/',char(alldir(i)),'/',char(filenames(list(k)))];
            copyfile(oldFileDir,newFileDir,'f');
        end
    end
end
