function maxFrame = getMaxFrame(DESTINATION_DIR,code)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
filesDest = getFilenames(DESTINATION_DIR);
maxFrame = getTokens(filesDest, code);
if ~isempty(maxFrame)
    maxFrame = max(cellfun(@str2num,maxFrame))+1;
else
    maxFrame = -1;
end

end

