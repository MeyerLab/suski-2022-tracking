function [STmat] = nd2Map(time,sites)
%ND2MAP returns the map of sites and time frames read in by bio-formats
%   Takes in the number of timepoints "time" and sites "sites"
%   Returns a matrix STmat{s,t} = [i,j], where s and t are the sites and
%   timpoints you need the index of, and i is the encoded "Series" and j is
%   the encoded "Time" from bio-formats


STmat = {};
tTemp = 1;
sTemp = 1;
for t = 1:time
    for s = 1:sites
        STmat{s,t} = [sTemp, tTemp];
        if tTemp + 1 <= time
            tTemp = tTemp + 1;
        else
            tTemp = 1;
            sTemp = sTemp + 1;
        end
    end
end

end

