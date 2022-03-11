function maskwatershed=markershed_filter(mask,eroderadius, filtersize)
valleys=-bwdist(~mask);
basins=imerode(mask,strel('disk',eroderadius,0));
basins=bwareaopen(basins,filtersize);
bigvalleys=bwdist(mask);
outerridges=watershed(bigvalleys);
outerridges=outerridges==0;
finalvalleys=imimposemin(valleys,basins | outerridges);
finalridges=watershed(finalvalleys);
maskwatershed=mask;
maskwatershed(finalridges==0)=0;
end