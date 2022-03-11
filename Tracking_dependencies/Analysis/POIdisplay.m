function POIdisplay(traceIDs,sigstore, firstframe, lastframe, badtraces,ylims,POIs)
numvar=length(POIs); %number of different POIs to mark
for i=1:numel(traceIDs)
    id = traceIDs(i);
    startFrame = firstframe(id);
    endFrame = lastframe(id);
    
    if isnan(startFrame) | isnan(endFrame)
        continue;
    end
    
    figure(ceil(i/24)+1000); set(gcf,'color','w');
    subaxis(4,6,mod(i-1,24)+1,'ML',0.05,'MR',0.02,'MT',0.03,'MB',0.05,'SH',0.03);
    hold on;
    frames=startFrame:endFrame;
    plot(1:length(frames),sigstore(id,frames));
    ylim(ylims);
    
    title(num2str(id));

    if badtraces(traceIDs(i))==1
        midx=round(length(frames)/2);
        plot(midx,0.5,'rx','markersize',24);
        continue;
    end
    
    for vc=1:numvar
        POI = POIs{vc};
        point = POI(id);
        if ~isnan(point)
            absPoint = point - startFrame+1;
            plot(absPoint,sigstore(id,absPoint+startFrame-1), ...
                'go','markerfacecolor','g','markersize',6);
        end
    end
    hold off
end
end

