function mask=segmentdeflections_gated(mask,nucr,B, gate)
debug = 0;
obnum= length(B);
bordermask=zeros(size(mask));

if isempty(gate)
    gate = true(obnum,1);
end
%%% Debugging
if debug == 1
    figure(1),hold on
end
for ci=1:obnum
    orderedset=B{ci};
    orderedset=[orderedset(end:-1:1,2) orderedset(end:-1:1,1)]; % If using bwboundaries, reverse the order.
    
    if gate(ci)
        [bordermask,~,vIdx,bridgeEnds]=splitdeflections_4_bwboundaries_test_global(orderedset,bordermask,nucr);
        %%% Debugging
        if debug == 1            
            plot(orderedset(:,1),orderedset(:,2),'r')
            scatter(orderedset(vIdx,1),orderedset(vIdx,2))
            for i = 1:length(bridgeEnds)
                plot(bridgeEnds{i}(:,1),bridgeEnds{i}(:,2),'r');
            end
        end
    else
        %%% Debugging
        if debug == 1           
            plot(orderedset(:,1),orderedset(:,2),'k')
        end
    end
end
mask=mask & ~bordermask;
mask=~bwmorph(~mask,'diag');


end