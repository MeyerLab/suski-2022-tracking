function mask=segmentdeflections_bwboundaries_adaptive(mask, nucr, gate)
[B,L]=bwboundaries(mask,'noholes');
%[B,L]=bwboundaries(mask);
obnum=max(L(:));
bordermask=zeros(size(mask));

if isempty(gate)
    gate = true(obnum,1);
end
% %%% Debugging
% figure(1),hold on
% %tic;
for ci=1:obnum
    orderedset=B{ci};
    orderedset=[orderedset(end:-1:1,2) orderedset(end:-1:1,1)]; % If using bwboundaries, reverse the order.
    if gate(ci)
        [bordermask,~,vIdx,bridgeEnds]=splitdeflections_4_bwboundaries_test_global(orderedset,bordermask,nucr);   
    end
% %%% Debugging
%     plot(orderedset(:,1),orderedset(:,2))
%     scatter(orderedset(vIdx,1),orderedset(vIdx,2))
%     for i = 1:length(bridgeEnds)
%         plot(bridgeEnds{i}(:,1),bridgeEnds{i}(:,2));
%     end
    
end
mask=mask & ~bordermask;
mask=~bwmorph(~mask,'diag');

%%% Debugging
% toc;
% hold off
% h1 = gca;
% 
% figure(2),hold on
% tic;
% for ci=1:obnum
%     orderedset=B{ci};
%     % If using bwboundaries, reverse the order.
%     orderedset=[orderedset(end:-1:1,2) orderedset(end:-1:1,1)];
%     
%     [bordermask,~,vIdx,bridgeEnds]=splitdeflections_4_bwboundaries_test(orderedset,bordermask,nucr);    
%     plot(orderedset(:,1),orderedset(:,2))
%     
%     scatter(orderedset(vIdx,1),orderedset(vIdx,2))
%     for i = 1:length(bridgeEnds)
%         plot(bridgeEnds{i}(:,1),bridgeEnds{i}(:,2));
%     end
% end
% toc;
% hold off
% h2 = gca;
% linkaxes([h1 h2]);

end