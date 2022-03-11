%%
    FolderName = ['Figures\'];   % Your destination folder
if ~exist(FolderName)
    mkdir(FolderName);
end
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
  FigHandle = FigList(iFig);
  FigName   = get(FigHandle, 'Number');
  saveas(FigHandle, [FolderName, num2str(FigName) '.fig']);
   saveas(FigHandle, [FolderName,  num2str(FigName) '.epsc']);
    saveas(FigHandle, [FolderName,  num2str(FigName) '.png']);
end