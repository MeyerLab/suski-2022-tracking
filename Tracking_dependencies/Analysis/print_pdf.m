function print_pdf(name)

set(gcf,'Units','Inches');
set(gcf,'renderer','Painters')
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(name,'-dpdf');

end

