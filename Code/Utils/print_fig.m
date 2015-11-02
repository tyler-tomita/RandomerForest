function print_fig(h,F)
% h: figure handle
% fname: name of fig file & script
% wh: set paper width and height
if ~isfield(F,'renderer'), F.renderer='painters'; end
if isfield(F,'wh')
    F.PaperSize=F.wh;
    F.PaperPosition=[0 0 F.wh];
end
if ~isfield(F,'PaperSize'), F.PaperSize=[2 2]; end
if ~isfield(F,'PaperPosition'), F.PaperPosition=[0 0 F.PaperSize]; end
if ~isfield(F,'fname'), F.fname='temp_fig'; end
set(h,'PaperSize',F.PaperSize,'PaperPosition',F.PaperPosition,'color','w');
set(h, 'InvertHardCopy', 'off');
set(h,'renderer',F.renderer)
print(h,F.fname,'-dpdf')
print(h,F.fname,'-dpng','-r300')
saveas(h,F.fname,'fig')