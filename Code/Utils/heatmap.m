function h = heatmap(Z,XTickLabel,YTickLabel,ColorMap,IsBinned)
    imagesc(Z);
    
%     ColorMap = custom_colors(11,'purple2green');
%     [nr,nc] = size(Z);
%     xmin = min(X(:));
%     xmax = max(X(:));
%     ymin = min(Y(:));
%     ymax = max(Y(:));
%     h = surf(X,Y,Z);
%     C = h.CData;
%     h.ZData = zeros(nr,nc);
%     h.CData = C;
    colormap(ColorMap)
%     h.LineStyle = 'none';
%     axis([xmin xmax ymin ymax])
%     view(2)
%     grid off
    h = gca;
    h.XTick = 1:size(Z,2);
    if IsBinned
        h.YTick = 0.5:size(Z,1)+1;
    else
        h.YTick = 1:size(Z,1);
    end
    h.XTickLabel = XTickLabel;
    h.YTickLabel = YTickLabel;
    h.TickDir = 'out';
end