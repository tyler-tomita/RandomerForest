function h = heatmap(Z,XTickLabel,YTickLabel,ColorMap,IsBinned,varargin)

    imagesc(Z);
    colormap(ColorMap)
    h = gca;
    if nargin > 5
        Orientation = varargin(1);
    else
        Orientation = 'vertical';
    end
    if strcmp(varargin(1),'horizontal')
        if IsBinned
            h.XTick = 0.5:size(Z,2)+1;
            h.XTickLabelRotation = 90;
        else
            h.XTick = 1:size(Z,2);
        end
        h.YTick = 1:size(Z,1);
    else
        h.XTick = 1:size(Z,2);
        if IsBinned
            h.YTick = 0.5:size(Z,1)+1;
        else
            h.YTick = 1:size(Z,1);
        end
    end
    h.XTickLabel = XTickLabel;
%     h.XAxisLocation = 'top';
    h.YTickLabel = YTickLabel;
    h.TickDir = 'out';
    
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
    %     h.LineStyle = 'none';
    %     axis([xmin xmax ymin ymax])
    %     view(2)
    %     grid off
end