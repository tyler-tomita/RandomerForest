function h = heat_map(X,Y,f,ColorMap,cmin,cmax,isColorBar)

        h = surf(X,Y,f);
        C = h.CData;
        h.ZData = zeros(size(f));
        h.CData = C;
        colormap(ColorMap)
        h.LineStyle = 'none';
        view(2)
        grid off
        caxis([cmin,cmax])
        xlim([min(X(:)) max(X(:))])
        ylim([min(Y(:)) max(Y(:))])
        if isColorBar
            colorbar;
        end
        
end