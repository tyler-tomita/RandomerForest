function h = heat_map(X,Y,f,ColorMap)

        h = surf(X,Y,f);
        C = h.CData;
        h.ZData = zeros(size(f));
        h.CData = C;
        colormap(ColorMap)
        h.LineStyle = 'none';
        view(2)
        grid off
end