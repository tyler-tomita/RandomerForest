function p = posterior_map(forest,xmin,xmax,ymin,ymax,npoints)

    [xgv,ygv] = meshgrid(linspace(xmin,xmax,npoints),linspace(ymin,ymax,npoints));
    X = xgv(:);
    Y = ygv(:);
    posteriors = rerf_classprob(forest,[X Y zeros(npoints^2,18)]);
    p = surf(xgv,ygv,reshape(posteriors(:,2),npoints,npoints));
    C = p.CData;
    p.ZData = zeros(npoints,npoints);
    p.CData = C;
    colormap(gray)
    p.FaceAlpha = 0.25;
    p.LineStyle = 'none';
    axis([xmin xmax ymin ymax])
    view(2)
    colorbar
    grid off
    
end