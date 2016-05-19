function p = posterior_map(X,Y,posteriors)

    npoints = sqrt(length(posteriors));
    xmin = min(X);
    xmax = max(X);
    ymin = min(Y);
    ymax = max(Y);
    p = surf(reshape(X,npoints,npoints),reshape(Y,npoints,npoints),reshape(posteriors(:,2),npoints,npoints));
    C = p.CData;
    p.ZData = zeros(npoints,npoints);
    p.CData = C;
    colormap(gray)
    p.LineStyle = 'none';
    axis([xmin xmax ymin ymax])
    view(2)
    grid off
    
end