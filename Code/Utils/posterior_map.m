function p = posterior_map(X,Y,posteriors,binarize,ColorMap)

    npoints = sqrt(length(posteriors));
    xmin = min(X);
    xmax = max(X);
    ymin = min(Y);
    ymax = max(Y);
    if ~binarize
        p = surf(reshape(X,npoints,npoints),reshape(Y,npoints,npoints),reshape(posteriors(:,2),npoints,npoints));
        C = p.CData;
        p.ZData = zeros(npoints,npoints);
        p.CData = C;
        colormap(ColorMap)
        p.LineStyle = 'none';
        axis([xmin xmax ymin ymax])
        view(2)
        grid off
    else
        posteriors = posteriors(:,2);
        posteriors(posteriors>=0.5) = 1;
        posteriors(posteriors<0.5) = 0;
        p = surf(reshape(X,npoints,npoints),reshape(Y,npoints,npoints),reshape(posteriors,npoints,npoints));
        C = p.CData;
        p.ZData = zeros(npoints,npoints);
        p.CData = C;
        colormap(ColorMap)
        p.LineStyle = 'none';
        axis([xmin xmax ymin ymax])
        view(2)
        grid off
    end
end