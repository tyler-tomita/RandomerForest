function XY = bvn_level_curve(Mu,Sigma,Level,npoints)
    if length(Mu) > 2 || ~all(size(Sigma) == [2,2])
        error('Number of dimensions must be less than or equal to 2')
    else
        C = chol(Sigma);
        angle = linspace(0,2*pi,npoints)';
        xy = Level*[sin(angle) cos(angle)];
        XY = xy*C;
        XY = XY + repmat(Mu,npoints,1);
    end
end