function Xnorm = normalize(X)
    nrows = size(X,1);
    has_variance = range(X) ~= 0;
    Xnorm(:,has_variance) = (X(:,has_variance) - repmat(min(X(:,has_variance)),nrows,1))./(repmat(max(X(:,has_variance)),nrows,1) - repmat(min(X(:,has_variance)),nrows,1));
    Xnorm(:,~has_variance) = 0;
end
