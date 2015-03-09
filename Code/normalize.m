function Xnorm = normalize(X)
    nrows = size(X,1);
    allzero = ~any(X);
    Xnorm(:,~allzero) = (X(:,~allzero) - repmat(min(X(:,~allzero)),nrows,1))./(repmat(max(X(:,~allzero)),nrows,1) - repmat(min(X(:,~allzero)),nrows,1));
    Xnorm(:,allzero) = 0;
end