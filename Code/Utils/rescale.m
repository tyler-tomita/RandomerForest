function Xscale = rescale(Xtrain,Xtest,Method)
    if strcmp(Method,'normalize')
        Mn = min(Xtrain);
        Mx = max(Xtrain);
        if isempty(Xtest)
            [n,~] = size(Xtrain);
            Xscale = (Xtrain - repmat(Mn,n,1))./(repmat(Mx,n,1) - repmat(Mn,n,1));
        else
            [n,~] = size(Xtest);
            Xscale = (Xtest - repmat(Mn,n,1))./(repmat(Mx,n,1) - repmat(Mn,n,1));
        end
        Xscale(:,~any(diff(Xscale))) = 0;
    elseif strcmp(Method,'rank')
        if isempty(Xtest)
            Xscale = tiedrank(Xtrain);
        else
            Xscale = interpolate_rank(Xtrain,Xtest);
        end
    elseif strcmp(Method,'zscore')
        Mu_hat = mean(Xtrain);
        Sigma_hat = std(Xtrain);
        if isempty(Xtest)
            [n,~] = size(Xtrain);
            Xscale = (Xtrain - repmat(Mu_hat,n,1))./repmat(Sigma_hat,n,1);
        else
            [n,~] = size(Xtest);
            Xscale = (Xtest - repmat(Mu_hat,n,1))./repmat(Sigma_hat,n,1);
        end
        Xscale(:,~any(diff(Xscale))) = 0;
    else
        error('Method must be either ''normalize'', ''rank'', or ''zscore''');
    end
end