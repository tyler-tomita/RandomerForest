function [W,lambdas] = discriminant_projection(X,Y)
%Elements of Y must be cell strings
    labels = unique(Y);
    d = size(X,2);
    MU = zeros(length(labels),size(X,2));    %class-conditional mean vectors, where row i is the mean vector for the ith class
    mu = mean(X);    %overall mean vector
    Sw = zeros(d,d,length(labels));
    Sb = zeros(d,d,length(labels));

    for k = 1:length(labels)
        MU(k,:) = mean(X(strcmp(Y,labels(k)),:));
        nk = sum(strcmp(Y,labels(k)));
        Xk = X(strcmp(Y,labels(k)),:);
        s = zeros(d,d,nk);
        for i = 1:nk
            dff = Xk(i,:) - MU(k);
            s(:,:,i) = dff'*dff;
        end
        Sw(:,:,k) = sum(s,3);
        dff = MU(k,:) - mu;
        Sb(:,:,k) = nk*(dff'*dff);
    end
    Sw = sum(Sw,3);
    Sb = sum(Sb,3);
    if issparse(X)
        Wpca = pca(full(X));
    else
        Wpca = pca(X);
    end
    %[W,lambdas] = eig((Wpca'*Sw*Wpca)\(Wpca'*Sb*Wpca));
    if rcond(Sw) >= 1e-12
        [W,lambdas] = eig(Sw\Sb);
    else
        warning('Sw is poorly conditioned. Using pseudoinverse instead')
        [W,lambdas] = eig(pinv(Sw)*Sb);
    end
    lambdas = sum(lambdas);
    W = W(:,1:length(labels)-1);
end