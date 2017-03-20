function g = stump_ensemble_train(X,method)

if strcmp(method,'rf')
    g.method = 'rf';
else
    g.method = 'rerf';
    X = [X sum(X,2) -diff(X,1,2)];
end

g.mu_hat = mean(X,1);

end