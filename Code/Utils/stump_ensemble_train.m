function g = stump_ensemble_train(X,Y,method)

if strcmp(method,'rf')
    g.method = 'rf';
else
    g.method = 'rerf';
    X = [X sum(X,2) -diff(X,1,2)];
end

g.mu0_hat = mean(X(Y==0,:),1);
g.mu1_hat = mean(X(Y==1,:),1);
g.mu_ave_hat = mean([g.mu0_hat;g.mu1_hat],1);

end