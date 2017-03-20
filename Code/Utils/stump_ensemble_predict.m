function Yhats = stump_ensemble_predict(g,Xtest)

if strcmp(g.method,'rf')
    Yhats = 1/2*((Xtest(:,1)>g.mu_hat(1)) + (Xtest(:,2)>g.mu_hat(2))) > 1/2;
else
    Xtest = [Xtest sum(Xtest,2) -diff(Xtest,1,2)];
    Yhats = 1/4*((Xtest(:,1)>g.mu_hat(1)) + (Xtest(:,2)>g.mu_hat(2)) + ...
        (Xtest(:,3)>g.mu_hat(3)) + (Xtest(:,4)>g.mu_hat(4))) > 1/2;
end

end