function Yhats = stump_ensemble_predict(g,Xtest)

if strcmp(g.method,'rf')
    Yhats = 1/2*(((g.mu1_hat(1)-g.mu0_hat(1))*(Xtest(:,1)-g.mu_ave_hat(1))>0) + ...
        ((g.mu1_hat(2)-g.mu0_hat(2))*(Xtest(:,2)-g.mu_ave_hat(2))>0)) > 1/2;
else
    Xtest = [Xtest sum(Xtest,2) -diff(Xtest,1,2)];
    Yhats = 1/4*(((g.mu1_hat(1)-g.mu0_hat(1))*(Xtest(:,1)-g.mu_ave_hat(1))>0) + ...
        ((g.mu1_hat(2)-g.mu0_hat(2))*(Xtest(:,2)-g.mu_ave_hat(2))>0) + ...
        ((g.mu1_hat(3)-g.mu0_hat(3))*(Xtest(:,3)-g.mu_ave_hat(3))>0) + ...
        ((g.mu1_hat(4)-g.mu0_hat(4))*(Xtest(:,4)-g.mu_ave_hat(4))>0)) > 1/2;
end

end