%% 2-component GMM, p = 50, plotting error vs n for RF and RerF

close all
clear
clc

rng(1)

ns = [10,50,100,500];
p =50;
Classes = [0;1];
ntrials = 20;

Mu0 = zeros(1,p);
Mu1 = zeros(1,p);
Mu1(1) = 1;
Sigma = eye(p)/4;
obj = gmdistribution([Mu0;Mu1],Sigma);


for trial = 1:ntrials
    fprintf('trial %d\n',trial)
    for i = 1:length(ns)
        n = ns(i);
        fprintf('n = %d\n',n)
        
        [X,idx] = random(obj,n);
        Y = cellstr(num2str(Classes(idx)));
        
        %Compute Bayes error
        Posteriors = gmm_class_post(X,[Mu0;Mu1],repmat(Sigma,1,1,2));
        YBayes = cellstr(num2str(double(Posteriors(:,2)>0.5)));
        Lhat.bayes(trial,i) = mean(~strcmp(YBayes,Y));

        ntrees = 500;
        mtrys.rf = ceil(p.^[0,1/4,1/2,3/4,1]);
        mtrys.rerf = ceil(p.^[0,1/4,1/2,3/4,1,1.5,2]);
        NWorkers = 2;
        oobError.rf = [];
        oobError.rerf = [];
        
        %Compute forest errors

        for m = 1:length(mtrys.rf)
            Forest = rpclassificationforest(ntrees,X,Y,'RandomForest',true,'nvartosample',mtrys.rf(m),'NWorkers',NWorkers,'Stratified',true);
            Yhats = oobpredict(Forest,X,Y);
            oobError.rf(m) = oob_error(Yhats,Y);
        end

        for m = 1:length(mtrys.rerf)
            Forest = rpclassificationforest(ntrees,X,Y,'sparsemethod','sparse','nvartosample',mtrys.rerf(m),'NWorkers',NWorkers,'Stratified',true);
            Yhats = oobpredict(Forest,X,Y);
            oobError.rerf(m) = oob_error(Yhats,Y);
        end
        Lhat.rf(trial,i) = min(oobError.rf);
        Lhat.rerf(trial,i) = min(oobError.rerf);
    end
end

sem.bayes = std(Lhat.bayes)/sqrt(ntrials);
sem.rf = std(Lhat.rf)/sqrt(ntrials);
sem.rerf = std(Lhat.rerf)/sqrt(ntrials);
mn.bayes = mean(Lhat.bayes);
mn.rf = mean(Lhat.rf);
mn.rerf = mean(Lhat.rerf);

errorbar(ns,mn.rf,sem.rf,'m','LineWidth',2)
hold on
errorbar(ns,mn.rerf,sem.rerf,'g','LineWidth',2)
errorbar(ns,mn.bayes,sem.bayes,'r--','LineWidth',2)
xlabel('n')
ylabel('oob error')
title('2-component GMM, one informative dimension, p = 2')
legend('RF','RerF','Bayes')
ax = gca;
ax.XScale = 'log';