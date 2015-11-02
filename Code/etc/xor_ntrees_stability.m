clear
close all
clc

task.D = 1000;
task.n=100;
ntrees = 1000;
ntrials = 10;

ntrees_b = NaN(ntrials,1);
ntrees_f1 = NaN(ntrials,1);
ntrees_f2 = NaN(ntrials,1);
ntrees_f3 = NaN(ntrials,1);

winlength = 20;

embeddims = 4*sqrt(task.D);
for trial = 1:ntrials
    trial
    mu0 = zeros(task.D,1);
    mu1 = repmat([1;0],task.D/2,1);
    Sigma = 10^(-1/2)*ones(1,task.D);
    gmm = gmdistribution([mu0,mu1]',Sigma);
    [X0,Y0] = random(gmm,task.n*0.5);

    mu0 = ones(task.D,1);
    mu1 = repmat([0;1],task.D/2,1);
    gmm = gmdistribution([mu0,mu1]',Sigma);
    [X1,Y1] = random(gmm,task.n*0.5);
    X = [X0;X1];
    Y = [Y0;Y1];
    Ystr = cellstr(num2str(Y));

    b = TreeBagger(ntrees,X,Y,'oobpred','on','NVarToSample',embeddims);
    berr = oobError(b);
    j = 1;
    while var(berr(j:j+winlength-1)) ~= 0 && j+winlength-1 < ntrees
        j = j + 1;
    end
    ntrees_b(trial) = j;
    f1 = rpclassificationforest2(ntrees,X,Ystr,'s',3,'nvartosample',embeddims,'mdiff','off','sparsemethod','new');
    f1err = oobpredict(f1,X,Ystr,'every');
    j = 1;
    while var(f1err(j:j+winlength-1)) ~= 0 && j+winlength-1 < ntrees
        j = j + 1;
    end
    ntrees_f1(trial) = j;
    f2 = rpclassificationforest2(ntrees,X,Ystr,'s',3,'nvartosample',embeddims,'mdiff','off','sparsemethod','old');
    f2err = oobpredict(f2,X,Ystr,'every');
    j = 1;
    while var(f2err(j:j+winlength-1)) ~= 0 && j+winlength-1 < ntrees
        j = j + 1;
    end
    ntrees_f2(trial) = j;
    f3 = rpclassificationforest2(ntrees,X,Ystr,'s',3,'nvartosample',embeddims,'mdiff','on','sparsemethod','new');
    f3err = oobpredict(f3,X,Ystr,'every');
    j = 1;
    while var(f3err(j:j+winlength-1)) ~= 0 && j+winlength-1 < ntrees
        j = j + 1;
    end
    ntrees_f3(trial) = j;
    figure(trial)
    plot(1:ntrees,berr,'b',1:ntrees,f1err,'r',1:ntrees,f2err,'g',1:ntrees,f3err,'k')
    xlabel('# of Trees')
    ylabel('OOB Error')
    legend('RandomForest','TylerForest+','TylerForest','TylerForest+meandiff')
end
bmean = mean(ntrees_b);
f1mean = mean(ntrees_f1);
f2mean = mean(ntrees_f2);
f3mean = mean(ntrees_f3);

figure(trial+1)
bar([bmean f1mean f2mean f3mean;0 0 0 0])
ylabel('Number of trees to stabilize')
xlim([0.5 1.5]); % trick that hides second category
legend('RandomForest','TylerForest+','TylerForest','TylerForest+meandiff')
set(gca, 'XTickLabel', ''); % hide the '1' on the X axis

for i = 1:10
    h = figure(i);
    fname = sprintf('XOR/xor_ooberror_vs_ntrees_n100_d1000_trial%d',i);
    save_fig(h,fname)
end

h = figure(11);
fname = 'XOR/xor_#trees_to_stabilize';
save_fig(h,fname)
