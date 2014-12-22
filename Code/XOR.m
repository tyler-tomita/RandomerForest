clear
close all
clc

task.D = 1000;
task.n=100;
ntrees = 500;
ntrials = 10;

Sigmas = logspace(-1,1,5);

cumberr = NaN(ntrials,length(Sigmas));
cumf1err = NaN(ntrials,length(Sigmas));
cumf2err = NaN(ntrials,length(Sigmas));
cumf3err = NaN(ntrials,length(Sigmas));

embeddims = 4*sqrt(task.D);
for trial = 1:ntrials
    trial
    for i = 1:length(Sigmas)

        mu0 = zeros(task.D,1);
        mu1 = repmat([1;0],task.D/2,1);
        Sigma = Sigmas(i)*ones(1,task.D);
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
        cumberr(trial,i) = berr(end);
        f1 = rpclassificationforest2(ntrees,X,Ystr,'s',3,'nvartosample',embeddims,'mdiff','off','sparsemethod','new');
        f1err = oobpredict(f1,X,Ystr,'every');
        cumf1err(trial,i) = f1err(end);
        f2 = rpclassificationforest2(ntrees,X,Ystr,'s',3,'nvartosample',embeddims,'mdiff','off','sparsemethod','old');
        f2err = oobpredict(f2,X,Ystr,'every');
        cumf2err(trial,i) = f2err(end);
        f3 = rpclassificationforest2(ntrees,X,Ystr,'s',3,'nvartosample',embeddims,'mdiff','on','sparsemethod','new');
        f3err = oobpredict(f3,X,Ystr,'every');
        cumf3err(trial,i) = f3err(end);
        %figure(i)
        %plot(1:ntrees,berr,'b',1:ntrees,f1err,'r',1:ntrees,f2err,'g',1:ntrees,f3err,'k')
        %xlabel('# of Trees')
        %ylabel('OOB Error')
        %legend('RandomForest','TylerForest+','TylerForest','TylerForest+meandiff')
    end
end
bsem = std(cumberr)/sqrt(ntrials);
f1sem = std(cumf1err)/sqrt(ntrials);
f2sem  = std(cumf2err)/sqrt(ntrials);
f3sem = std(cumf3err)/sqrt(ntrials);
cumberr = mean(cumberr);
cumf1err = mean(cumf1err);
cumf2err = mean(cumf2err);
cumf3err = mean(cumf3err);
Ynames = {'cumberr' 'cumf1err' 'cumf2err' 'cumf3err'};
Enames = {'bsem' 'f1sem' 'f2sem' 'f3sem'};
lspec = {'-bo','-rx','-gd','-ks'};
hold on
for i = 1:length(Ynames)
    errorbar(Sigmas,eval(Ynames{i}),eval(Enames{i}),lspec{i});
end
set(gca,'XScale','log')
xlabel('variance')
ylabel('oob error')
legend('RandomForest','TylerForest+','TylerForest','TylerForest+meandiff')