close all
clear
clc

n = 100;
d = 1000;
ntrials = 10;
dims = 1:d;
mu1 = 1./dims;
mu0 = -1*mu1;
mu = cat(1,mu1,mu0);
sigma = ones(1,d);
obj = gmdistribution(mu,sigma);
[X,idx] = random(obj,n);
Y = idx;
Y(idx==2) = 0;
Ystr = cellstr(num2str(Y));

embeddims = [0.5 1 2 4 8]*sqrt(d);
cumberr = NaN(ntrials,length(embeddims));
cumf1err = NaN(ntrials,length(embeddims));
cumf2err = NaN(ntrials,length(embeddims));
cumf3err = NaN(ntrials,length(embeddims));

ntrees = 500;

for trial = 1:ntrials
    trial
    for i = 1:length(embeddims)
        b = TreeBagger(ntrees,X,Y,'oobpred','on','NVarToSample',sqrt(d));
        berr = oobError(b);
        cumberr(trial,i) = berr(end);
        f1 = rpclassificationforest2(ntrees,X,Ystr,'s',3,'nvartosample',sqrt(d),'mdiff','off','sparsemethod','new');
        f1err = oobpredict(f1,X,Ystr,'every');
        cumf1err(trial,i) = f1err(end);
        f2 = rpclassificationforest2(ntrees,X,Ystr,'s',3,'nvartosample',sqrt(d),'mdiff','off','sparsemethod','old');
        f2err = oobpredict(f2,X,Ystr,'every');
        cumf2err(trial,i) = f2err(end);
        f3 = rpclassificationforest2(ntrees,X,Ystr,'s',3,'nvartosample',sqrt(d),'mdiff','on','sparsemethod','new');
        f3err = oobpredict(f3,X,Ystr,'every');
        cumf3err(trial,i) = f3err(end);
        %figure(i)
        %plot(1:ntrees,berr,'b',1:ntrees,f1err,'r',1:ntrees,f2err,'g',1:ntrees,f3err,'k')
        %xlabel('# of Trees')
        %ylabel('OOB Error')
        %legend('RandomForest','TylerForest+','TylerForest','TylerForest+meandiff')
    end
end

%figure(i+1)
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
lspec = {'bo','rx','gd','ks'};
hold on
for i = 1:length(Ynames)
    errorbar(embeddims,eval(Ynames{i}),eval(Enames{i}),lspec{i});
end
set(gca,'XScale','log')
xlabel('# embedded dimensions')
ylabel('oob error')
title('Trunk')
legend('RandomForest','TylerForest+','TylerForest','TylerForest+meandiff')