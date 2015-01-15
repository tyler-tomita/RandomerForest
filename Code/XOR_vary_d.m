clear
close all
clc

dims = round(logspace(log10(2),3,10));
n=100;
ntrees = 1000;
ntrials = 10;

cumberr = NaN(ntrials,length(dims));
cumf1err = NaN(ntrials,length(dims));
cumf2err = NaN(ntrials,length(dims));
cumf3err = NaN(ntrials,length(dims));

hold on
for trial = 1:ntrials
    trial
    for i = 1:length(dims)
        
        d = dims(i);
        d = d - mod(d,2);
        embeddims = ceil(d^(2/3));

        mu1 = transpose(-1*ones(d,1));
        mu2 = -1*mu1;
        mu3 = repmat([-1 1],1,d/2);
        mu4 = repmat([1 -1],1,d/2);
        sigma = 8*ones(1,d);
        X = cat(1,mvnrnd(mu1,sigma,n/4),mvnrnd(mu2,sigma,n/4),mvnrnd(mu3,sigma,n/4),mvnrnd(mu4,sigma,n/4));
        Y = cat(1,zeros(n/2,1),ones(n/2,1));
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
    end
end
save('xor_vary_d.mat','cumberr','cumf1err','cumf2err','cumf3err')
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
facespec = {'b','r','g','k'};
set(gcf,'visible','off')
h = gcf;
for i = 1:length(Ynames)
    errorbar(dims,eval(Ynames{i}),eval(Enames{i}),lspec{i});
end
set(gca,'XScale','log')
xlabel('# ambient dimensions')
ylabel('oob error')
title('XOR')
legend('RandomForest','TylerForest+','TylerForest','TylerForest+meandiff')
F.fname = 'XOR/xor_ooberror_vs_d'
save_fig(h,F)