close all
clear
clc

fpath = mfilename('fullpath');
findex = strfind(fpath,'/');
rootDir=fpath(1:findex(end-1));
p = genpath(rootDir);
gits=strfind(p,'.git');
colons=strfind(p,':');
for i=0:length(gits)-1
endGit=find(colons>gits(end-i),1);
p(colons(endGit-1):colons(endGit)-1)=[];
end
addpath(p);

n = 100;
d = 1000;
ntrials = 10;
ntrees = 1000;
dims = 1:d;
nvartosample = ceil(8*sqrt(d));
mu1 = 1./dims;
mu0 = -1*mu1;
mu = cat(1,mu1,mu0);
Sigmas = logspace(-1,1,5);
cumberr = NaN(ntrials,length(Sigmas));
cumf1err = NaN(ntrials,length(Sigmas));
cumf2err = NaN(ntrials,length(Sigmas));
cumf3err = NaN(ntrials,length(Sigmas));
tb = NaN(ntrials,length(Sigmas));
tf1 = NaN(ntrials,length(Sigmas));
tf2 = NaN(ntrials,length(Sigmas));
tf3 = NaN(ntrials,length(Sigmas));

for trial = 1:ntrials
    fprintf('trial %d\n',trial)
    for i = 1:length(Sigmas)
        fprintf('variance = %d\n',Sigmas(i))
        Sigma = Sigmas(i)*ones(1,d);
        obj = gmdistribution(mu,Sigma);
        [X,idx] = random(obj,n);
        Y = idx;
        Y(idx==2) = 0;
        Ystr = cellstr(num2str(Y));
        
        tic
        b = rpclassificationforest2(ntrees,X,Ystr,'nvartosample',nvartosample,'mdiff','off','RandomForest',true);
        tb(trial,i) = toc;
        berr = oobpredict(b,X,Ystr,'every');
        cumberr(trial,i) = berr(end);
        
        fprintf('Random Forest complete\n')
        
        tic
        f1 = rpclassificationforest2(ntrees,X,Ystr,'s',3,'nvartosample',nvartosample,'mdiff','off','sparsemethod','old');
        tf1(trial,i) = toc;
        f1err = oobpredict(f1,X,Ystr,'every');
        cumf1err(trial,i) = f1err(end);
        
        fprintf('TylerForest complete\n')
        
        tic
        f2 = rpclassificationforest2(ntrees,X,Ystr,'s',3,'nvartosample',nvartosample,'mdiff','off','sparsemethod','new');
        tf2(trial,i) = toc;
        f2err = oobpredict(f2,X,Ystr,'every');
        cumf2err(trial,i) = f2err(end);
        
        fprintf('TylerForest+ complete\n')
        
        tic
        f3 = rpclassificationforest2(ntrees,X,Ystr,'s',3,'nvartosample',nvartosample,'mdiff','on','sparsemethod','new');
        tf3(trial,i) = toc;
        f3err = oobpredict(f3,X,Ystr,'every');
        cumf3err(trial,i) = f3err(end);
        
        fprintf('TylerForest+meandiff complete\n')
    end
end

save('Trunk_vary_sigma.mat','cumberr','cumf1err','cumf2err','cumf3err','tb','tf1','tf2','tf3')
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
title('Trunk')
legend('RandomForest','TylerForest','TylerForest+','TylerForest+meandiff')
fname = sprintf('Trunk/trunk_ooberror_vs_variance_n%d_d%d_embed%d_ntrees%d_ntrials%d',n,d,nvartosample,ntrees,ntrials);
save_fig(gcf,fname)

figure(2)
bsem = std(tb)/sqrt(ntrials);
f1sem = std(tf1)/sqrt(ntrials);
f2sem  = std(tf2)/sqrt(ntrials);
f3sem = std(tf3)/sqrt(ntrials);
tb = mean(tb);
tf1 = mean(tf1);
tf2 = mean(tf2);
tf3 = mean(tf3);
Ynames = {'tb' 'tf1' 'tf2' 'tf3'};
Enames = {'bsem' 'f1sem' 'f2sem' 'f3sem'};
lspec = {'-bo','-rx','-gd','-ks'};
hold on
for i = 1:length(Ynames)
    errorbar(Sigmas,eval(Ynames{i}),eval(Enames{i}),lspec{i});
end
set(gca,'XScale','log')
xlabel('Variance')
ylabel('Training time (sec)')
title('Trunk')
legend('RandomForest','TylerForest','TylerForest+','TylerForest+meandiff')
fname = sprintf('Trunk/trunk_time_vs_variance_n%d_d%d_embed%d_ntrees%d_ntrials%d',n,d,nvartosample,ntrees,ntrials);
save_fig(gcf,fname)