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
dims = round(logspace(log10(2),3,5));
ntrials = 10;
ntrees = 1000;
cumberr = NaN(ntrials,length(dims));
cumf1err = NaN(ntrials,length(dims));
cumf2err = NaN(ntrials,length(dims));
cumf3err = NaN(ntrials,length(dims));
tb = NaN(ntrials,length(dims));
tf1 = NaN(ntrials,length(dims));
tf2 = NaN(ntrials,length(dims));
tf3 = NaN(ntrials,length(dims));

for trial = 1:ntrials
    fprintf('trial %d\n',trial)
    for i = 1:length(dims)
        d = dims(i);
        fprintf('dimensions = %d\n',dims(i))
        d_idx = 1:d;
        nvartosample = ceil(d^(2/3));
        mu1 = 1./d_idx;
        mu0 = -1*mu1;
        mu = cat(1,mu1,mu0);
        Sigma = ones(1,d);
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

save('Trunk_vary_d.mat','cumberr','cumf1err','cumf2err','cumf3err','tb','tf1','tf2','tf3')
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
    errorbar(dims,eval(Ynames{i}),eval(Enames{i}),lspec{i});
end
set(gca,'XScale','log')
xlabel('# ambient dimensions')
ylabel('oob error')
title('Trunk')
legend('RandomForest','TylerForest','TylerForest+','TylerForest+meandiff')
fname = sprintf('Trunk/trunk_ooberror_vs_d_n%d_var%d_embed%d_ntrees%d_ntrials%d',n,Sigma(1),nvartosample,ntrees,ntrials);
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
    errorbar(dims,eval(Ynames{i}),eval(Enames{i}),lspec{i});
end
set(gca,'XScale','log')
xlabel('# ambient dimensions')
ylabel('Training time (sec)')
title('Trunk')
legend('RandomForest','TylerForest','TylerForest+','TylerForest+meandiff')
fname = sprintf('Trunk/trunk_time_vs_d_n%d_var%d_embed%d_ntrees%d_ntrials%d',n,Sigma(1),nvartosample,ntrees,ntrials);
save_fig(gcf,fname)