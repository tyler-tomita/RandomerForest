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
dims = round(logspace(log10(2),3,10));
ntrials = 10;
ntrees = 1500;
NWorkers = 2;
cumrerferr = NaN(ntrials,length(dims));
cumrerf2err = NaN(ntrials,length(dims));
cumrerf3err = NaN(ntrials,length(dims));
Class = [0;1];

poolobj = gcp('nocreate');
if isempty(poolobj)
    parpool('local',NWorkers,'IdleTimeout',360);
end

for trial = 1:ntrials
    fprintf('trial %d\n',trial)
    for i = 1:length(dims)
        d = dims(i);
        fprintf('dimensions = %d\n',dims(i))
        d_idx = 1:d;
        nvartosample = ceil(d^(2/3));
        mu1 = 1./sqrt(d_idx);
        mu0 = -1*mu1;
        Mu = cat(1,mu0,mu1);
        Sigma = ones(1,d);
        obj = gmdistribution(Mu,Sigma);
        [X,idx] = random(obj,n);
        Y = cellstr(num2str(Class(idx)));
        
        rerf = rpclassificationforest(ntrees,X,Y,'sparsemethod','sparse','nvartosample',ceil(sqrt(d)),'NWorkers',NWorkers);
        cumrerferr(trial,i) = oobpredict(rerf,X,Y,'last');
        
        fprintf('RF complete\n')
        
        rerf2 = rpclassificationforest(ntrees,X,Y,'sparsemethod','sparse','nvartosample',ceil(d^(5/8)),'NWorkers',NWorkers);
        cumrerf2err(trial,i) = oobpredict(rerf2,X,Y,'last');
        
        fprintf('RF2 complete\n')
        
        rerf3 = rpclassificationforest(ntrees,X,Y,'sparsemethod','sparse','nvartosample',ceil(d^(3/4)),'NWorkers',NWorkers);
        cumrerf3err(trial,i) = oobpredict(rerf3,X,Y,'last');
        
        fprintf('RF3 complete\n')
    end
end

rerfsem = std(cumrerferr)/sqrt(ntrials);
rerf2sem  = std(cumrerf2err)/sqrt(ntrials);
rerf3sem  = std(cumrerf3err)/sqrt(ntrials);
cumrerferr = mean(cumrerferr);
cumrerf2err = mean(cumrerf2err);
cumrerf3err = mean(cumrerf3err);
Ynames = {'cumrerferr' 'cumrerf2err' 'cumrerf3err'};
Enames = {'rerfsem' 'rerf2sem' 'rerf3sem'};
lspec = {'-bo','-rx','-gs'};
hold on
for i = 1:length(Ynames)
    errorbar(dims,eval(Ynames{i}),eval(Enames{i}),lspec{i});
end
set(gca,'XScale','log')
xlabel('# ambient dimensions')
ylabel('oob error')
title('Trunk')
legend('mtry = d^1^/^2','mtry = d^5^/^8','mtry = d^3^/^4')

fname = 'Trunk_cross-validation_RerF';
save_fig(gcf,fname)