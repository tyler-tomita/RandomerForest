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

ns = logspace(2,6,5);    %# samples
d = 10;
nvartosample = ceil(d^(2/3));
d_idx = 1:d;
ntrees = 1000;
ntrials = 10;
NWorkers = 2;
Class = [0;1];

rf_err = NaN(ntrials,length(ns));
f1_err = NaN(ntrials,length(ns));
rf_rot_err = NaN(ntrials,length(ns));
f1_rot_err = NaN(ntrials,length(ns));

for trial = 1:ntrials
    
    fprintf('trial %d\n',trial)
    
    for i = 1:length(ns)
    
        n = ns(i);
        fprintf('n = %d\n',n)
        X = zeros(n,d);
        Sigma = 1/32*ones(1,d);
        Mu = zeros(n,d);
        for j = 1:n
            Mu(j,:) = binornd(1,0.5,1,d);
            X(j,:) = mvnrnd(Mu(j,:),Sigma);
        end
        
        nones = sum(Mu,2);
        Y = cellstr(num2str(mod(nones,2)));
        
        R = random_rotation(d);
        X_rot = X*R;

        rf = rpclassificationforest(ntrees,X,Y,'nvartosample',nvartosample,'mdiff','off','RandomForest',true,'NWorkers',NWorkers);
        rf_err(trial,i) = oobpredict(rf,X,Y,'last');
        clear rf

        f1 = rpclassificationforest(ntrees,X,Y,'nvartosample',nvartosample,'mdiff','off','sparsemethod','dense','s',3,'NWorkers',NWorkers);
        f1_err(trial,i) = oobpredict(f1,X,Y,'last');

        rf_rot = rpclassificationforest(ntrees,X_rot,Y,'nvartosample',nvartosample,'mdiff','off','RandomForest',true,'NWorkers',NWorkers);
        rf_rot_err(trial,i) = oobpredict(rf_rot,X_rot,Y,'last');
        clear rf_rot

        f1_rot = rpclassificationforest(ntrees,X_rot,Y,'nvartosample',nvartosample,'mdiff','off','sparsemethod','dense','s',3,'NWorkers',NWorkers);
        f1_rot_err(trial,i) = oobpredict(f1_rot,X_rot,Y,'last');
        clear f1_rot
    end
end

mean_rf_err = mean(rf_err);
mean_f1_err = mean(f1_err);
mean_rf_rot_err = mean(rf_rot_err);
mean_f1_rot_err = mean(f1_rot_err);

sem_rf = std(rf_err)/sqrt(ntrials);
sem_f1 = std(f1_err)/sqrt(ntrials);
sem_rf_rot = std(rf_rot_err)/sqrt(ntrials);
sem_f1_rot = std(f1_rot_err)/sqrt(ntrials);

Ynames = {'mean_rf_err' 'mean_rf_rot_err' 'mean_f1_err' 'mean_f1_rot_err'};
Enames = {'sem_rf' 'sem_rf_rot' 'sem_f1' 'sem_f1_rot'};
lspec = {'-bs','-rs','-gs' '-cs'};
facespec = {'b','r','g' 'c'};
figure(1)
for i = 1:length(Ynames)
    errorbar(ns,eval(Ynames{i}),eval(Enames{i}),lspec{i},'MarkerEdgeColor','k','MarkerFaceColor',facespec{i});
    hold on
end
set(gca,'XScale','log')
xlabel('n')
ylabel('Misclassification Rate')
legend('Untransformed RF','Rotated RF','Untransformed RerF','Rotated RerF')
title('Parity (d = 10)')

%filename = 'Parity_vary_n';
%save_fig(gcf,filename)
