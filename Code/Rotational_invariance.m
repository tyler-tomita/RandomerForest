%TF+meandiff is trained on data where each of two classes is sampled from a
%multivariate guassian with no correlation and then trained again on the
%same data but randomly rotated. The error is compared over 10 trials.

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
dims = 1:d;
nvartosample = d^(2/3);
ntrials = 10;
ntrees = 1000;
R = zeros(d,d,ntrials);
cum_rf_err = NaN(ntrials,1);
cum_rfrot_err = NaN(ntrials,1);
cum_tfp_err = NaN(ntrials,1);
cum_tfprot_err = NaN(ntrials,1);
cum_tfpmd_err = NaN(ntrials,1);
cum_tfpmdrot_err = NaN(ntrials,1);

Mu1 = 1./dims;
Mu0 = -1*Mu1;
Mu = cat(1,Mu1,Mu0);
Sigma = ones(1,d);
obj = gmdistribution(Mu,Sigma);

for trial = 1:ntrials
    
    fprintf('trial %d\n',trial)
    
    [X,idx] = random(obj,n);
    Y = idx;
    Y(idx==2) = 0;
    Ystr = cellstr(num2str(Y));
    R(:,:,trial) = random_rotation(d);
    Xrot = X*R(:,:,trial);
    
    rf = rpclassificationforest2(ntrees,X,Ystr,'nvartosample',nvartosample,'RandomForest',true);
    cum_rf_err(trial) = oobpredict(rf,X,Ystr,'last');
    
    rfrot = rpclassificationforest2(ntrees,Xrot,Ystr,'nvartosample',nvartosample,'RandomForest',true);
    cum_rfrot_err(trial) = oobpredict(rfrot,Xrot,Ystr,'last');
    
    tfp = rpclassificationforest(ntrees,X,Ystr,'s',3,'nvartosample',nvartosample,'mdiff','off','sparsemethod','gaussian');
    cum_tfp_err(trial) = oobpredict(tfp,X,Ystr,'last');
    
    tfprot = rpclassificationforest(ntrees,Xrot,Ystr,'s',3,'nvartosample',nvartosample,'mdiff','off','sparsemethod','gaussian');
    cum_tfprot_err(trial) = oobpredict(tfprot,Xrot,Ystr,'last');
    
    tfpmd = rpclassificationforest2(ntrees,X,Ystr,'s',3,'nvartosample',nvartosample,'mdiff','on','sparsemethod','new');
    cum_tfpmd_err(trial) = oobpredict(tfpmd,X,Ystr,'last');
    
    tfpmdrot = rpclassificationforest2(ntrees,Xrot,Ystr,'s',3,'nvartosample',nvartosample,'mdiff','on','sparsemethod','new');
    cum_tfpmdrot_err(trial) = oobpredict(tfpmdrot,Xrot,Ystr,'last');
end

save('Invariance/Rotational_invariance.mat','R','cum_rf_err','cum_rfrot_err','cum_tfpmd_err','cum_tfpmdrot_err')

mean_rf = mean(cum_rf_err);
mean_rfrot = mean(cum_rfrot_err);
mean_tfp = mean(cum_tfp_err);
mean_tfprot = mean(cum_tfprot_err);
mean_tfpmd = mean(cum_tfpmd_err);
mean_tfpmdrot = mean(cum_tfpmdrot_err);

sem_rf = std(cum_rf_err)/sqrt(ntrials);
sem_rfrot = std(cum_rfrot_err)/sqrt(ntrials);
sem_tfp = std(cum_tfp_err)/sqrt(ntrials);
sem_tfprot = std(cum_tfprot_err)/sqrt(ntrials);
sem_tfpmd = std(cum_tfpmd_err)/sqrt(ntrials);
sem_tfpmdrot = std(cum_tfpmdrot_err)/sqrt(ntrials);

errorbar_groups([mean_rf mean_tfp mean_tfpmd;mean_rfrot mean_tfprot mean_tfpmdrot],[sem_rf sem_tfp sem_tfpmd;sem_rfrot sem_tfprot sem_tfpmdrot],'bar_names',{'RF' 'TF+' 'TF+md'});

ylabel('OOB Error')
title('Trunk')
legend('Untransformed','Rotated')
fname = sprintf('Invariance/Rotational_invariance_n%.0f_d%.0f_var%.0f_embed%.0f_ntrees%.0f_ntrials%.0f',n,d,Sigma(1),nvartosample,ntrees,ntrials);
save_fig(gcf,fname)