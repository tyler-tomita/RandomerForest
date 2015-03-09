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

ntrees = 2000;
ntrials = 10;

%read in unnormalized data
[data,txt,raw] = xlsread('Qing Dataset 70 Samples Unnormalized (12-16-2014).xlsx');
data = data(:,cat(2,1:5,7:46));   %exclude NORM_2361916 from training set
X = data';
Xrank = passtorank(X);
mean_rank = mean(Xrank,2);
Y = zeros(size(X,1),1);
Y(1:25) = 0;
Y(26:45) = 1;
Ystr = cellstr(num2str(Y));

%read in normalized data
data_norm = xlsread('Qing Dataset (Complete).xlsx');
data_norm = data_norm(:,cat(2,5:9,11:50));  %exclude NORM_2361916 from training set
X_norm = data_norm';
Xrank_norm = passtorank(X_norm);
mean_rank_norm = mean(Xrank_norm,2);

d = size(X,2);
nvartosample = ceil(d^(2/3));

err_rf = NaN(ntrials,1);
err_tfp = NaN(ntrials,1);
err_tfpmd = NaN(ntrials,1);
err_tfpmd_rank = NaN(ntrials,1);
err_rf_norm = NaN(ntrials,1);
err_tfp_norm = NaN(ntrials,1);
err_tfpmd_norm = NaN(ntrials,1);
err_tfpmd_rank_norm = NaN(ntrials,1);

for trial = 1:ntrials
    fprintf('trial %d\n',trial)
    
    rf = rpclassificationforest2(ntrees,X,Ystr,'nvartosample',nvartosample,'mdiff','off','RandomForest',true);
    err_rf(trial) = oobpredict(rf,X,Ystr,'last');
    clear rf
    
    tfp = rpclassificationforest2(ntrees,X,Ystr,'nvartosample',nvartosample,'mdiff','off','sparsemethod','new');
    err_tfp(trial) = oobpredict(tfp,X,Ystr,'last');
    clear tfp
    
    tfpmd = rpclassificationforest2(ntrees,X,Ystr,'nvartosample',nvartosample,'mdiff','on','sparsemethod','new');
    err_tfpmd(trial) = oobpredict(tfpmd,X,Ystr,'last');
    clear tfpmd
    
    tfpmd_rank = rpclassificationforest2(ntrees,Xrank,Ystr,'nvartosample',nvartosample,'mdiff','on','sparsemethod','new');
    err_tfpmd_rank(trial) = oobpredict(tfpmd_rank,Xrank,Ystr,'last');
    clear tfpmd_rank
    
    rf_norm = rpclassificationforest2(ntrees,X_norm,Ystr,'nvartosample',nvartosample,'mdiff','off','RandomForest',true);
    err_rf_norm(trial) = oobpredict(rf_norm,X_norm,Ystr,'last');
    clear rf_norm
    
    tfp_norm = rpclassificationforest2(ntrees,X_norm,Ystr,'nvartosample',nvartosample,'mdiff','off','sparsemethod','new');
    err_tfp_norm(trial) = oobpredict(tfp_norm,X_norm,Ystr,'last');
    clear tfp_norm
    
    tfpmd_norm = rpclassificationforest2(ntrees,X_norm,Ystr,'nvartosample',nvartosample,'mdiff','on','sparsemethod','new');
    err_tfpmd_norm(trial) = oobpredict(tfpmd_norm,X_norm,Ystr,'last');
    clear tfpmd_norm
    
    tfpmd_rank_norm = rpclassificationforest2(ntrees,Xrank_norm,Ystr,'nvartosample',nvartosample,'mdiff','on','sparsemethod','new');
    err_tfpmd_rank_norm(trial) = oobpredict(tfpmd_rank_norm,Xrank_norm,Ystr,'last');
    clear tfpmd_rank_norm
end

mean_rf = mean(err_rf);
mean_tfp = mean(err_tfp);
mean_tfpmd = mean(err_tfpmd);
mean_tfpmd_rank = mean(err_tfpmd_rank);
mean_rf_norm = mean(err_rf_norm);
mean_tfp_norm = mean(err_tfp_norm);
mean_tfpmd_norm = mean(err_tfpmd_norm);
mean_tfpmd_rank_norm = mean(err_tfpmd_rank_norm);

sem_rf = std(err_rf)/sqrt(ntrials);
sem_tfp = std(err_tfp)/sqrt(ntrials);
sem_tfpmd = std(err_tfpmd)/sqrt(ntrials);
sem_tfpmd_rank = std(err_tfpmd_rank)/sqrt(ntrials);
sem_rf_norm = std(err_rf_norm)/sqrt(ntrials);
sem_tfp_norm = std(err_tfp_norm)/sqrt(ntrials);
sem_tfpmd_norm = std(err_tfpmd_norm)/sqrt(ntrials);
sem_tfpmd_rank_norm = std(err_tfpmd_rank_norm)/sqrt(ntrials);

mr_classifier = fitcdiscr(mean_rank,Ystr);
resuberr_mr = resubLoss(mr_classifier);

svm = svmtrain(X,Ystr);
Y_hat = svmclassify(svm,X);
resuberr_svm = sum(~strcmp(Y_hat,Ystr))/length(Y_hat);

mr_norm_classifier = fitcdiscr(mean_rank_norm,Ystr);
resuberr_mr_norm = resubLoss(mr_norm_classifier);

svm_norm = svmtrain(X_norm,Ystr);
Y_hat_norm = svmclassify(svm,X_norm);
resuberr_svm_norm = sum(~strcmp(Y_hat_norm,Ystr))/length(Y_hat_norm);

save('CancerAnalysis/Scripts/Qing_Comprehensive.mat','err_rf','err_tfp','err_tfpmd','err_tfpmd_rank','err_rf_norm','err_tfp_norm','err_tfpmd_norm',...
    'err_tfpmd_rank_norm','resuberr_mr','resuberr_svm','resuberr_mr_norm','resuberr_svm_norm')

errorbar_groups([mean_rf mean_rf_norm;mean_tfp mean_tfp_norm;mean_tfpmd mean_tfpmd_norm;mean_tfpmd_rank mean_tfpmd_rank_norm;...
    resuberr_mr resuberr_mr_norm;resuberr_svm resuberr_svm_norm],[sem_rf sem_rf_norm;sem_tfp sem_tfp_norm;sem_tfpmd sem_tfpmd_norm;...
    sem_tfpmd_rank sem_tfpmd_rank_norm;0 0;0 0])

legend('RF','TF+','TF+md','TF+md_rank','Mean rank','Linear SVM')
ylabel('O.O.B. Error or Resubstution Error')

save_fig(gcf,'CancerAnalysis/Qing_error_barplot')