clear
close
clc

[data,txt,raw] = xlsread('Qing Dataset 70 Samples Unnormalized (12-16-2014).xlsx');
Xtrain = data(:,1:46)';
Xtrain = passtorank(Xtrain);
Xtest = data(:,47:end)';
coeff = pca(Xtrain);
Xtrain_pca = Xtrain*coeff;
Ytrain = zeros(size(Xtrain,1),1);
Ytrain(27:46) = 1;
Y_hat = NaN(size(Xtest,1),1);
Ystr = cellstr(num2str(Ytrain));
parms.types={'DENL';'NENL'};
parms.ks=1:10;
[Lhat,Yhats] = LOL_loocv(Xtrain_pca,Ystr,parms);

plot(parms.ks,Lhat(1,:),'-b',parms.ks,Lhat(2,:),'-r')
xlabel('k')
ylabel('Lhat')
legend(parms.types{1},parms.types{2})
save_fig(gcf,'CancerAnalysis/Qing_rank_LOL_loocv')