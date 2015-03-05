clear
close
clc

[data,txt,raw] = xlsread('Qing Dataset 70 Samples Unnormalized (12-16-2014).xlsx');
Xtrain = data(:,1:46)';
Xtrain_rank = passtorank(Xtrain);
Xtest = data(:,47:end)';
coeff = pca(Xtrain_rank);
Xtrain_pca = Xtrain_rank*coeff;
Ytrain = zeros(size(Xtrain_rank,1),1);
Ytrain(27:46) = 1;

Xtest_rank = interpolate_rank(Xtrain,Xtest);
Xtest_pca = Xtest_rank*coeff;

plot3(Xtrain_pca(Ytrain==0,1),Xtrain_pca(Ytrain==0,2),Xtrain_pca(Ytrain==0,3),'bo',Xtrain_pca(Ytrain==1,1),Xtrain_pca(Ytrain==1,2),Xtrain_pca(Ytrain==1,3),'rx',Xtest_pca(:,1),Xtest_pca(:,2),Xtest_pca(:,3),'k.')
set(gcf,'Visible','On')
legend('Normal','Cancer','Unknown')
fname = '~/Documents/MATLAB/CancerAnalysis/Plots/Qing_rank_pca';
save_fig(gcf,fname)