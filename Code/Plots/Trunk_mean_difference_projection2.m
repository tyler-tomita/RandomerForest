%Trunk mean difference projection

clear
close all
clc
n = 100;
d = 1000;
d_idx = 1:d;
mu1 = 1./d_idx;
mu0 = -1*mu1;
Mu = cat(1,mu0,mu1);
Mu_diff = mu1' - mu0';
Sigma = 1*speye(d);
Class = [0;1];
obj = gmdistribution(Mu,Sigma);
[Xtrain,idx] = random(obj,n);
Ytrain = Class(idx);
md = transpose(mean(Xtrain(Ytrain==1,:))) - transpose(mean(Xtrain(Ytrain==0,:)));
[Xtest,idx] = random(obj,n);
Ytest = Class(idx);
Xtrain_proj_hat = Xtrain*md;
Xtrain_proj = Xtrain*Mu_diff;
Xtest_proj_hat = Xtest*md;
subplot(1,3,1)
plot(Xtrain_proj_hat(Ytrain==0),zeros(sum(Ytrain==0),1),'bo',Xtrain_proj_hat(Ytrain==1),zeros(sum(Ytrain==1,1)),'rx')
xlabel('Projection onto mdhat')
title('Xtrain')
ax = gca;
ax.YTick = [];
subplot(1,3,2)
plot(Xtrain_proj(Ytrain==0),zeros(sum(Ytrain==0),1),'bo',Xtrain_proj(Ytrain==1),zeros(sum(Ytrain==1,1)),'rx')
xlabel('Projection onto md_0')
title('Xtrain')
ax = gca;
ax.YTick = [];
subplot(1,3,3)
plot(Xtest_proj_hat(Ytest==0),zeros(sum(Ytest==0),1),'bo',Xtest_proj_hat(Ytest==1),zeros(sum(Ytest==1,1)),'rx')
xlabel('Projection onto mdhat')
title('Xtest')
ax = gca;
ax.YTick = [];

fname = '~/LOVEFest/Figures/Trunk_mean_difference_projection_v2';
save_fig(gcf,fname)