close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

rng(1);

ps = [10,100,1000];
ns = {[10,100,1000,10000], [10,100,1000,10000], [10,100,1000,10000]};
ntrials = 10;
OutlierProportion = 0.2;

for j = 1:length(ps)
    p = ps(j);
    Xtest = dlmread(sprintf('~/R/Data/Trunk/dat/Raw/Test/Trunk_test_set_p%d.dat',p));
    Ytest = Xtest(:,end);
    Xtest(:,end) = [];
    ntest = size(Xtest,1);
    for i = 1:length(ns{j})
        ntrain = ns{j}(i);
        for trial = 1:ntrials
            R = random_rotation(p);
            S = 10.^random_scaling(ntrain+ntest,p,-5,5);
            Xtrain = dlmread(sprintf('~/R/Data/Trunk/dat/Raw/Train/Trunk_train_set_n%d_p%d_trial%d.dat',ntrain,p,trial));
            Ytrain = Xtrain(:,end);
            Xtrain(:,end) = [];
            dlmwrite(sprintf('~/R/Data/Trunk/dat/Rotated/Train/Trunk_rotated_train_set_n%d_p%d_trial%d.dat',ntrain,p,trial),...
                [Xtrain*R,Ytrain],'delimiter','\t','precision','%0.15f');
            dlmwrite(sprintf('~/R/Data/Trunk/dat/Rotated/Test/Trunk_rotated_test_set_n%d_p%d_trial%d.dat',ntest,p,trial),...
                [Xtest*R,Ytest],'delimiter','\t','precision','%0.15f');
            
            dlmwrite(sprintf('~/R/Data/Trunk/dat/Scaled/Trunk_scaled_train_set_n%d_p%d_trial%d.dat',ntrain,p,trial),...
                [Xtrain.*S(1:ntrain,:),Ytrain],'delimiter','\t','precision','%0.15f');
            dlmwrite(sprintf('~/R/Data/Trunk/dat/Scaled/Test/Trunk_scaled_test_set_n%d_p%d_trial%d.dat',ntest,p,trial),...
                [Xtest.*S(ntrain+1:end,:),Ytest],'delimiter','\t','precision','%0.15f');
            
            dlmwrite(sprintf('~/R/Data/Trunk/dat/Affine/Train/Trunk_affine_train_set_n%d_p%d_trial%d.dat',ntrain,p,trial),...
                [(Xtrain*R).*S(1:ntrain,:),Ytrain],'delimiter','\t','precision','%0.15f');
            dlmwrite(sprintf('~/R/Data/Trunk/dat/Affine/Test/Trunk_affine_test_set_n%d_p%d_trial%d.dat',ntest,p,trial),...
                [(Xtest*R).*S(ntrain+1:end,:),Ytest],'delimiter','\t','precision','%0.15f');
            
            nOutlier = ceil(OutlierProportion*numel(Xtrain));
            Exponents = randsample(2:5,nOutlier,true)';
            OutlierIdx = randperm(numel(Xtrain),nOutlier)';
            Xtrain_out = Xtrain;
            Xtrain_out(OutlierIdx) = Xtrain_out(OutlierIdx).*10.^(Exponents);
            dlmwrite(sprintf('~/R/Data/Trunk/dat/Corrupted/Train/Trunk_corrupted_train_set_n%d_p%d_trial%d.dat',ntrain,p,trial),...
                [Xtrain_out,Ytrain],'delimiter','\t','precision','%0.15f');
            dlmwrite(sprintf('~/R/Data/Trunk/dat/Corrupted/Test/Trunk_corrupted_test_set_n%d_p%d_trial%d.dat',ntest,p,trial),...
                [Xtest,Ytest],'delimiter','\t','precision','%0.15f');
        end
    end   
end