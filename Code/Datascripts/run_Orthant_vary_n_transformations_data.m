close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

rng(1);

ps = [2,4,6];
ns = {[20,200,400], [80,400,4000], [400,2000,4000]};
ntrials = 10;
OutlierProportion = 0.2;

for j = 1:length(ps)
    p = ps(j);
    fprintf('p = %d\n',p)
    Xtest = dlmread(sprintf('~/R/Data/Orthant/dat/Raw/Test/Orthant_test_p%d.dat',p));
    Ytest = Xtest(:,end);
    Xtest(:,end) = [];
    ntest = size(Xtest,1);
    for i = 1:length(ns{j})
        ntrain = ns{j}(i);
        fprintf('n = %d\n',ntrain)        
        for trial = 1:ntrials
            R = random_rotation(p);
            S = 10.^random_scaling(ntrain+ntest,p,-5,5);
            Xtrain = dlmread(sprintf('~/R/Data/Orthant/dat/Raw/Train/Orthant_train_p%d_n%d_trial%d.dat',p,ntrain,trial));
            Ytrain = Xtrain(:,end);
            Xtrain(:,end) = [];
            
            dlmwrite(sprintf('~/R/Data/Orthant/dat/Rotated/Train/Orthant_rotated_train_p%d_n%d_trial%d.dat',p,ntrain,trial),...
                [Xtrain*R,Ytrain],'delimiter','\t','precision','%0.15f');
            dlmwrite(sprintf('~/R/Data/Orthant/dat/Rotated/Test/Orthant_rotated_test_p%d_n%d_trial%d.dat',ntest,p,trial),...
                [Xtest*R,Ytest],'delimiter','\t','precision','%0.15f');
            
            dlmwrite(sprintf('~/R/Data/Orthant/dat/Scaled/Orthant_scaled_train_p%d_n%d_trial%d.dat',p,ntrain,trial),...
                [Xtrain.*S(1:ntrain,:),Ytrain],'delimiter','\t','precision','%0.15f');
            dlmwrite(sprintf('~/R/Data/Orthant/dat/Scaled/Test/Orthant_scaled_test_p%d_n%d_trial%d.dat',ntest,p,trial),...
                [Xtest.*S(ntrain+1:end,:),Ytest],'delimiter','\t','precision','%0.15f');
            
            dlmwrite(sprintf('~/R/Data/Orthant/dat/Affine/Train/Orthant_affine_train_p%d_n%d_trial%d.dat',p,ntrain,trial),...
                [(Xtrain*R).*S(1:ntrain,:),Ytrain],'delimiter','\t','precision','%0.15f');
            dlmwrite(sprintf('~/R/Data/Orthant/dat/Affine/Test/Orthant_affine_test_p%d_n%d_trial%d.dat',ntest,p,trial),...
                [(Xtest*R).*S(ntrain+1:end,:),Ytest],'delimiter','\t','precision','%0.15f');
            
            nOutlier = ceil(OutlierProportion*numel(Xtrain));
            Exponents = randsample(2:5,nOutlier,true)';
            OutlierIdx = randperm(numel(Xtrain),nOutlier)';
            Xtrain_out = Xtrain;
            Xtrain_out(OutlierIdx) = Xtrain_out(OutlierIdx).*10.^(Exponents);
            dlmwrite(sprintf('~/R/Data/Orthant/dat/Corrupted/Train/Orthant_corrupted_train_p%d_n%d_trial%d.dat',p,ntrain,trial),...
                [Xtrain_out,Ytrain],'delimiter','\t','precision','%0.15f');
            dlmwrite(sprintf('~/R/Data/Orthant/dat/Corrupted/Test/Orthant_corrupted_test_p%d_n%d_trial%d.dat',ntest,p,trial),...
                [Xtest,Ytest],'delimiter','\t','precision','%0.15f');
        end
    end   
end