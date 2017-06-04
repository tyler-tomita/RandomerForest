% Generate 2d spherical Gaussian binary classification datasets

close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

rng(1);

ps = 2;                         % numbers of dimensions
ns = [10,50,100,500,1000];      % numbers of train samples for ps(1)
ntest = 10000;                  % number of test samples
ntrials = 50;                   % number of replicate experiments
Class = [0;1];
theta = pi/4;
R = [cos(theta) sin(theta);-sin(theta) cos(theta)];

% generate data
for j = 1:length(ps)
    p = ps(j);
    fprintf('p = %d\n',p)
    p_idx = 1:p;
    mu0 = [-1 0];
    mu1 = [1 0];
    Mu = cat(1,mu0,mu1);
    Sigma = eye(p);
    obj = gmdistribution(Mu,Sigma);
    for i = 1:length(ns)
        ntrain = ns(i);
        fprintf('n = %d\n',ntrain)
        for trial = 1:ntrials
            if ntrain == 10
                go = true;
                while go
                    [Xtrain,idx] = random(obj,ntrain);
                    Ytrain = Class(idx);
                    if mean(Ytrain==1) == 0.5
                        go = false;
                    end
                end
            else
                [Xtrain,idx] = random(obj,ntrain);
                Ytrain = Class(idx);
            end
            Xtrain_rot = Xtrain*R;
            dlmwrite(sprintf('~/R/Data/Gaussian/dat/Train/Gaussian_0_deg_train_set_n%d_p%d_trial%d.dat',ntrain,p,trial),...
                [Xtrain,Ytrain],'delimiter',',','precision','%0.15f');
            dlmwrite(sprintf('~/R/Data/Gaussian/dat/Train/Gaussian_45_deg_train_set_n%d_p%d_trial%d.dat',ntrain,p,trial),...
                [Xtrain_rot,Ytrain],'delimiter',',','precision','%0.15f');
        end
    end
    [Xtest,idx] = random(obj,ntest);
    Ytest = Class(idx);
    Xtest_rot = Xtest*R;
    dlmwrite(sprintf('~/R/Data/Gaussian/dat/Test/Gaussian_0_deg_test_set_p%d.dat',p),...
    [Xtest,Ytest],'delimiter',',','precision','%0.15f');
    dlmwrite(sprintf('~/R/Data/Gaussian/dat/Test/Gaussian_45_deg_test_set_p%d.dat',p),...
    [Xtest_rot,Ytest],'delimiter',',','precision','%0.15f');
    
    ClassPosteriors = posterior(obj,Xtest);
    dlmwrite(sprintf('~/R/Data/Gaussian/dat/Test/Gaussian_test_set_posteriors_p%d.dat',p),...
    ClassPosteriors,'delimiter',',','precision','%0.15f');
end