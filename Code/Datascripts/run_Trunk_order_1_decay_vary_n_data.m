% Generate Trunk binary classification datasets for various values of n and p

close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

rng(1);

ps = [10,100,1000];             % numbers of dimensions
ns{1} = [10,100,1000,10000];    % numbers of train samples for ps(1)
ns{2} = [10,100,1000,10000];    % numbers of train samples for ps(2)
ns{3} = [10,100,1000,10000];    % numbers of train samples for ps(3)
ntest = 10000;                  % number of test samples
ntrials = 10;                   % number of replicate experiments
Class = [0;1];

% generate data
for j = 1:length(ps)
    p = ps(j);
    fprintf('p = %d\n',p)
    p_idx = 1:p;
    mu1 = 1./p_idx;
    mu0 = -1*mu1;
    Mu = cat(1,mu0,mu1);
    Sigma = eye(p);
    obj = gmdistribution(Mu,Sigma);
    for i = 1:length(ns{j})
        ntrain = ns{j}(i);
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
            dlmwrite(sprintf('~/Documents/R/Data/Trunk/dat/Train/Trunk_order_1_decay_train_set_n%d_p%d_trial%d.dat',ntrain,p,trial),...
                [Xtrain,Ytrain],'delimiter','\t','precision','%0.15f');
        end
    end
    [Xtest,idx] = random(obj,ntest);
    Ytest = Class(idx);
    dlmwrite(sprintf('~/Documents/R/Data/Trunk/dat/Test/Trunk_order_1_decay_test_set_p%d.dat',p),...
    [Xtest,Ytest],'delimiter','\t','precision','%0.15f');
    
    ClassPosteriors = posterior(obj,Xtest);
    dlmwrite(sprintf('~/Documents/R/Data/Trunk/dat/Test/Trunk_order_1_decay_test_set_posteriors_p%d.dat',p),...
    ClassPosteriors,'delimiter','\t','precision','%0.15f');
end