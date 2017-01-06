% Generate Trunk datasets for various values of n and p

close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

rng(1);

ps = [3,10,20];                 % numbers of dimensions
ns{1} = [10,100,1000];          % numbers of train samples for ps(1)
ns{2} = [100,1000,10000];       % numbers of train samples for ps(2)
ns{3} = [1000,5000,10000];     % numbers of train samples for ps(3)
ntest = 10000;                  % numbers of test samples
ntrials = 10;                   % number of replicate experiments

% generate data
for j = 3:length(ps)
    p = ps(j);
    fprintf('p = %d\n',p)
    p_prime = min(3,p);
    for i = 1:length(ns{j})
        ntrain = ns{j}(i);
        fprintf('n = %d\n',ntrain)
        if ntrain <= 10
            for trial = 1:ntrials
                go = true;
                while go
                    Xtrain = rand(ntrain,p)*2 - 1;
                    Ytrain = mod(sum(Xtrain(:,1:p_prime)>0,2),2);
                    if mean(Ytrain==1) == 0.5
                        go = false;
                    end
                end
                dlmwrite(sprintf('~/Documents/R/Data/Sparse_parity/dat/Train/Sparse_parity_train_set_n%d_p%d_trial%d.dat',ntrain,p,trial),...
                    [Xtrain,Ytrain],'delimiter','\t','precision','%0.15f');
            end
        else
            for trial = 1:ntrials
                Xtrain = rand(ntrain,p)*2 - 1;
                Ytrain = mod(sum(Xtrain(:,1:p_prime)>0,2),2);
                dlmwrite(sprintf('~/Documents/R/Data/Sparse_parity/dat/Train/Sparse_parity_train_set_n%d_p%d_trial%d.dat',ntrain,p,trial),...
                    [Xtrain,Ytrain],'delimiter','\t','precision','%0.15f');
            end
        end
    end
    
    Xtest = rand(ntest,p)*2 - 1;
    Ytest = mod(sum(Xtest(:,1:p_prime)>0,2),2);
    dlmwrite(sprintf('~/Documents/R/Data/Sparse_parity/dat/Test/Sparse_parity_test_set_p%d.dat',p),...
    [Xtest,Ytest],'delimiter','\t','precision','%0.15f');

    ClassPosteriors = zeros(ntest,2);
    ClassPosteriors(:,2) = cellfun(@str2double,Ytest);
    ClassPosteriors(:,1) = 1 - ClassPosteriors(:,2);
    dlmwrite(sprintf('~/Documents/R/Data/Sparse_parity/dat/Test/Sparse_parity_test_set_posteriors_p%d.dat',p),...
    ClassPosteriors,'delimiter','\t','precision','%0.15f');
end