% Generate Trunk datasets for various values of n and p

close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

rng(1);

ps = [3,5,10];
ns{1} = [10,50,100,500,1000];
ns{2} = [50,100,500,1000,5000];
ns{3} = [100,500,1000,5000,10000];
ntest = 10000;
ntrials = 10;

Xtrain = cell(length(ns),length(ps));
Ytrain = cell(length(ns),length(ps));
Xtest = cell(1,length(ps));
Ytest = cell(1,length(ps));
% sample training data
for j = 1:length(ps)
    p = ps(j);
    fprintf('p = %d\n',p)
    p_prime = min(3,p);
    for i = 1:length(ns{j})
        ntrain = ns{j}(i);
        fprintf('n = %d\n',ntrain)
        Xtrain{i,j} = rand(ntrain,p,ntrials)*2 - 1;
        Ytrain{i,j} = cell(ntrain,ntrials);
        for trial = 1:ntrials
            Ytrain{i,j}(:,trial) = cellstr(num2str(mod(sum(Xtrain{i,j}(:,1:p_prime,trial)>0,2),2)));
        end
    end
    Xtest{j} = rand(ntest,p)*2 - 1;
    Ytest{j} = cellstr(num2str(mod(sum(Xtest{j}(:,1:p_prime)>0,2),2)));
end

save('~/Documents/MATLAB/Data/Sparse_parity_vary_n_data.mat','Xtrain','Ytrain',...
    'Xtest','Ytest','ns','ntest','ps','ntrials','-v7.3')