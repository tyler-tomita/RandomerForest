% Generate Trunk datasets for various values of n and p

close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

rng(1);

ps = [10,100,500];
ns{1} = [10,100,1000,10000];
ns{2} = [10,100,1000,10000];
ns{3} = [10,100,1000,10000];
ntest = 10000;
ntrials = 10;
Class = [0;1];

Xtrain = cell(length(ns),length(ps));
Ytrain = cell(length(ns),length(ps));
Xtest = cell(1,length(ps));
Ytest = cell(1,length(ps));
% sample training data
for j = 1:length(ps)
    p = ps(j);
    fprintf('p = %d\n',p)
    p_idx = 1:p;
    mu1 = 1./sqrt(p_idx);
    mu0 = -1*mu1;
    Mu = cat(1,mu0,mu1);
    Sigma = ones(1,p);
    obj = gmdistribution(Mu,Sigma);
    for i = 1:length(ns{j})
        ntrain = ns{j}(i);
        fprintf('n = %d\n',ntrain)
        Xtrain{i,j} = zeros(ntrain,p,ntrials);
        Ytrain{i,j} = cell(ntrain,ntrials);
        for trial = 1:ntrials
            [Xtrain{i,j}(:,:,trial),idx] = random(obj,ntrain);
            Ytrain{i,j}(:,trial) = cellstr(num2str(Class(idx)));
        end
    end
    [Xtest{j},idx] = random(obj,ntest);
    Ytest{j} = cellstr(num2str(Class(idx)));
end

save('~/Documents/MATLAB/Data/Trunk_vary_n_data.mat','Xtrain','Ytrain',...
    'Xtest','Ytest','ns','ntest','ps','ntrials','-v7.3')