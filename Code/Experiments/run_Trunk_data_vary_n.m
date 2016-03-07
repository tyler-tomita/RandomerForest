close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

rng(1);

ns = [10 100 1000 10000];
d = 100;
ntrials = 5;
Class = [0;1];
X = cell(1,length(ns));
Y = cell(1,length(ns));
d_idx = 1:d;
mu1 = 1./sqrt(d_idx);
mu0 = -1*mu1;
Mu = cat(1,mu0,mu1);
Sigma = ones(1,d);
obj = gmdistribution(Mu,Sigma);

for i = 1:length(ns)
    n = ns(i);
    x = zeros(n,d,ntrials);
    y = cell(n,ntrials);
    for trial = 1:ntrials
        [x(:,:,trial),idx] = random(obj,n);
        y(:,trial) = cellstr(num2str(Class(idx)));
    end
    X{i} = x;
    Y{i} = y;
end

save([rerfPath 'RandomerForest/Results/Trunk_data_vary_n.mat'],'X','Y','ns','d','ntrials')