close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

rng(1);

n = 100;
dims = [2 10 50 100 500 1000];
ndims = length(dims);
ntrials = 25;
Class = [0;1];
X = cell(1,ndims);
Y = cell(1,ndims);

for i = 1:ndims
    d = dims(i);
    d_idx = 1:d;
    mu1 = 1./sqrt(d_idx);
    mu0 = -1*mu1;
    Mu = cat(1,mu0,mu1);
    Sigma = ones(1,d);
    obj = gmdistribution(Mu,Sigma);
    x = zeros(n,d,ntrials);
    y = cell(n,ntrials);
    for trial = 1:ntrials
        [x(:,:,trial),idx] = random(obj,n);
        y(:,trial) = cellstr(num2str(Class(idx)));
    end
    X{i} = x;
    Y{i} = y;
end

save([rerfPath 'RandomerForest/Results/Trunk_data.mat'],'X','Y','n','dims','ntrials')