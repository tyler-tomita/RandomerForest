close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

rng(1);

n = 1000;
dims = [2 5 10 25 50 100];
ndims = length(dims);
ntrials = 25;
Class = [0;1];
X = cell(1,ndims);
Y = cell(1,ndims);

for i = 1:ndims
    d = dims(i);
    dgood = min(3,d);
    x = zeros(n,d,ntrials);
    y = cell(n,ntrials);
    Sigma = 1/32*ones(1,d);
    for trial = 1:ntrials
        Mu = sparse(n,d);
        for jj = 1:n
            Mu(jj,:) = binornd(1,0.5,1,d);
            x(jj,:,trial) = mvnrnd(Mu(jj,:),Sigma);
        end
        nones = sum(Mu(:,1:dgood),2);
        y(:,trial) = cellstr(num2str(mod(nones,2)));
    end
    X{i} = x;
    Y{i} = y;
end

save([rerfPath 'RandomerForest/Results/Sparse_parity_data.mat'],'X','Y','n','dims','ntrials')