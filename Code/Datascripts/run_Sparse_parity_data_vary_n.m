close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

rng(1);

ns = [10 100 1000 10000];
d = 5;
dgood = 3;
ntrials = 5;
X = cell(1,length(ns));
Y = cell(1,length(ns));
Sigma = 1/32*ones(1,d);

for i = 1:length(ns)
    n = ns(i);
    x = zeros(n,d,ntrials);
    y = cell(n,ntrials);
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

save([rerfPath 'RandomerForest/Results/Sparse_parity_data_vary_n.mat'],'X','Y','ns','d','ntrials')