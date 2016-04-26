close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

rng(1);

n.train = 1000;
n.test = 1000;
dims = [2 5 10 25 50 100];
ndims = length(dims);
n.trainSets = 20;
X.train = cell(1,ndims);
Y.train = cell(1,ndims);
X.test = cell(1,ndims);
Y.test = cell(1,ndims);
TestPosteriors = cell(1,ndims);

for i = 1:ndims
    d = dims(i);
    dgood = min(3,d);
    Sigma = 1/32*eye(d);
    x = zeros(n.train,d,n.trainSets);
    y = cell(n.train,n.trainSets);
    for trainSet = 1:n.trainSets
        Mu = sparse(n.train,d);
        for j = 1:n.train
            Mu(j,:) = binornd(1,0.5,1,d);
            x(j,:,trainSet) = mvnrnd(Mu(j,:),Sigma);
        end
        nones = sum(Mu(:,1:dgood),2);
        y(:,trainSet) = cellstr(num2str(mod(nones,2)));
    end
    X.train{i} = x;
    Y.train{i} = y;
    Mu = zeros(n.test,d);
    X.test{i} = zeros(n.test,d);
    for j = 1:n.test
        Mu(j,:) = binornd(1,0.5,1,d);
        X.test{i}(j,:) = mvnrnd(Mu(j,:),Sigma);
    end
    nones = sum(Mu(:,1:dgood),2);
    Y.test{i} = cellstr(num2str(mod(nones,2)));
    Mu = unique(Mu,'rows');
    nones = sum(Mu(:,1:dgood),2);
    Labels = mod(nones,2)';
    Sigma = repmat(Sigma,1,1,size(Mu,1));
    ComponentPosteriors = gmm_class_post(X.test{i},Mu,Sigma);
    UniqueLabels = unique(Labels);
    for l = 1:length(UniqueLabels)
        TestPosteriors{i}(:,l) = sum(ComponentPosteriors(:,Labels==UniqueLabels(l)),2);
    end
end

save([rerfPath 'RandomerForest/Data/Sparse_parity_partitioned_data.mat'],'X','Y',...
    'TestPosteriors','n','dims')