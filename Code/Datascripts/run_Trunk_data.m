close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

rng(1);

ntrain = 100;
ntest = 10000;
dims = [2 10 50 100 500 1000];
ndims = length(dims);
ntrials = 10;
OutlierProportion = 0.2;
Class = [0;1];

for i = 1:ndims
    p = dims(i);
    p_idx = 1:p;
    mu1 = 1./sqrt(p_idx);
    mu0 = -1*mu1;
    Mu = cat(1,mu0,mu1);
    Sigma = ones(1,p);
    obj = gmdistribution(Mu,Sigma);
    Xtrain(i).Untransformed = zeros(ntrain,p,ntrials);
    Ytrain(i).Untransformed = cell(ntrain,ntrials);
    for trial = 1:ntrials
        [Xtrain(i).Untransformed(:,:,trial),idx] = random(obj,ntrain);
        Ytrain(i).Untransformed(:,trial) = cellstr(num2str(Class(idx)));
    end
    [Xtest(i).Untransformed,idx] = random(obj,ntest);
    Xtest(i).Untransformed = repmat(Xtest(i).Untransformed,1,1,ntrials);
    for trial = 1:ntrials
        R = random_rotation(p);
        S = 10.^random_scaling(ntrain+ntest,p,-5,5);
        Ytest(i).Untransformed(:,trial) = cellstr(num2str(Class(idx)));
        Xtrain(i).Rotated(:,:,trial) = Xtrain(i).Untransformed(:,:,trial)*R;
        Xtest(i).Rotated(:,:,trial) = Xtest(i).Untransformed(:,:,trial)*R;
        Ytrain(i).Rotated(:,trial) = Ytrain(i).Untransformed(:,trial);
        Ytest(i).Rotated(:,trial) = Ytest(i).Untransformed(:,trial);
        Xtrain(i).Scaled(:,:,trial) = Xtrain(i).Untransformed(:,:,trial).*S(1:ntrain,:);
        Xtest(i).Scaled(:,:,trial) = Xtest(i).Untransformed(:,:,trial).*S(ntrain+1:end,:);
        Ytrain(i).Scaled(:,trial) = Ytrain(i).Untransformed(:,trial);
        Ytest(i).Scaled(:,trial) = Ytest(i).Untransformed(:,trial);
        Xtrain(i).Affine(:,:,trial) = (Xtrain(i).Untransformed(:,:,trial)*R).*S(1:ntrain,:);
        Xtest(i).Affine(:,:,trial) = (Xtest(i).Untransformed(:,:,trial)*R).*S(ntrain+1:end,:);
        Ytrain(i).Affine(:,trial) = Ytrain(i).Untransformed(:,trial);
        Ytest(i).Affine(:,trial) = Ytest(i).Untransformed(:,trial);
        nOutlier = ceil(OutlierProportion*numel(Xtrain(i).Untransformed(:,:,trial)));
        Xtrain(i).Outlier(:,:,trial) = Xtrain(i).Untransformed(:,:,trial);
        Ytrain(i).Outlier(:,trial) = Ytrain(i).Untransformed(:,trial);
        Exponents = randsample(2:5,nOutlier,true)';
        OutlierIdx = randperm(numel(Xtrain(i).Untransformed(:,:,trial)),nOutlier)' + numel(Xtrain(i).Untransformed(:,:,trial))*(trial-1);
        Xtrain(i).Outlier(OutlierIdx) = Xtrain(i).Untransformed(OutlierIdx).*10.^(Exponents);
        Xtest(i).Outlier(:,:,trial) = Xtest(i).Untransformed(:,:,trial);
        Ytest(i).Outlier(:,trial) = Ytest(i).Untransformed(:,trial);
    end
end

save('~/Documents/MATLAB/Data/Trunk_data.mat','Xtrain','Ytrain',...
    'Xtest','Ytest','ntrain','ntest','dims','ntrials','-v7.3')