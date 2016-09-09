close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

rng(1);

ntrain = 1000;
ntest = 10000;
dims = [2 5 10 20 40];
ndims = length(dims);
ntrials = 10;
OutlierProportion = 0.2;

for i = 1:ndims
    p = dims(i);
    fprintf('d = %d\n',p)
    p_prime = min(3,p);
    Sigma = 1/32*ones(1,p);
    Xtrain(i).Untransformed = rand(ntrain,p,ntrials)*2 - 1;
    Ytrain(i).Untransformed = cell(ntrain,ntrials);
    for trial = 1:ntrials
        Ytrain(i).Untransformed(:,trial) = cellstr(num2str(mod(sum(Xtrain(i).Untransformed(:,1:p_prime,trial)>0,2),2)));
    end
    Xtest(i).Untransformed = repmat(rand(ntest,p)*2 - 1,1,1,ntrials);
    for trial = 1:ntrials
        R = random_rotation(p);
        S = 10.^zeros(ntrain+ntest,p);
        Exponents = [-5 0 5];
        for j = 1:p_prime
            S(:,j) = 10^Exponents(j);
        end
        Ytest(i).Untransformed(:,trial) = cellstr(num2str(mod(sum(Xtest(i).Untransformed(:,1:p_prime,trial)>0,2),2)));
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

save('~/Documents/MATLAB/Data/Sparse_parity_data.mat','Xtrain','Ytrain',...
    'Xtest','Ytest','ntrain','ntest','dims','ntrials','p_prime','-v7.3')