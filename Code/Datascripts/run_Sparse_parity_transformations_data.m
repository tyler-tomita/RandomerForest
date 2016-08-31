close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

rng(1);

load Sparse_parity_data

for i = 1:length(Xtrain)
    S(i).Untransformed = Xtrain{i};
    SS(i).Untransformed = Xtest{i};
    T(i).Untransformed = Ytrain{i};
    TT(i).Untransformed = Ytest{i};
end
Xtrain = S;
Xtest = SS;
Ytrain = T;
Ytest = TT;

n_out = round(0.2*ntrain);

for i = 1:length(Xtrain)
    dgood = min(3,dims(i));
    for trial = 1:ntrials
        R = random_rotation(dims(i));
        S = 10.^zeros(ntrain+ntest,dims(i));
        Exponents = [-5 0 5];
        for j = 1:dgood
            S(:,j) = 10^Exponents(j);
        end
        Xtrain(i).Rotated(:,:,trial) = Xtrain(i).Untransformed(:,:,trial)*R;
        Xtest(i).Rotated(:,:,trial) = Xtest(i).Untransformed*R;
        Ytrain(i).Rotated(:,trial) = Ytrain(i).Untransformed(:,trial);
        Ytest(i).Rotated(:,trial) = Ytest(i).Untransformed;
        Xtrain(i).Scaled(:,:,trial) = Xtrain(i).Untransformed(:,:,trial).*S(1:ntrain,:);
        Xtest(i).Scaled(:,:,trial) = Xtest(i).Untransformed.*S(ntrain+1:end,:);
        Ytrain(i).Scaled(:,trial) = Ytrain(i).Untransformed(:,trial);
        Ytest(i).Scaled(:,trial) = Ytest(i).Untransformed;
        Xtrain(i).Affine(:,:,trial) = (Xtrain(i).Untransformed(:,:,trial)*R).*S(1:ntrain,:);
        Xtest(i).Affine(:,:,trial) = (Xtest(i).Untransformed*R).*S(ntrain+1:end,:);
        Ytrain(i).Affine(:,trial) = Ytrain(i).Untransformed(:,trial);
        Ytest(i).Affine(:,trial) = Ytest(i).Untransformed;
        OutlierProportion = 0.1;
        nOutlier = ceil(OutlierProportion*numel(Xtrain(i).Untransformed(:,:,trial)));
        Xtrain(i).Outlier(:,:,trial) = Xtrain(i).Untransformed(:,:,trial);
        Ytrain(i).Outlier(:,trial) = Ytrain(i).Untransformed(:,trial);
        Exponents = randsample(2:5,nOutlier,true)';
        OutlierIdx = randperm(numel(Xtrain(i).Untransformed(:,:,trial)),nOutlier)' + numel(Xtrain(i).Untransformed(:,:,trial))*(trial-1);
        Xtrain(i).Outlier(OutlierIdx) = Xtrain(i).Untransformed(OutlierIdx).*10.^(Exponents);
        Xtest(i).Outlier(:,:,trial) = Xtest(i).Untransformed;
        Ytest(i).Outlier(:,trial) = Ytest(i).Untransformed;
    end
end

save('~/Documents/MATLAB/Data/Sparse_parity_transformations_data.mat','Xtrain',...
    'Ytrain','Xtest','Ytest','-v7.3')