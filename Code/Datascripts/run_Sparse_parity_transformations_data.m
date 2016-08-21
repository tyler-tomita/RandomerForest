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
    Mu = [-1./sqrt(1:dims(i));1./sqrt(1:dims(i))];
    Sigma = 1/32*ones(1,dims(i));
    Sigma_outlier = 4*Sigma;
    for trial = 1:ntrials
        R = random_rotation(dims(i));
        S = 10.^random_scaling(ntrain+ntest,dims(i),-5,5);
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
        x_out = zeros(n_out,dims(i));
        Mu_out = zeros(n_out,dims(i));
        for jj = 1:n_out
            Mu_out(jj,:) = binornd(1,0.5,1,dims(i));
            x_out(jj,:) = mvnrnd(Mu_out(jj,:),Sigma_outlier);
        end
        Xtrain(i).Outlier(:,:,trial) = [Xtrain(i).Untransformed(:,:,trial);x_out];
        Xtest(i).Outlier(:,:,trial) = Xtest(i).Untransformed;
        nones_out = sum(Mu_out(:,1:dgood),2);
        y_out = cellstr(num2str(mod(nones_out,2)));
        Ytrain(i).Outlier(:,trial) = [Ytrain(i).Untransformed(:,trial);y_out];
        Ytest(i).Outlier(:,trial) = Ytest(i).Untransformed;
    end
end

save('~/Documents/MATLAB/Data/Sparse_parity_transformations_data.mat','Xtrain',...
    'Ytrain','Xtest','Ytest','-v7.3')