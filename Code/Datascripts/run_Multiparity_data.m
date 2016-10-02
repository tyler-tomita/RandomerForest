close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

rng(1);

ntrain = 1000;
ntest = 10000;
dims = [5 10 20 40];
ndims = length(dims);
ntrials = 10;
OutlierProportion = 0.2;
p_prime1 = 1:3;
p_prime2 = 4:5;
p_prime = [p_prime1 p_prime2];

for i = 1:ndims
    p = dims(i);
    fprintf('p = %d\n',p)
    Xtrain(i).Raw = rand(ntrain,p,ntrials)*2 - 1;
    Ytrain(i).Raw = cell(ntrain,ntrials);
    for trial = 1:ntrials
        y = zeros(ntrain,1);
        y(mod(sum(Xtrain(i).Raw(:,p_prime1,trial)>0,2),2)==0 & mod(sum(Xtrain(i).Raw(:,p_prime2,trial)>0,2),2)==0) = 1;
        y(mod(sum(Xtrain(i).Raw(:,p_prime1,trial)>0,2),2)==1 & mod(sum(Xtrain(i).Raw(:,p_prime2,trial)>0,2),2)==0) = 2;
        y(mod(sum(Xtrain(i).Raw(:,p_prime2,trial)>0,2),2)==1) = 3;
        Ytrain(i).Raw(:,trial) = cellstr(num2str(y));
    end
    Xtest(i).Raw = repmat(rand(ntest,p)*2 - 1,1,1,ntrials);
    S1{i} = 10.^zeros(ntrain+ntest,p);
    Exponents = [-5 -2.5 0 2.5 5];
    for j = 1:p_prime
        S1{i}(:,j) = 10^Exponents(j);
    end
    for trial = 1:ntrials
        R{i}(:,:,trial) = random_rotation(p);
        S2{i}(:,:,trial) = 10.^random_scaling(ntrain+ntest,p,-5,5);
        y = zeros(ntest,1);
        y(mod(sum(Xtest(i).Raw(:,p_prime1,trial)>0,2),2)==0 & mod(sum(Xtest(i).Raw(:,p_prime2,trial)>0,2),2)==0) = 1;
        y(mod(sum(Xtest(i).Raw(:,p_prime1,trial)>0,2),2)==1 & mod(sum(Xtest(i).Raw(:,p_prime2,trial)>0,2),2)==0) = 2;
        y(mod(sum(Xtest(i).Raw(:,p_prime2,trial)>0,2),2)==1) = 3;
        Ytest(i).Raw(:,trial) = cellstr(num2str(y));
        Xtrain(i).Rotated(:,:,trial) = Xtrain(i).Raw(:,:,trial)*R{i}(:,:,trial);
        Xtest(i).Rotated(:,:,trial) = Xtest(i).Raw(:,:,trial)*R{i}(:,:,trial);
        Ytrain(i).Rotated(:,trial) = Ytrain(i).Raw(:,trial);
        Ytest(i).Rotated(:,trial) = Ytest(i).Raw(:,trial);
        Xtrain(i).Scaled(:,:,trial) = Xtrain(i).Raw(:,:,trial).*S1{i}(1:ntrain,:);
        Xtest(i).Scaled(:,:,trial) = Xtest(i).Raw(:,:,trial).*S1{i}(ntrain+1:end,:);
        Ytrain(i).Scaled(:,trial) = Ytrain(i).Raw(:,trial);
        Ytest(i).Scaled(:,trial) = Ytest(i).Raw(:,trial);
        Xtrain(i).Affine(:,:,trial) = (Xtrain(i).Raw(:,:,trial)*R{i}(:,:,trial)).*S2{i}(1:ntrain,:,trial);
        Xtest(i).Affine(:,:,trial) = (Xtest(i).Raw(:,:,trial)*R{i}(:,:,trial)).*S2{i}(ntrain+1:end,:,trial);
        Ytrain(i).Affine(:,trial) = Ytrain(i).Raw(:,trial);
        Ytest(i).Affine(:,trial) = Ytest(i).Raw(:,trial);
        nOutlier = ceil(OutlierProportion*numel(Xtrain(i).Raw(:,:,trial)));
        Xtrain(i).Outlier(:,:,trial) = Xtrain(i).Raw(:,:,trial);
        Ytrain(i).Outlier(:,trial) = Ytrain(i).Raw(:,trial);
        Exponents = randsample(2:5,nOutlier,true)';
        OutlierIdx = randperm(numel(Xtrain(i).Raw(:,:,trial)),nOutlier)' + numel(Xtrain(i).Raw(:,:,trial))*(trial-1);
        Xtrain(i).Outlier(OutlierIdx) = Xtrain(i).Raw(OutlierIdx).*10.^(Exponents);
        Xtest(i).Outlier(:,:,trial) = Xtest(i).Raw(:,:,trial);
        Ytest(i).Outlier(:,trial) = Ytest(i).Raw(:,trial);
    end
end

save('~/Documents/MATLAB/Data/Multiparity_data.mat','Xtrain','Ytrain',...
    'Xtest','Ytest','ntrain','ntest','dims','ntrials','p_prime',...
    'p_prime1','p_prime2','S1','S2','R','-v7.3')