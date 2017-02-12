close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

rng(1);

ntrain = 500;

Xtrain = rand(ntrain,2)*2 - 1;
Ytrain = zeros(ntrain,1);
Ytrain(Xtrain(:,1) > 0 & Xtrain(:,2) > 0) = 1;
Ytrain(Xtrain(:,1) < 0 & Xtrain(:,2) > 0) = 2;
Ytrain(Xtrain(:,1) < 0 & Xtrain(:,2) < 0) = 3;
Ytrain(Xtrain(:,1) > 0 & Xtrain(:,2) < 0) = 4;
Ytrain = cellstr(num2str(Ytrain));

save('~/Documents/MATLAB/Data/Orthant_data.mat','Xtrain','Ytrain')