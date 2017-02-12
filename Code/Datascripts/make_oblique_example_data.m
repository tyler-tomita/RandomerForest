close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

rng(1);

ns = [25,1000];
p = 2;

for i = 1:length(ns)
    ntrain = ns(i);
    
    Xtrain{i} = rand(ntrain,p)*2 - 1;
    Ytrain{i} = cellstr(num2str(Xtrain{i}(:,2) - Xtrain{i}(:,1) > 0));
end

save('~/Documents/MATLAB/Data/Oblique_example_data.mat','Xtrain','Ytrain',...
    'ns','p')