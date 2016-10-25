% Multiclass problem where a different class lies in each quadrant of the
% hypercube

close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

rng(1);

ns = [100,200,400,800,1600];
ntest = 10000;
p = 6;
ntrials = 10;
XSign = double(all_binary_sets(p));
XSign(XSign==0) = -1;
nLabels = size(XSign,1);
Sigma = ones(1,p);
Lambda = 1/nLabels; % uniform class priors

% sample training data
for i = 1:length(ns)
    ntrain = ns(i);
    nSample = ceil(ntrain*Lambda);  % size of random sample per class
    fprintf('ntrain = %d\n',ntrain)
    Xtrain{i} = zeros(ntrain,p,ntrials);
    Ytrain{i} = cell(ntrain,ntrials);
    for trial = 1:ntrials
        x = zeros(nSample*nLabels,p);
        y = zeros(nSample*nLabels,1);
        for j = 1:nLabels
            x((1:nSample)+(j-1)*nSample,:) = rand(nSample,p).*repmat(XSign(j,:),nSample,1);
            y((1:nSample)+(j-1)*nSample) = j;
        end
        KeepIdx = randperm(nSample*nLabels,ntrain);
        Xtrain{i}(:,:,trial) = x(KeepIdx,:);
        Ytrain{i}(:,trial) = cellstr(num2str(y(KeepIdx)));
    end
end

% sample test data
nSample = ceil(ntest*Lambda);
x = zeros(nSample*nLabels,p);
y = zeros(nSample*nLabels,1);
for j = 1:nLabels
    x((1:nSample)+(j-1)*nSample,:) = rand(nSample,p).*repmat(XSign(j,:),nSample,1);
    y((1:nSample)+(j-1)*nSample) = j;
end
KeepIdx = randperm(nSample*nLabels,ntrain);
Xtest = x(KeepIdx,:);
Ytest = cellstr(num2str(y(KeepIdx)));

save('~/Documents/MATLAB/Data/Multiclass_hypercube_data.mat','Xtrain','Ytrain',...
    'Xtest','Ytest','ns','ntest','p','ntrials')