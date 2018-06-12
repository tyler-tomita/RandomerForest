clear
close all
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

rng(1);

% load data
dataSet = 'abalone';
datapath = '~/tmp/uci/processed/';
fold = getFolds([datapath,'cv_partitions/',dataSet,'_partitions.txt']);
nFolds = length(fold);
X = dlmread([datapath,'data/',dataSet,'.csv'],',');
Y = cellstr(num2str(X(:,end)));
X(:,end) = [];
X = zscore(X); % standardize columns
[n,p] = size(X);
if exist([datapath,'categorical_map/',dataSet,'_catmap.txt'],'file')
    iFeatureNum = getCatMap([datapath,'categorical_map/',dataSet,'_catmap.txt'],X);
else
    iFeatureNum = [];
end


% set up parameters
nTrees = 500;
options = optionsClassCCF;

if p <= 4
    lambdas = 1:p;
else
    lambdas = ceil(p.^([1/4, 1/2, 3/4, 1]));
end

nModels = length(lambdas);

testError = NaN(1,nFolds);

% iterate over partitions
for k = 1:nFolds
    fprintf('fold %d\n',k)
    trainIdx = cell2mat(fold((1:nFolds)~=k));
    testIdx = fold{k};
    ntrain = length(trainIdx);
    ntest = n - ntrain;

    % estimate classification performance of  lambda values using CV
    if ntrain < 10
        c = cvpartition(n,'LeaveOut');
    else
        c = cvpartition(Y(trainIdx),'KFold',10,'Stratify',true);
    end
    
    
    cvError = zeros(1,nModels);
    
    for i = 1:nModels
        options.lambda = lambdas(i);
        cvfun = @(xtrain,ytrain,xtest,ytest)...
            (mean(~strcmp(CCFPredict(xtrain,ytrain,xtest,nTrees,options,iFeatureNum),ytest)));
        cvError(i) = mean(crossval(cvfun,X(trainIdx,:),Y(trainIdx),'partition',c));
    end
    
    minError = min(cvError);
    minIdx = find(cvError==minError);
    nmin = length(minIdx);
    if nmin > 1
        minIdx = minIdx(randperm(nmin,1));
    end
    options.lambda = lambdas(minIdx);
    
    yhats = CCFPredict(X(trainIdx,:),Y(trainIdx),X(testIdx,:),nTrees,options,iFeatureNum);
    testError(k) = mean(~strcmp(yhats,Y(testIdx)));
end

save([rerfPath 'RandomerForest/Results/2018.06.12/' dataSet '_ccf_2018_06_12.mat'],'testError')
