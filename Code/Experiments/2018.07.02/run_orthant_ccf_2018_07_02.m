clear
close all
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

rng(1);

ps = [2 4 6];
ns = [20 200 400;80 400 4000;400 2000 4000]';
ntest = 10000;
nTrials = 10;
testError = NaN(length(ns(:,1)),length(ps),nTrials);

datapath= '~/work/tyler/R/Data/Orthant/dat/Raw/';

% set up parameters
nTrees = 500;
options = optionsClassCCF;
iFeatureNum = [];

for j = 1:length(ps)
    
    p = ps(j);
    
    Xtest = dlmread([datapath 'Test/Orthant_raw_test_set_p' p '.dat']);
    Ytest = cellstr(num2str(Xtest(:,end)));
    Xtest(:,end) = [];

    if p <= 4
        lambdas = 1:p;
    else
        lambdas = ceil(p.^([1/4, 1/2, 3/4, 1]));
    end

    nModels = length(lambdas);
    
    for i = 1:length(ns(:,j))
        ntrain = ns(i,j);
        
        for trial = 1:nTrials
            
            Xtrain = dlmread([datapath 'Train/Orthant_raw_train_set_p' p '_n' ntrain '_trial' trial '.dat']);
            Ytrain = cellstr(num2str(Xtrain(:,end)));
            Xtrain(:,end) = [];

        
            % initialize parallel pool
            poolobj = gcp('nocreate');
            if isempty(poolobj)
                parpool('local',20,'IdleTimeout',360);
            end

            % estimate classification performance of lambda values using CV
            if ntrain <= 10
                c = cvpartition(ntrain,'LeaveOut');
            else
                c = cvpartition(Ytrain,'KFold',10,'Stratify',true);
            end
    
            cvError = zeros(1,nModels);
    
            for k = 1:nModels
                options.lambda = lambdas(k);
                cvfun = @(xtrain,ytrain,xtest,ytest)...
                    (mean(~strcmp(CCFPredict(xtrain,ytrain,xtest,nTrees,options,iFeatureNum),ytest)));
                cvError(k) = mean(crossval(cvfun,X(trainIdx,:),Y(trainIdx),'partition',c));
            end
    
            minError = min(cvError);
            minIdx = find(cvError==minError);
            nmin = length(minIdx);
            if nmin > 1
                minIdx = minIdx(randperm(nmin,1));
            end
            options.lambda = lambdas(minIdx);
    
            Yhats = CCFPredict(Xtrain,Ytrain,Xtest,nTrees,options,iFeatureNum);
            testError(i,j,trial) = mean(~strcmp(Yhats,Ytest));
            save([rerfPath 'RandomerForest/Results/2018.07.02/Orthant_ccf_2018_07_02.mat'],'testError')
        end
    end
end

delete(gcp);
