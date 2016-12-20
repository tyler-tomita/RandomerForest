close all
clear
clc

Xtrain = csvread('~/RandomerForest/Data/Parity/Parity_Xtrain.dat');
Ytrain = cellstr(num2str(csvread('~/RandomerForest/Data/Parity/Parity_Ytrain.dat')));
Xtest = csvread('~/RandomerForest/Data/Parity/Parity_Xtest.dat');
Ytest = csvread('~/RandomerForest/Data/Parity/Parity_Ytest.dat');

[ntest,p] = size(Xtest);

ntrials = 1;

rng(1);

Params.nTrees = 10000;
Params.NWorkers = 2;
Stratified = true;
Params.Rescale = 'off';
Params.mdiff = 'off';
Params.ForestMethod = 'rerf2';
Params.RandomMatrix = 'sparse-unadjusted';
Params.d = round(sqrt(p));
Params.dx = p;

for trial = 1:ntrials
    fprintf('Trial %d\n',trial)

    % train classifier
    poolobj = gcp('nocreate');
    if isempty(poolobj)
        parpool('local',Params.NWorkers,...
            'IdleTimeout',360);
    end

    Forest = RerF_train(Xtrain,Ytrain,Params);
    OOBScores = rerf_oob_classprob(Forest{1},Xtrain,'every');
    for k = 1:Forest{1}.nTrees
        Predictions = predict_class(OOBScores(:,:,k),Forest{1}.classname);
        OOBError(trial,k) = misclassification_rate(Predictions,Ytrain,false);
    end
    
%     TestScores = rerf_classprob(Forest{1},Xtest,'last');
%     TestPredictions(:,trial) = predict_class(TestScores,Forest{1}.classname);

end

plot(1:Params.nTrees,OOBError(1,:),'LineWidth',2)

xlabel('# of Trees')
ylabel('OOB Error')
title({'Parity';'RerF (A = mdim x mdim, mtry = mdim^{0.5})'})

% save_fig(gcf,'~/Parity_RerF2_OOB_Error')

% MajorityVote = mean(reshape(str2num(char(Predictions)),ntest,ntrials),2);

% save('~/RandomerForest/Data/Parity_predictons_100_round_average.mat','MajorityVote')