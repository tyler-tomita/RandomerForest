close all
clear
clc

Xtrain = csvread('~/RandomerForest/Data/Parity/Parity_Xtrain.dat');
Ytrain = csvread('~/RandomerForest/Data/Parity/Parity_Ytrain.dat');
Xtest = csvread('~/RandomerForest/Data/Parity/Parity_Xtest.dat');
Ytest = csvread('~/RandomerForest/Data/Parity/Parity_Ytest.dat');

[ntest,p] = size(Xtest);

ntrials = 100;

rng(1);

Params.nTrees = 1000;
Params.NWorkers = 2;
Stratified = true;
Params.Rescale = 'off';
Params.mdiff = 'off';
Params.ForestMethod = 'rerf2';
Params.d = round(sqrt(p));
Params.dx = round(sqrt(p));

for trial = 1:ntrials
    fprintf('Trial %d\n',trial)

    % train classifier
    poolobj = gcp('nocreate');
    if isempty(poolobj)
        parpool('local',Params.NWorkers,...
            'IdleTimeout',360);
    end

    tic;
    Forest = RerF_train(Xtrain,Ytrain,Params);
    Scores = rerf_classprob(Forest{1},Xtest,'last');
    Predictions(:,trial) = predict_class(Scores,Forest{1}.classname);

    clear Forest
end

MajorityVote = mean(reshape(str2num(char(Predictions)),ntest,ntrials),2);

save('~/RandomerForest/Data/Parity_predictons_100_round_average.mat','MajorityVote')