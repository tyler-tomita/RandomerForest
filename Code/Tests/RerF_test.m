% train and test classifiers on abalone dataset

clear
close all
clc

nTrials = 10;

% FilePath = '/scratch/groups/jvogels3/tyler/Data/tests/';
OutFile = '/scratch/groups/jvogels3/tyler/RandomerForest/Results/test/test_results_2017_08_01.mat';
FilePath = '~/tests/';

Datasets = {'Sparse_parity';'Trunk';'Orthant';'mnist'};

for d = 1:length(Datasets)
    D = Datasets{d};

    TrainFile = [FilePath D '_train.csv'];
    TestFile = [FilePath D '_test.csv'];

    Xtrain = dlmread(TrainFile);
    Ytrain = cellstr(num2str(Xtrain(:,end)));
    Xtrain(:,end) = [];

    Xtest = dlmread(TestFile);
    Ytest = cellstr(num2str(Xtest(:,end)));
    Xtest(:,end) = [];

    [ntrain,p] = size(Xtrain);
    ntest = length(Ytest);

    Labels = unique([Ytrain;Ytest]);
    nClasses = length(Labels);

    mtrys = ceil(p.^[1/2 2]);

    Params.(D).nTrees = 500;
    Params.(D).Stratified = true;
    Params.(D).NWorkers = 2;
    Params.(D).Rescale = 'off';
    Params.(D).ForestMethod = 'rerf';
    Params.(D).RandomMatrix = 'binary-redundant';
    Params.(D).d = mtrys;
    Params.(D).rho = 1/p;
    Params.(D).Rotate = false;

    Train_Time = NaN(nTrials,length(Params.(D).d)*length(Params.(D).rho));
    OOB_Error = NaN(nTrials,length(Params.(D).d)*length(Params.(D).rho));
    Num_Nodes = NaN(nTrials,length(Params.(D).d)*length(Params.(D).rho));
    Test_Error = NaN(nTrials,length(Params.(D).d)*length(Params.(D).rho));

    for trial = 1:nTrials

        % train classifier
        poolobj = gcp('nocreate');
        if isempty(poolobj)
            parpool('local',Params.(D).NWorkers,...
                'IdleTimeout',360);
        end

        [Forest,~,Train_Time(trial,:)] = ...
            RerF_train(Xtrain,Ytrain,Params.(D));

        fprintf('Training complete\n')

        % evaluate classifier

        for k = 1:length(Forest)
            Labels = Forest{k}.classname;
            nClasses = length(Labels);
            Scores = rerf_oob_classprob(Forest{k},Xtrain,'last');
            Predictions = predict_class(Scores,Labels);
            OOB_Error(trial,k) = misclassification_rate(Predictions,Ytrain,false);
            NN = NaN(Forest{k}.nTrees,1);
            Trees = Forest{k}.Tree;
            parfor kk = 1:Forest{k}.nTrees
                NN(kk) = Trees{kk}.numnodes;
            end
            Num_Nodes(trial,k) = mean(NN);

            if ~strcmp(Forest{k}.Rescale,'off')
                TestScores = rerf_classprob(Forest{k},Xtest,'last',Xtrain);
            else
                TestScores = rerf_classprob(Forest{k},Xtest,'last');
            end
            TestPredictions = predict_class(TestScores,Forest{k}.classname);
            Test_Error(trial,k) = misclassification_rate(TestPredictions,Ytest,false);
        end

        clear Forest
    end

    TrainTime.(D).min = min(Train_Time);
    TrainTime.(D).max = max(Train_Time);
    TrainTime.(D).mean = mean(Train_Time);
    TrainTime.(D).std = std(Train_Time);

    OOBError.(D).min = min(OOB_Error);
    OOBError.(D).max = max(OOB_Error);
    OOBError.(D).mean = mean(OOB_Error);
    OOBError.(D).std = std(OOB_Error);

    TestError.(D).min = min(Test_Error);
    TestError.(D).max = max(Test_Error);
    TestError.(D).mean = mean(Test_Error);
    TestError.(D).std = std(Test_Error);

    NumNodes.(D).min = min(Num_Nodes);
    NumNodes.(D).max = max(Num_Nodes);
    NumNodes.(D).mean = mean(Num_Nodes);
    NumNodes.(D).std = std(Num_Nodes);
end

% save(OutFile,'Params','TrainTime','OOBError','TestError','NumNodes')

delete(gcp);
