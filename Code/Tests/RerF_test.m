% train and test clas

clear
close all
clc

nTrials = 20;

FilePath = '/scratch/groups/jvogels3/tyler/Data/tests/';
OutFile = '/scratch/groups/jvogels3/tyler/RandomerForest/Results/Tests/test_results_2017_08_01.mat';

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

    if strcmp(D, 'mnist')
        mtrys = ceil([p^(1/2) 10*p]);
    else
        mtrys = ceil(p.^[1/2 2]);
    end

    Params.(D).nTrees = 500;
    Params.(D).Stratified = true;
    Params.(D).NWorkers = 16;
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

    TrainTime.(D) = Train_Time;
    OOBError.(D) = OOB_Error;
    TestError.(D) = Test_Error;
    NumNodes.(D) = Num_Nodes;
end

save(OutFile,'Params','TrainTime','OOBError','TestError','NumNodes')

delete(gcp);
