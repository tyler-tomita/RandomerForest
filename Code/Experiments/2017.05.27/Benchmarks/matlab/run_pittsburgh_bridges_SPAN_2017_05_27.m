% train and test classifiers on pittsburgh_bridges_SPAN dataset

close all
clear
clc

Classifiers = {'rf','rerf','rerfr','rerfp','rerfpr','frc','frcr','rr_rf','rr_rfr'};

TrainFile = '/scratch/groups/jvogels3/tyler/Benchmarks/Data/uci/processed/pittsburgh_bridges_SPAN.train.csv';
TestFile = '/scratch/groups/jvogels3/tyler/Benchmarks/Data/uci/processed/pittsburgh_bridges_SPAN.test.csv';
OutFile = '/scratch/groups/jvogels3/tyler/RandomerForest/Results/2017.05.27/Benchmarks/pittsburgh_bridges_SPAN_2017_05_27.mat';

rng(1);

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

if p <= 5
    mtrys = [1:p p.^[2 3]];
elseif p <= 10 && ntrain <= 10000
    mtrys = ceil(p.^[1/4 1/2 3/4 1 2 3]);
elseif p <= 100 && ntrain <= 2000
    mtrys = ceil(p.^[1/4 1/2 3/4 1 2 2.5]);
elseif p <= 100 && ntrain <= 10000
    mtrys = ceil(p.^[1/4 1/2 3/4 1 2]);
else
    mtrys = [ceil(p.^[1/4 1/2 3/4 1]) 20*p];
end
mtrys_rf = mtrys(mtrys<=p);

for c = 1:length(Classifiers)
    cl = Classifiers{c};
    fprintf('%s start\n',cl)

    if ntrain <= 2000
        Params.(cl).nTrees = 1000;
    else
        Params.(cl).nTrees = 500;
    end
    Params.(cl).Stratified = true;
    Params.(cl).NWorkers = 12;
    if strcmp(cl,'rerfr') || strcmp(cl,'rerfpr') || strcmp(cl,'frcr') || ...
            strcmp(cl,'rr_rfr')
        Params.(cl).Rescale = 'rank';
    else
        Params.(cl).Rescale = 'off';
    end
    Params.(cl).mdiff = 'off';
    if strcmp(cl,'rf')
        Params.(cl).ForestMethod = 'rf';
        Params.(cl).d = mtrys_rf;
    elseif strcmp(cl,'rerf') || strcmp(cl,'rerfr')
        Params.(cl).ForestMethod = 'rerf';
        Params.(cl).RandomMatrix = 'binary';
        Params.(cl).d = mtrys;
        Params.(cl).rho = 1/p;
%     elseif strcmp(cl,'rerfb') || strcmp(cl,'rerfbr')
%         Params.(cl).ForestMethod = 'rerf';
%         Params.(cl).RandomMatrix = 'uniform-nnzs-binary';
%         Params.(cl).d = mtrys;
%         Params.(cl).nmix = 1:5;
%         Params.(cl).rho = 1/p;
%     elseif strcmp(cl,'rerfc') || strcmp(cl,'rerfcr')
%         Params.(cl).ForestMethod = 'rerf';
%         Params.(cl).RandomMatrix = 'uniform-nnzs-continuous';
%         Params.(cl).d = mtrys;
%         Params.(cl).nmix = 1:5;
%         Params.(cl).rho = 1/p;
    elseif strcmp(cl,'rerfp') || strcmp(cl,'rerfpr')
        Params.(cl).ForestMethod = 'rerf';
        Params.(cl).RandomMatrix = 'poisson';
        Params.(cl).d = mtrys;
        Params.(cl).lambda = 2:5;            
    elseif strcmp(cl,'frc') || strcmp(cl,'frcr')
        Params.(cl).ForestMethod = 'rerf';
        Params.(cl).RandomMatrix = 'frc';
        Params.(cl).d = mtrys;
        Params.(cl).L = 2:5;       
    elseif strcmp(cl,'rr_rf') || strcmp(cl,'rr_rfr')
        Params.(cl).ForestMethod = 'rf';   
        Params.(cl).Rotate = true;
        Params.(cl).d = mtrys_rf;
    end
    
    if strcmp(Params.(cl).ForestMethod,'rf')
        OOBError.(cl) = NaN(1,length(Params.(cl).d));
        OOBAUC.(cl) = NaN(1,length(Params.(cl).d));
        TrainTime.(cl) = NaN(1,length(Params.(cl).d));
        Depth.(cl) = NaN(Params.(cl).nTrees,length(Params.(cl).d));
        NumNodes.(cl) = NaN(Params.(cl).nTrees,length(Params.(cl).d));
        NumSplitNodes.(cl) = NaN(Params.(cl).nTrees,length(Params.(cl).d));
        TreeStrength.(cl) = NaN(1,length(Params.(cl).d));
        TreeDiversity.(cl) = NaN(1,length(Params.(cl).d));
        TestScores.(cl) = NaN(ntest,nClasses,length(Params.(cl).d));
        TestError.(cl) = NaN(1,length(Params.(cl).d));
    else
        if strcmp(Params.(cl).RandomMatrix,'frc')
            OOBError.(cl) = NaN(1,length(Params.(cl).d)*length(Params.(cl).L));
            OOBAUC.(cl) = NaN(1,length(Params.(cl).d)*length(Params.(cl).L));
            TrainTime.(cl) = NaN(1,length(Params.(cl).d)*length(Params.(cl).L));
            Depth.(cl) = NaN(Params.(cl).nTrees,length(Params.(cl).d)*length(Params.(cl).L));
            NumNodes.(cl) = NaN(Params.(cl).nTrees,length(Params.(cl).d)*length(Params.(cl).L));
            NumSplitNodes.(cl) = NaN(Params.(cl).nTrees,length(Params.(cl).d)*length(Params.(cl).L));
            TreeStrength.(cl) = NaN(1,length(Params.(cl).d)*length(Params.(cl).L));
            TreeDiversity.(cl) = NaN(1,length(Params.(cl).d)*length(Params.(cl).L));
            TestScores.(cl) = NaN(ntest,nClasses,length(Params.(cl).d)*length(Params.(cl).L));
            TestError.(cl) = NaN(1,length(Params.(cl).d)*length(Params.(cl).L));
        elseif strcmp(Params.(cl).RandomMatrix,'poisson')
            OOBError.(cl) = NaN(1,length(Params.(cl).d)*length(Params.(cl).lambda));
            OOBAUC.(cl) = NaN(1,length(Params.(cl).d)*length(Params.(cl).lambda));
            TrainTime.(cl) = NaN(1,length(Params.(cl).d)*length(Params.(cl).lambda));
            Depth.(cl) = NaN(Params.(cl).nTrees,length(Params.(cl).d)*length(Params.(cl).lambda));
            NumNodes.(cl) = NaN(Params.(cl).nTrees,length(Params.(cl).d)*length(Params.(cl).lambda));
            NumSplitNodes.(cl) = NaN(Params.(cl).nTrees,length(Params.(cl).d)*length(Params.(cl).lambda));
            TreeStrength.(cl) = NaN(1,length(Params.(cl).d)*length(Params.(cl).lambda));
            TreeDiversity.(cl) = NaN(1,length(Params.(cl).d)*length(Params.(cl).lambda));
            TestScores.(cl) = NaN(ntest,nClasses,length(Params.(cl).d)*length(Params.(cl).lambda));
            TestError.(cl) = NaN(1,length(Params.(cl).d)*length(Params.(cl).lambda));
        else
            OOBError.(cl) = NaN(1,length(Params.(cl).d)*length(Params.(cl).rho));
            OOBAUC.(cl) = NaN(1,length(Params.(cl).d)*length(Params.(cl).rho));
            TrainTime.(cl) = NaN(1,length(Params.(cl).d)*length(Params.(cl).rho));
            Depth.(cl) = NaN(Params.(cl).nTrees,length(Params.(cl).d)*length(Params.(cl).rho));
            NumNodes.(cl) = NaN(Params.(cl).nTrees,length(Params.(cl).d)*length(Params.(cl).rho));
            NumSplitNodes.(cl) = NaN(Params.(cl).nTrees,length(Params.(cl).d)*length(Params.(cl).rho));
            TreeStrength.(cl) = NaN(1,length(Params.(cl).d)*length(Params.(cl).rho));
            TreeDiversity.(cl) = NaN(1,length(Params.(cl).d)*length(Params.(cl).rho));
            TestScores.(cl) = NaN(ntest,nClasses,length(Params.(cl).d)*length(Params.(cl).rho));
            TestError.(cl) = NaN(1,length(Params.(cl).d)*length(Params.(cl).rho));
        end
    end

    % train classifier
    poolobj = gcp('nocreate');
    if isempty(poolobj)
        parpool('local',Params.(cl).NWorkers,...
            'IdleTimeout',360);
    end

    [Forest,~,TrainTime.(cl)] = ...
        RerF_train(Xtrain,Ytrain,Params.(cl));

    fprintf('Training complete\n')

    % compute oob auc, oob error, and tree stats

    for k = 1:length(Forest)
        Labels = Forest{k}.classname;
        nClasses = length(Labels);
        Scores = rerf_oob_classprob(Forest{k},...
            Xtrain,'last');
        Predictions = predict_class(Scores,Labels);
        OOBError.(cl)(k) = ...
            misclassification_rate(Predictions,Ytrain,...
        false);        
        if nClasses > 2
            Yb = binarize_labels(Ytrain,Labels);
            [~,~,~,OOBAUC.(cl)(k)] = ... 
                perfcurve(Yb(:),Scores(:),'1');
        else
            [~,~,~,OOBAUC.(cl)(k)] = ...
                perfcurve(Ytrain,Scores(:,2),'1');
        end
        Depth.(cl)(:,k) = forest_depth(Forest{k})';
        NN = NaN(Forest{k}.nTrees,1);
        NS = NaN(Forest{k}.nTrees,1);
        Trees = Forest{k}.Tree;
        parfor kk = 1:Forest{k}.nTrees
            NN(kk) = Trees{kk}.numnodes;
            NS(kk) = sum(Trees{kk}.isbranch);
        end
        NumNodes.(cl)(:,k) = NN;
        NumSplitNodes.(cl)(:,k) = NS;

        if ~strcmp(Forest{k}.Rescale,'off')
            Scores = rerf_classprob(Forest{k},Xtest,'individual',true,Xtrain);
        else
            Scores = rerf_classprob(Forest{k},Xtest,'individual',true);
        end
        PredCell = cell(ntest,Params.(cl).nTrees);
        parfor kk = 1:Params.(cl).nTrees
            PredCell(:,kk) = predict_class(Scores(:,:,kk),Labels);
        end

        TreeStrength.(cl)(k) = 1 - misclassification_rate(PredCell,Ytest,true);
        TreeDiversity.(cl)(k) = classifier_variance(PredCell);
        
        if ~strcmp(Forest{k}.Rescale,'off')
            TestScores.(cl)(:,:,k) = rerf_classprob(Forest{k},Xtest,'last',true,Xtrain);
        else
        TestScores.(cl)(:,:,k) = rerf_classprob(Forest{k},Xtest,true,'last');
        end
        TestPredictions = predict_class(TestScores.(cl)(:,:,k),Forest{k}.classname);

        TestError.(cl)(k) = ...
            misclassification_rate(TestPredictions,Ytest,false);
    end

    % select best model based on OOB errors and AUCs
    BI = hp_optimize(OOBError.(cl)(end,:),...
        OOBAUC.(cl)(end,:));
    BestIdx.(cl) = BI(randperm(length(BI),1));

    Forest = Forest{BestIdx.(cl)};
    
    save(OutFile,'Params','OOBError','OOBAUC','TestError',...
        'TrainTime','Depth','NumNodes','NumSplitNodes','TreeStrength',...
        'TreeDiversity','BestIdx','TestScores','Forest')
    
    clear Forest
end  
