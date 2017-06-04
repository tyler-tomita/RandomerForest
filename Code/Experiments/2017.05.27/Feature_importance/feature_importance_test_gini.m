% train and test classifiers on sparse parity dataset

close all
clear
clc

Classifiers = {'rerfp'};

TrainFile = '~/R/Data/Sparse_parity/dat/Raw/Train/Sparse_parity_train_set_n1000_p10_trial1.dat';
TestFile = '~/R/Data/Sparse_parity/dat/Raw/Test/Sparse_parity_test_set_p10.dat';
OutFile = '~/RandomerForest/Results/2017.05.27/Feature_importance/Sparse_parity_feature_importance_gini_2017_05_27.mat';

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

% mtrys = ceil(p.^[1/2 1]);
mtrys = ceil(p.^[1/4 1/2 3/4 1 2]);
mtrys_rf = mtrys(mtrys<=p);

for c = 1:length(Classifiers)
    fprintf('%s start\n',Classifiers{c})

    Params.(Classifiers{c}).nTrees = 500;
    Params.(Classifiers{c}).Stratified = true;
    Params.(Classifiers{c}).NWorkers = 2;
    if strcmp(Classifiers{c},'rerfr') || strcmp(Classifiers{c},'rr_rfr')
        Params.(Classifiers{c}).Rescale = 'rank';
    else
        Params.(Classifiers{c}).Rescale = 'off';
    end
    Params.(Classifiers{c}).mdiff = 'off';
    if strcmp(Classifiers{c},'rf')
        Params.(Classifiers{c}).ForestMethod = 'rf';
        Params.(Classifiers{c}).d = mtrys_rf;
    elseif strcmp(Classifiers{c},'rerf') || strcmp(Classifiers{c},'rerfr')
        Params.(Classifiers{c}).ForestMethod = 'rerf';
        Params.(Classifiers{c}).RandomMatrix = 'binary';
        Params.(Classifiers{c}).d = mtrys;
        Params.(Classifiers{c}).rho = 1/p;
    elseif strcmp(Classifiers{c},'rerfb') || strcmp(Classifiers{c},'rerfbr')
        Params.(Classifiers{c}).ForestMethod = 'rerf';
        Params.(Classifiers{c}).RandomMatrix = 'uniform-nnzs-binary';
        Params.(Classifiers{c}).d = mtrys;
        Params.(Classifiers{c}).nmix = 1:5;
        Params.(Classifiers{c}).rho = 1/p;
    elseif strcmp(Classifiers{c},'rerfc') || strcmp(Classifiers{c},'rerfcr')
        Params.(Classifiers{c}).ForestMethod = 'rerf';
        Params.(Classifiers{c}).RandomMatrix = 'uniform-nnzs-continuous';
        Params.(Classifiers{c}).d = mtrys;
        Params.(Classifiers{c}).nmix = 1:5;
        Params.(Classifiers{c}).rho = 1/p;
    elseif strcmp(Classifiers{c},'rerfp') || strcmp(Classifiers{c},'rerfpr')
        Params.(Classifiers{c}).ForestMethod = 'rerf';
        Params.(Classifiers{c}).RandomMatrix = 'poisson';
        Params.(Classifiers{c}).d = mtrys;
        Params.(Classifiers{c}).lambda = 3;            
    elseif strcmp(Classifiers{c},'frc') || strcmp(Classifiers{c},'frcr')
        Params.(Classifiers{c}).ForestMethod = 'rerf';
        Params.(Classifiers{c}).RandomMatrix = 'frc';
        Params.(Classifiers{c}).d = mtrys;
        Params.(Classifiers{c}).L = 2:5;       
    elseif strcmp(Classifiers{c},'rr_rf') || strcmp(Classifiers{c},'rr_rfr')
        Params.(Classifiers{c}).ForestMethod = 'rf';   
        Params.(Classifiers{c}).Rotate = true;
        Params.(Classifiers{c}).d = mtrys_rf;
    end
    
    if strcmp(Params.(Classifiers{c}).ForestMethod,'rf')
        OOBError.(Classifiers{c}) = NaN(1,length(Params.(Classifiers{c}).d));
        OOBAUC.(Classifiers{c}) = NaN(1,length(Params.(Classifiers{c}).d));
        TrainTime.(Classifiers{c}) = NaN(1,length(Params.(Classifiers{c}).d));
        Depth.(Classifiers{c}) = NaN(Params.(Classifiers{c}).nTrees,length(Params.(Classifiers{c}).d));
        NumNodes.(Classifiers{c}) = NaN(Params.(Classifiers{c}).nTrees,length(Params.(Classifiers{c}).d));
        NumSplitNodes.(Classifiers{c}) = NaN(Params.(Classifiers{c}).nTrees,length(Params.(Classifiers{c}).d));
        TreeStrength.(Classifiers{c}) = NaN(1,length(Params.(Classifiers{c}).d));
        TreeDiversity.(Classifiers{c}) = NaN(1,length(Params.(Classifiers{c}).d));
        TestScores.(Classifiers{c}) = NaN(ntest,nClasses,length(Params.(Classifiers{c}).d));
        TestError.(Classifiers{c}) = NaN(1,length(Params.(Classifiers{c}).d));
    else
        if strcmp(Params.(Classifiers{c}).RandomMatrix,'frc')
            OOBError.(Classifiers{c}) = NaN(1,length(Params.(Classifiers{c}).d)*length(Params.(Classifiers{c}).L));
            OOBAUC.(Classifiers{c}) = NaN(1,length(Params.(Classifiers{c}).d)*length(Params.(Classifiers{c}).L));
            TrainTime.(Classifiers{c}) = NaN(1,length(Params.(Classifiers{c}).d)*length(Params.(Classifiers{c}).L));
            Depth.(Classifiers{c}) = NaN(Params.(Classifiers{c}).nTrees,length(Params.(Classifiers{c}).d)*length(Params.(Classifiers{c}).L));
            NumNodes.(Classifiers{c}) = NaN(Params.(Classifiers{c}).nTrees,length(Params.(Classifiers{c}).d)*length(Params.(Classifiers{c}).L));
            NumSplitNodes.(Classifiers{c}) = NaN(Params.(Classifiers{c}).nTrees,length(Params.(Classifiers{c}).d)*length(Params.(Classifiers{c}).L));
            TreeStrength.(Classifiers{c}) = NaN(1,length(Params.(Classifiers{c}).d)*length(Params.(Classifiers{c}).L));
            TreeDiversity.(Classifiers{c}) = NaN(1,length(Params.(Classifiers{c}).d)*length(Params.(Classifiers{c}).L));
            TestScores.(Classifiers{c}) = NaN(ntest,nClasses,length(Params.(Classifiers{c}).d)*length(Params.(Classifiers{c}).L));
            TestError.(Classifiers{c}) = NaN(1,length(Params.(Classifiers{c}).d)*length(Params.(Classifiers{c}).L));
        elseif strcmp(Params.(Classifiers{c}).RandomMatrix,'poisson')
            OOBError.(Classifiers{c}) = NaN(1,length(Params.(Classifiers{c}).d)*length(Params.(Classifiers{c}).lambda));
            OOBAUC.(Classifiers{c}) = NaN(1,length(Params.(Classifiers{c}).d)*length(Params.(Classifiers{c}).lambda));
            TrainTime.(Classifiers{c}) = NaN(1,length(Params.(Classifiers{c}).d)*length(Params.(Classifiers{c}).lambda));
            Depth.(Classifiers{c}) = NaN(Params.(Classifiers{c}).nTrees,length(Params.(Classifiers{c}).d)*length(Params.(Classifiers{c}).lambda));
            NumNodes.(Classifiers{c}) = NaN(Params.(Classifiers{c}).nTrees,length(Params.(Classifiers{c}).d)*length(Params.(Classifiers{c}).lambda));
            NumSplitNodes.(Classifiers{c}) = NaN(Params.(Classifiers{c}).nTrees,length(Params.(Classifiers{c}).d)*length(Params.(Classifiers{c}).lambda));
            TreeStrength.(Classifiers{c}) = NaN(1,length(Params.(Classifiers{c}).d)*length(Params.(Classifiers{c}).lambda));
            TreeDiversity.(Classifiers{c}) = NaN(1,length(Params.(Classifiers{c}).d)*length(Params.(Classifiers{c}).lambda));
            TestScores.(Classifiers{c}) = NaN(ntest,nClasses,length(Params.(Classifiers{c}).d)*length(Params.(Classifiers{c}).lambda));
            TestError.(Classifiers{c}) = NaN(1,length(Params.(Classifiers{c}).d)*length(Params.(Classifiers{c}).lambda));
        else
            OOBError.(Classifiers{c}) = NaN(1,length(Params.(Classifiers{c}).d)*length(Params.(Classifiers{c}).rho));
            OOBAUC.(Classifiers{c}) = NaN(1,length(Params.(Classifiers{c}).d)*length(Params.(Classifiers{c}).rho));
            TrainTime.(Classifiers{c}) = NaN(1,length(Params.(Classifiers{c}).d)*length(Params.(Classifiers{c}).rho));
            Depth.(Classifiers{c}) = NaN(Params.(Classifiers{c}).nTrees,length(Params.(Classifiers{c}).d)*length(Params.(Classifiers{c}).rho));
            NumNodes.(Classifiers{c}) = NaN(Params.(Classifiers{c}).nTrees,length(Params.(Classifiers{c}).d)*length(Params.(Classifiers{c}).rho));
            NumSplitNodes.(Classifiers{c}) = NaN(Params.(Classifiers{c}).nTrees,length(Params.(Classifiers{c}).d)*length(Params.(Classifiers{c}).rho));
            TreeStrength.(Classifiers{c}) = NaN(1,length(Params.(Classifiers{c}).d)*length(Params.(Classifiers{c}).rho));
            TreeDiversity.(Classifiers{c}) = NaN(1,length(Params.(Classifiers{c}).d)*length(Params.(Classifiers{c}).rho));
            TestScores.(Classifiers{c}) = NaN(ntest,nClasses,length(Params.(Classifiers{c}).d)*length(Params.(Classifiers{c}).rho));
            TestError.(Classifiers{c}) = NaN(1,length(Params.(Classifiers{c}).d)*length(Params.(Classifiers{c}).rho));
        end
    end

    % train classifier
    poolobj = gcp('nocreate');
    if isempty(poolobj)
        parpool('local',Params.(Classifiers{c}).NWorkers,...
            'IdleTimeout',360);
    end

    [Forest,~,TrainTime.(Classifiers{c})] = ...
        RerF_train(Xtrain,Ytrain,Params.(Classifiers{c}));

    fprintf('Training complete\n')

    % compute oob auc, oob error, and tree stats

    for k = 1:length(Forest)
        Labels = Forest{k}.classname;
        nClasses = length(Labels);
        Scores = rerf_oob_classprob(Forest{k},...
            Xtrain,'last');
        Predictions = predict_class(Scores,Labels);
        OOBError.(Classifiers{c})(k) = ...
            misclassification_rate(Predictions,Ytrain,...
        false);        
        if nClasses > 2
            Yb = binarize_labels(Ytrain,Labels);
            [~,~,~,OOBAUC.(Classifiers{c})(k)] = ... 
                perfcurve(Yb(:),Scores(:),'1');
        else
            [~,~,~,OOBAUC.(Classifiers{c})(k)] = ...
                perfcurve(Ytrain,Scores(:,2),'1');
        end
        Depth.(Classifiers{c})(:,k) = forest_depth(Forest{k})';
        NN = NaN(Forest{k}.nTrees,1);
        NS = NaN(Forest{k}.nTrees,1);
        Trees = Forest{k}.Tree;
        parfor kk = 1:Forest{k}.nTrees
            NN(kk) = Trees{kk}.numnodes;
            NS(kk) = sum(Trees{kk}.isbranch);
        end
        NumNodes.(Classifiers{c})(:,k) = NN;
        NumSplitNodes.(Classifiers{c})(:,k) = NS;

        if ~strcmp(Forest{k}.Rescale,'off')
            Scores = rerf_classprob(Forest{k},Xtest,'individual',true,Xtrain);
        else
            Scores = rerf_classprob(Forest{k},Xtest,'individual',true);
        end
        PredCell = cell(ntest,Params.(Classifiers{c}).nTrees);
        parfor kk = 1:Params.(Classifiers{c}).nTrees
            PredCell(:,kk) = predict_class(Scores(:,:,kk),Labels);
        end

        TreeStrength.(Classifiers{c})(k) = 1 - misclassification_rate(PredCell,Ytest,true);
        TreeDiversity.(Classifiers{c})(k) = classifier_variance(PredCell);
        
        if ~strcmp(Forest{k}.Rescale,'off')
            TestScores.(Classifiers{c})(:,:,k) = rerf_classprob(Forest{k},Xtest,'last',true,Xtrain);
        else
        TestScores.(Classifiers{c})(:,:,k) = rerf_classprob(Forest{k},Xtest,true,'last');
        end
        TestPredictions = predict_class(TestScores.(Classifiers{c})(:,:,k),Forest{k}.classname);

        TestError.(Classifiers{c})(k) = ...
            misclassification_rate(TestPredictions,Ytest,false);
    end

    % select best model based on OOB errors and AUCs
    BI = hp_optimize(OOBError.(Classifiers{c})(end,:),...
        OOBAUC.(Classifiers{c})(end,:));
    BestIdx.(Classifiers{c}) = BI(randperm(length(BI),1));

    tic;
    [importance,features] = feature_importance(Forest{BestIdx.(Classifiers{c})},'gini');
    Time.importance = toc;
    
    save(OutFile,'Params','OOBError','OOBAUC','TestError',...
        'TrainTime','Depth','NumNodes','NumSplitNodes','TreeStrength',...
        'TreeDiversity','BestIdx','TestScores','importance','features','Time')
end  