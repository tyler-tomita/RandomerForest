% train and test classifiers on arrhythmia dataset

TrainFile = '/scratch/groups/jvogels3/tyler/Benchmarks/Data/dat/Corrupted/arrhythmia_train.dat';
TestFile = '/scratch/groups/jvogels3/tyler/Benchmarks/Data/dat/Corrupted/arrhythmia_test.dat';
OutFile = '/scratch/groups/jvogels3/tyler/RandomerForest/Results/2017.04.01/Benchmarks/Corrupted/arrhythmia.mat';

Classifiers = {'rf','rerf','rerfr','rr_rf','rr_rfr'};

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
    fprintf('%s start\n',Classifiers{c})

    if ntrain <= 2000
        Params.(Classifiers{c}).nTrees = 1000;
    else
        Params.(Classifiers{c}).nTrees = 500;
    end
    Params.(Classifiers{c}).Stratified = true;
    Params.(Classifiers{c}).NWorkers = 12;
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
        Params.(Classifiers{c}).rho = (1:min(p,3))/p;
    elseif strcmp(Classifiers{c},'rr_rf') || strcmp(Classifiers{c},'rr_rfr')
        Params.(Classifiers{c}).ForestMethod = 'rf';   
        Params.(Classifiers{c}).Rotate = true;
        Params.(Classifiers{c}).d = mtrys_rf;
    end
    
    if strcmp(Params.(Classifiers{c}).ForestMethod,'rf')
        OOBError.(Classifiers{c}) = NaN(Params.(Classifiers{c}).nTrees,length(Params.(Classifiers{c}).d));
        OOBAUC.(Classifiers{c}) = NaN(Params.(Classifiers{c}).nTrees,length(Params.(Classifiers{c}).d));
        TrainTime.(Classifiers{c}) = NaN(1,length(Params.(Classifiers{c}).d));
        Depth.(Classifiers{c}) = NaN(Params.(Classifiers{c}).nTrees,length(Params.(Classifiers{c}).d));
        NumNodes.(Classifiers{c}) = NaN(Params.(Classifiers{c}).nTrees,length(Params.(Classifiers{c}).d));
        NumSplitNodes.(Classifiers{c}) = NaN(Params.(Classifiers{c}).nTrees,length(Params.(Classifiers{c}).d));
        TreeStrength.(Classifiers{c}) = NaN(1,length(Params.(Classifiers{c}).d));
        TreeDiversity.(Classifiers{c}) = NaN(1,length(Params.(Classifiers{c}).d));
        TestScores.(Classifiers{c}) = NaN(ntest,nClasses,length(Params.(Classifiers{c}).d));
        TestError.(Classifiers{c}) = NaN(1,length(Params.(Classifiers{c}).d));
    else
        OOBError.(Classifiers{c}) = NaN(Params.(Classifiers{c}).nTrees,length(Params.(Classifiers{c}).d)*length(Params.(Classifiers{c}).rho));
        OOBAUC.(Classifiers{c}) = NaN(Params.(Classifiers{c}).nTrees,length(Params.(Classifiers{c}).d)*length(Params.(Classifiers{c}).rho));
        TrainTime.(Classifiers{c}) = NaN(1,length(Params.(Classifiers{c}).d)*length(Params.(Classifiers{c}).rho));
        Depth.(Classifiers{c}) = NaN(Params.(Classifiers{c}).nTrees,length(Params.(Classifiers{c}).d)*length(Params.(Classifiers{c}).rho));
        NumNodes.(Classifiers{c}) = NaN(Params.(Classifiers{c}).nTrees,length(Params.(Classifiers{c}).d)*length(Params.(Classifiers{c}).rho));
        NumSplitNodes.(Classifiers{c}) = NaN(Params.(Classifiers{c}).nTrees,length(Params.(Classifiers{c}).d)*length(Params.(Classifiers{c}).rho));
        TreeStrength.(Classifiers{c}) = NaN(1,length(Params.(Classifiers{c}).d)*length(Params.(Classifiers{c}).rho));
        TreeDiversity.(Classifiers{c}) = NaN(1,length(Params.(Classifiers{c}).d)*length(Params.(Classifiers{c}).rho));
        TestScores.(Classifiers{c}) = NaN(ntest,nClasses,length(Params.(Classifiers{c}).d)*length(Params.(Classifiers{c}).rho));
        TestError.(Classifiers{c}) = NaN(1,length(Params.(Classifiers{c}).d)*length(Params.(Classifiers{c}).rho));
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
            Xtrain,'every');
        OE = zeros(Forest{k}.nTrees,1);
        OA = zeros(Forest{k}.nTrees,1);
        if nClasses > 2
            parfor t = 1:Forest{k}.nTrees
                Scores_t = Scores(:,:,t);
                Predictions = predict_class(Scores_t,Labels);
                OE(t) = ...
                    misclassification_rate(Predictions,Ytrain,...
                false);
                Yb = binarize_labels(Ytrain,Labels);
                [~,~,~,OA(t)] = ... 
                    perfcurve(Yb(:),Scores_t(:),'1');
            end
        else
            PositiveScores = Scores(:,2,:);
            parfor t = 1:Forest{k}.nTrees
                Predictions = predict_class(Scores(:,:,t),Labels);
                OE(t) = ...
                    misclassification_rate(Predictions,Ytrain,...
                false);
                [~,~,~,OA(t)] = ...
                    perfcurve(Ytrain,PositiveScores(:,1,t),'1');
            end
        end
        OOBError.(Classifiers{c})(:,k) = OE;
        OOBAUC.(Classifiers{c})(:,k) = OA;
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

        Scores = rerf_classprob(Forest{k},Xtest,'individual');
        PredCell = cell(ntest,Params.(Classifiers{c}).nTrees);
        parfor kk = 1:Params.(Classifiers{c}).nTrees
            PredCell(:,kk) = predict_class(Scores(:,:,kk),Labels);
        end

        TreeStrength.(Classifiers{c})(k) = 1 - misclassification_rate(PredCell,Ytest,true);
        TreeDiversity.(Classifiers{c})(k) = classifier_variance(PredCell);
        
        TestScores.(Classifiers{c})(:,:,k) = rerf_classprob(Forest{k},Xtest,'last');
        TestPredictions = predict_class(TestScores.(Classifiers{c})(:,:,k),Forest{k}.classname);

        TestError.(Classifiers{c})(k) = ...
            misclassification_rate(TestPredictions,Ytest,false);
    end

    % select best model based on OOB errors and AUCs
    BI = hp_optimize(OOBError.(Classifiers{c})(end,:),...
        OOBAUC.(Classifiers{c})(end,:));
    BestIdx.(Classifiers{c}) = BI(end);

    clear Forest

    fprintf('%s complete\n',Classifiers{c})

    save(OutFile,'Params','OOBError','OOBAUC','TestError',...
        'TrainTime','Depth','NumNodes','NumSplitNodes','TreeStrength',...
        'TreeDiversity','BestIdx','TestScores')
end  

delete(gcp);
