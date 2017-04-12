close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

rng(1);

ps = [10,100,1000];
ns = {[10,100,1000,10000], [10,100,1000,10000], [10,100,1000,10000]};
ntrials = 10;
ntest = 10e3;

Classifiers = {'rr_rf'};

OOBError = cell(length(ns{1}),length(ps));
OOBAUC = cell(length(ns{1}),length(ps));
TrainTime = cell(length(ns{1}),length(ps));
Depth = cell(length(ns{1}),length(ps));
NumNodes = cell(length(ns{1}),length(ps));
NumSplitNodes = cell(length(ns{1}),length(ps));
TestError = cell(length(ns{1}),length(ps));
Bias = cell(length(ns{1}),length(ps));
Variance = cell(length(ns{1}),length(ps));
MR = cell(length(ns{1}),length(ps));
BestIdx = cell(length(ns{1}),length(ps));
Noise = zeros(1,length(ps));
Labels = {'0';'1'};

for j = 1:length(ps)
    p = ps(j);
    fprintf('p = %d\n',p)
    
    Xtest = dlmread(sprintf('/scratch/groups/jvogels3/tyler/R/Data/Trunk/dat/Scaled/Test/Trunk_test_set_p%d.dat',p));
    Ytest = cellstr(num2str(Xtest(:,end)));
    Xtest(:,end) = [];
    ClassPosteriors = dlmread(sprintf('/scratch/groups/jvogels3/tyler/R/Data/Trunk/dat/Scaled/Test/Trunk_test_set_posteriors_p%d.dat',p));
    Noise(j) = 0.5*(1 - mean(sum(ClassPosteriors.^2,2)));   % noise from bias-variance-noise decomposition for 0-1 loss
    
    if p <= 10
        mtrys = ceil(p.^[1/4 1/2 3/4 1 2]);
    else
        mtrys = [ceil(p.^[1/4 1/2 3/4 1]) 20*p];
    end
    mtrys_rf = mtrys(mtrys<=p);
      
    for i = 1:length(ns{j})
        fprintf('n = %d\n',ns{j}(i))

        for c = 1:length(Classifiers)
            fprintf('%s start\n',Classifiers{c})
            
            if ns{j}(i) <= 1000
                Params{i,j}.(Classifiers{c}).nTrees = 1000;
            else
                Params{i,j}.(Classifiers{c}).nTrees = 500;
            end
            Params{i,j}.(Classifiers{c}).Stratified = true;
            Params{i,j}.(Classifiers{c}).NWorkers = 12;
            if strcmp(Classifiers{c},'rerfr') || strcmp(Classifiers{c},'frcr') || strcmp(Classifiers{c},'rr_rfr')
                Params{i,j}.(Classifiers{c}).Rescale = 'rank';
            else
                Params{i,j}.(Classifiers{c}).Rescale = 'off';
            end
            Params{i,j}.(Classifiers{c}).mdiff = 'off';
            if strcmp(Classifiers{c},'rf')
                Params{i,j}.(Classifiers{c}).ForestMethod = 'rf';
                Params{i,j}.(Classifiers{c}).d = mtrys_rf;
            elseif strcmp(Classifiers{c},'rerf') || strcmp(Classifiers{c},'rerfr')
                Params{i,j}.(Classifiers{c}).ForestMethod = 'rerf';
                Params{i,j}.(Classifiers{c}).RandomMatrix = 'binary';
                Params{i,j}.(Classifiers{c}).d = mtrys;
                Params{i,j}.(Classifiers{c}).rho = (1:min(p,3))/p;
            elseif strcmp(Classifiers{c},'frc') || strcmp(Classifiers{c},'frcr')
                Params{i,j}.(Classifiers{c}).ForestMethod = 'rerf';
                Params{i,j}.(Classifiers{c}).RandomMatrix = 'frc';
                Params{i,j}.(Classifiers{c}).d = mtrys;
                Params{i,j}.(Classifiers{c}).L = 2:4;
            elseif strcmp(Classifiers{c},'rr_rf') || strcmp(Classifiers{c},'rr_rfr')
                Params{i,j}.(Classifiers{c}).ForestMethod = 'rf';   
                Params{i,j}.(Classifiers{c}).Rotate = true;
                Params{i,j}.(Classifiers{c}).d = mtrys_rf;
            end

            if strcmp(Params{i,j}.(Classifiers{c}).ForestMethod,'rerf')
                if strcmp(Params{i,j}.(Classifiers{c}).RandomMatrix,'frc')
                    OOBError{i,j}.(Classifiers{c}) = NaN(ntrials,length(Params{i,j}.(Classifiers{c}).d)*length(Params{i,j}.(Classifiers{c}).L),Params{i,j}.(Classifiers{c}).nTrees);
                    OOBAUC{i,j}.(Classifiers{c}) = NaN(ntrials,length(Params{i,j}.(Classifiers{c}).d)*length(Params{i,j}.(Classifiers{c}).L),Params{i,j}.(Classifiers{c}).nTrees);
                    TrainTime{i,j}.(Classifiers{c}) = NaN(ntrials,length(Params{i,j}.(Classifiers{c}).d)*length(Params{i,j}.(Classifiers{c}).L));
                    Depth{i,j}.(Classifiers{c}) = NaN(ntrials,Params{i,j}.(Classifiers{c}).nTrees,length(Params{i,j}.(Classifiers{c}).d)*length(Params{i,j}.(Classifiers{c}).L));
                    NumNodes{i,j}.(Classifiers{c}) = NaN(ntrials,Params{i,j}.(Classifiers{c}).nTrees,length(Params{i,j}.(Classifiers{c}).d)*length(Params{i,j}.(Classifiers{c}).L));
                    NumSplitNodes{i,j}.(Classifiers{c}) = NaN(ntrials,Params{i,j}.(Classifiers{c}).nTrees,length(Params{i,j}.(Classifiers{c}).d)*length(Params{i,j}.(Classifiers{c}).L));
                    TreeStrength.(Classifiers{c}) = NaN(ntrials,length(Params.(Classifiers{c}).d)*length(Params.(Classifiers{c}).L));
                    TreeDiversity.(Classifiers{c}) = NaN(ntrials,length(Params.(Classifiers{c}).d)*length(Params.(Classifiers{c}).L));
                    Bias{i,j}.(Classifiers{c}) = NaN(1,length(Params{i,j}.(Classifiers{c}).d)*length(Params{i,j}.(Classifiers{c}).L));
                    Variance{i,j}.(Classifiers{c}) = NaN(1,length(Params{i,j}.(Classifiers{c}).d)*length(Params{i,j}.(Classifiers{c}).L));
                    MR{i,j}.(Classifiers{c}) = NaN(1,length(Params{i,j}.(Classifiers{c}).d)*length(Params{i,j}.(Classifiers{c}).L));
                    TestPredictions = cell(ntest,ntrials,length(Params{i,j}.(Classifiers{c}).d)*length(Params{i,j}.(Classifiers{c}).L));
                else
                    OOBError{i,j}.(Classifiers{c}) = NaN(ntrials,length(Params{i,j}.(Classifiers{c}).d)*length(Params{i,j}.(Classifiers{c}).rho),Params{i,j}.(Classifiers{c}).nTrees);
                    OOBAUC{i,j}.(Classifiers{c}) = NaN(ntrials,length(Params{i,j}.(Classifiers{c}).d)*length(Params{i,j}.(Classifiers{c}).rho),Params{i,j}.(Classifiers{c}).nTrees);
                    TrainTime{i,j}.(Classifiers{c}) = NaN(ntrials,length(Params{i,j}.(Classifiers{c}).d)*length(Params{i,j}.(Classifiers{c}).rho));
                    Depth{i,j}.(Classifiers{c}) = NaN(ntrials,Params{i,j}.(Classifiers{c}).nTrees,length(Params{i,j}.(Classifiers{c}).d)*length(Params{i,j}.(Classifiers{c}).rho));
                    NumNodes{i,j}.(Classifiers{c}) = NaN(ntrials,Params{i,j}.(Classifiers{c}).nTrees,length(Params{i,j}.(Classifiers{c}).d)*length(Params{i,j}.(Classifiers{c}).rho));
                    NumSplitNodes{i,j}.(Classifiers{c}) = NaN(ntrials,Params{i,j}.(Classifiers{c}).nTrees,length(Params{i,j}.(Classifiers{c}).d)*length(Params{i,j}.(Classifiers{c}).rho));
                    TreeStrength.(Classifiers{c}) = NaN(ntrials,length(Params.(Classifiers{c}).d)*length(Params.(Classifiers{c}).rho));
                    TreeDiversity.(Classifiers{c}) = NaN(ntrials,length(Params.(Classifiers{c}).d)*length(Params.(Classifiers{c}).rho));                    
                    Bias{i,j}.(Classifiers{c}) = NaN(1,length(Params{i,j}.(Classifiers{c}).d)*length(Params{i,j}.(Classifiers{c}).rho));
                    Variance{i,j}.(Classifiers{c}) = NaN(1,length(Params{i,j}.(Classifiers{c}).d)*length(Params{i,j}.(Classifiers{c}).rho));
                    MR{i,j}.(Classifiers{c}) = NaN(1,length(Params{i,j}.(Classifiers{c}).d)*length(Params{i,j}.(Classifiers{c}).rho));
                    TestPredictions = cell(ntest,ntrials,length(Params{i,j}.(Classifiers{c}).d)*length(Params{i,j}.(Classifiers{c}).rho));
                end
            else
                OOBError{i,j}.(Classifiers{c}) = NaN(ntrials,length(Params{i,j}.(Classifiers{c}).d),Params{i,j}.(Classifiers{c}).nTrees);
                OOBAUC{i,j}.(Classifiers{c}) = NaN(ntrials,length(Params{i,j}.(Classifiers{c}).d),Params{i,j}.(Classifiers{c}).nTrees);
                TrainTime{i,j}.(Classifiers{c}) = NaN(ntrials,length(Params{i,j}.(Classifiers{c}).d));
                Depth{i,j}.(Classifiers{c}) = NaN(ntrials,Params{i,j}.(Classifiers{c}).nTrees,length(Params{i,j}.(Classifiers{c}).d));
                NumNodes{i,j}.(Classifiers{c}) = NaN(ntrials,Params{i,j}.(Classifiers{c}).nTrees,length(Params{i,j}.(Classifiers{c}).d));
                NumSplitNodes{i,j}.(Classifiers{c}) = NaN(ntrials,Params{i,j}.(Classifiers{c}).nTrees,length(Params{i,j}.(Classifiers{c}).d));
                TreeStrength.(Classifiers{c}) = NaN(ntrials,length(Params.(Classifiers{c}).d));
                TreeDiversity.(Classifiers{c}) = NaN(ntrials,length(Params.(Classifiers{c}).d));                
                Bias{i,j}.(Classifiers{c}) = NaN(1,length(Params{i,j}.(Classifiers{c}).d));
                Variance{i,j}.(Classifiers{c}) = NaN(1,length(Params{i,j}.(Classifiers{c}).d));
                MR{i,j}.(Classifiers{c}) = NaN(1,length(Params{i,j}.(Classifiers{c}).d));
                TestPredictions = cell(ntest,ntrials,length(Params{i,j}.(Classifiers{c}).d));
            end
            TestError{i,j}.(Classifiers{c}) = NaN(ntrials,1);                
            BestIdx{i,j}.(Classifiers{c}) = NaN(ntrials,1);

            for trial = 1:ntrials
                fprintf('Trial %d\n',trial)
                
                Xtrain = dlmread(sprintf('/scratch/groups/jvogels3/tyler/R/Data/Trunk/dat/Scaled/Train/Trunk_train_set_n%d_p%d_trial%d.dat',ns{j}(i),p,trial));
                Ytrain = cellstr(num2str(Xtrain(:,end)));
                Xtrain(:,end) = [];

                % train classifier
                poolobj = gcp('nocreate');
                if isempty(poolobj)
                    parpool('local',Params{i,j}.(Classifiers{c}).NWorkers,...
                        'IdleTimeout',360);
                end

                [Forest,~,TrainTime{i,j}.(Classifiers{c})(trial,:)] = ...
                    RerF_train(Xtrain,Ytrain,Params{i,j}.(Classifiers{c}));
                
                fprintf('Training complete\n')

                % compute oob auc, oob error, and tree stats

                for k = 1:length(Forest)
                    Scores = rerf_oob_classprob(Forest{k},...
                        Xtrain,'every');
                    for t = 1:Forest{k}.nTrees
                        Predictions = predict_class(Scores(:,:,t),Forest{k}.classname);
                        OOBError{i,j}.(Classifiers{c})(trial,k,t) = ...
                            misclassification_rate(Predictions,Ytrain,...
                        false);
                        if size(Scores,2) > 2
                            Yb = binarize_labels(Ytrain,Forest{k}.classname);
                            [~,~,~,OOBAUC{i,j}.(Classifiers{c})(trial,k,t)] = ... 
                                perfcurve(Yb(:),Scores((t-1)*ns{j}(i)*Params{i,j}.(Classifiers{c}).d+(1:ns{j}(i)*Params{i,j}.(Classifiers{c}).d)),'1');
                        else
                            [~,~,~,OOBAUC{i,j}.(Classifiers{c})(trial,k,t)] = ...
                                perfcurve(Ytrain,Scores(:,2,t),'1');
                        end
                    end
                    Depth{i,j}.(Classifiers{c})(trial,:,k) = forest_depth(Forest{k})';
                    NN = NaN(1,Forest{k}.nTrees);
                    NS = NaN(1,Forest{k}.nTrees);
                    Trees = Forest{k}.Tree;
                    parfor kk = 1:Forest{k}.nTrees
                        NN(kk) = Trees{kk}.numnodes;
                        NS(kk) = sum(Trees{kk}.isbranch);
                    end
                    NumNodes{i,j}.(Classifiers{c})(trial,:,k) = NN;
                    NumSplitNodes{i,j}.(Classifiers{c})(trial,:,k) = NS;
                    
                    if ~strcmp(Forest{k}.Rescale,'off')
                        Scores = rerf_classprob(Forest{k},Xtest,'individual',Xtrain);
                    else
                        Scores = rerf_classprob(Forest{k},Xtest,'individual');
                    end
                    PredCell = cell(ntest,Params{i,j}.(Classifiers{c}).nTrees);
                    parfor kk = 1:Params{i,j}.(Classifiers{c}).nTrees
                        PredCell(:,kk) = predict_class(Scores(:,:,kk),Labels);
                    end

                    TreeStrength.(Classifiers{c})(trial,k) = 1 - misclassification_rate(PredCell,Ytest,true);
                    TreeDiversity.(Classifiers{c})(trial,k) = classifier_variance(PredCell);                    
                    
                    if ~strcmp(Forest{k}.Rescale,'off')
                        Scores = rerf_classprob(Forest{k},Xtest,'last',Xtrain);
                    else
                        Scores = rerf_classprob(Forest{k},Xtest,'last');
                    end
                    TestPredictions(:,trial,k) = predict_class(Scores,Forest{k}.classname);
                end
                
                % select best model for test predictions
                BI = hp_optimize(OOBError{i,j}.(Classifiers{c})(trial,:,end),...
                    OOBAUC{i,j}.(Classifiers{c})(trial,:,end));
                BestIdx{i,j}.(Classifiers{c})(trial) = BI(end);
                
                TestError{i,j}.(Classifiers{c})(trial) = ...
                    misclassification_rate(TestPredictions(:,trial,BestIdx{i,j}.(Classifiers{c})(trial)),Ytest,false);

                clear Forest
            end
            for k = 1:length(Params{i,j}.(Classifiers{c}).d)
                Bias{i,j}.(Classifiers{c})(k) = classifier_bias(TestPredictions(:,:,k),ClassPosteriors);
                Variance{i,j}.(Classifiers{c})(k) = classifier_variance(TestPredictions(:,:,k));
            end
            fprintf('%s complete\n',Classifiers{c})
            save([rerfPath 'RandomerForest/Results/2017.04.12/Trunk/Trunk_scaled_vary_n_' Classifiers{c} '.mat'],'ps',...
                'ns','Params','OOBError','OOBAUC','TestError',...
                'TrainTime','Depth','NumNodes','NumSplitNodes','Bias',...
                'Variance','Noise','BestIdx','TreeStrength','TreeDiversity')
        end
    end   
end
