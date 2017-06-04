%% Evaluate on Gaussian oblique binary classification problem

close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

rng(1);

ps = 2;
ns = 10;
ntrials = 5;
% ns = [10,50,100,500,1000];
% ntrials = 50;
ntest = 10e3;

Classifiers = {'rf','rerf','frc','rr_rf'};

OOBError = cell(length(ns),length(ps));
OOBAUC = cell(length(ns),length(ps));
TrainTime = cell(length(ns),length(ps));
Depth = cell(length(ns),length(ps));
NumNodes = cell(length(ns),length(ps));
NumSplitNodes = cell(length(ns),length(ps));
TestError = cell(length(ns),length(ps));
Bias = cell(length(ns),length(ps));
Variance = cell(length(ns),length(ps));
MR = cell(length(ns),length(ps));
BestIdx = cell(length(ns),length(ps));
Noise = zeros(1,length(ps));
Labels = {'0';'1'};

for j = 1:length(ps)
    p = ps(j);
    fprintf('p = %d\n',p)
    
    Xtest = dlmread(sprintf('~/R/Data/Gaussian/dat/Test/Gaussian_45_deg_test_set_p%d.dat',p));
    Ytest = cellstr(num2str(Xtest(:,end)));
    Xtest(:,end) = [];
    ClassPosteriors = dlmread(sprintf('~/R/Data/Gaussian/dat/Test/Gaussian_test_set_posteriors_p%d.dat',p));
    Noise(j) = 1 - mean(max(ClassPosteriors,[],2));   % noise from bias-variance-noise decomposition for 0-1 loss
    
    if p <= 10
        mtrys = ceil(p.^[1/4 1/2 3/4 1 2 3]);
    else
        mtrys = [ceil(p.^[1/4 1/2 3/4 1]) 20*p];
    end
    mtrys_rf = mtrys(mtrys<=p);
      
    for i = 1:length(ns)
        fprintf('n = %d\n',ns(i))

        for c = 1:length(Classifiers)
            cl = Classifiers{c};
            fprintf('%s start\n',cl)
            
            if ns(i) <= 1000
                Params{i,j}.(cl).nTrees = 1000;
            else
                Params{i,j}.(cl).nTrees = 500;
            end
            Params{i,j}.(cl).Stratified = true;
            Params{i,j}.(cl).NWorkers = 16;
            if strcmp(cl,'rerfr') || strcmp(cl,'frcr') || strcmp(cl,'rr_rfr')
                Params{i,j}.(cl).Rescale = 'rank';
            else
                Params{i,j}.(cl).Rescale = 'off';
            end
            Params{i,j}.(cl).mdiff = 'off';
            if strcmp(cl,'rf')
                Params{i,j}.(cl).ForestMethod = 'rf';
                Params{i,j}.(cl).d = mtrys_rf;
            elseif strcmp(cl,'rerf') || strcmp(cl,'rerfr')
                Params{i,j}.(cl).ForestMethod = 'rerf';
                Params{i,j}.(cl).RandomMatrix = 'binary';
                Params{i,j}.(cl).d = mtrys;
                Params{i,j}.(cl).rho = 1/p;
            elseif strcmp(cl,'frc') || strcmp(cl,'frcr')
                Params{i,j}.(cl).ForestMethod = 'rerf';
                Params{i,j}.(cl).RandomMatrix = 'frc';
                Params{i,j}.(cl).d = mtrys;
                Params{i,j}.(cl).L = 2;
            elseif strcmp(cl,'rr_rf') || strcmp(cl,'rr_rfr')
                Params{i,j}.(cl).ForestMethod = 'rf';   
                Params{i,j}.(cl).Rotate = true;
                Params{i,j}.(cl).d = mtrys_rf;
            end

            if strcmp(Params{i,j}.(cl).ForestMethod,'rerf')
                if strcmp(Params{i,j}.(cl).RandomMatrix,'frc')
                    OOBError{i,j}.(cl) = NaN(ntrials,length(Params{i,j}.(cl).d)*length(Params{i,j}.(cl).L));
                    OOBAUC{i,j}.(cl) = NaN(ntrials,length(Params{i,j}.(cl).d)*length(Params{i,j}.(cl).L));
                    TrainTime{i,j}.(cl) = NaN(ntrials,length(Params{i,j}.(cl).d)*length(Params{i,j}.(cl).L));
                    Bias{i,j}.(cl) = NaN(1,length(Params{i,j}.(cl).d)*length(Params{i,j}.(cl).L));
                    Variance{i,j}.(cl) = NaN(1,length(Params{i,j}.(cl).d)*length(Params{i,j}.(cl).L));
                    TestError{i,j}.(cl) = NaN(ntrials,length(Params{i,j}.(cl).d)*length(Params{i,j}.(cl).L));
                    TestPredictions = cell(ntest,ntrials,length(Params{i,j}.(cl).d)*length(Params{i,j}.(cl).L));
                else
                    OOBError{i,j}.(cl) = NaN(ntrials,length(Params{i,j}.(cl).d)*length(Params{i,j}.(cl).rho));
                    OOBAUC{i,j}.(cl) = NaN(ntrials,length(Params{i,j}.(cl).d)*length(Params{i,j}.(cl).rho));
                    TrainTime{i,j}.(cl) = NaN(ntrials,length(Params{i,j}.(cl).d)*length(Params{i,j}.(cl).rho));                  
                    Bias{i,j}.(cl) = NaN(1,length(Params{i,j}.(cl).d)*length(Params{i,j}.(cl).rho));
                    Variance{i,j}.(cl) = NaN(1,length(Params{i,j}.(cl).d)*length(Params{i,j}.(cl).rho));
                    TestError{i,j}.(cl) = NaN(ntrials,length(Params{i,j}.(cl).d)*length(Params{i,j}.(cl).rho));
                    TestPredictions = cell(ntest,ntrials,length(Params{i,j}.(cl).d)*length(Params{i,j}.(cl).rho));
                end
            else
                OOBError{i,j}.(cl) = NaN(ntrials,length(Params{i,j}.(cl).d));
                OOBAUC{i,j}.(cl) = NaN(ntrials,length(Params{i,j}.(cl).d));
                TrainTime{i,j}.(cl) = NaN(ntrials,length(Params{i,j}.(cl).d));               
                Bias{i,j}.(cl) = NaN(1,length(Params{i,j}.(cl).d));
                Variance{i,j}.(cl) = NaN(1,length(Params{i,j}.(cl).d));
                TestError{i,j}.(cl) = NaN(ntrials,length(Params{i,j}.(cl).d));
                TestPredictions = cell(ntest,ntrials,length(Params{i,j}.(cl).d));
            end                 
            BestIdx{i,j}.(cl) = NaN(ntrials,1);
            
            TP = cell(ntest,ntrials);
            for trial = 1:ntrials
                fprintf('Trial %d\n',trial)
                
                Xtrain = dlmread(sprintf('~/R/Data/Gaussian/dat/Train/Gaussian_45_deg_train_set_n%d_p%d_trial%d.dat',ns(i),p,trial));
                Ytrain = cellstr(num2str(Xtrain(:,end)));
                Xtrain(:,end) = [];

                % train classifier
                poolobj = gcp('nocreate');
                if isempty(poolobj)
                    parpool('local',Params{i,j}.(cl).NWorkers,...
                        'IdleTimeout',360);
                end

                [Forest,~,TrainTime{i,j}.(cl)(trial,:)] = ...
                    RerF_train(Xtrain,Ytrain,Params{i,j}.(cl));
                
                fprintf('Training complete\n')

                % compute oob auc, oob error, and tree stats

                for k = 1:length(Forest)
                    fprintf('Forest %d\n',k)
                    Labels = Forest{k}.classname;
                    nClasses = length(Labels);
                    Scores = rerf_oob_classprob(Forest{k},...
                        Xtrain,'last');
                    Predictions = predict_class(Scores,Labels);
                    OOBError{i,j}.(cl)(trial,k) = ...
                        misclassification_rate(Predictions,Ytrain,...
                    false);
                    if nClasses > 2
                        Yb = binarize_labels(Ytrain,Labels);
                        [~,~,~,OOBAUC{i,j}.(cl)(trial,k)] = ... 
                            perfcurve(Yb(:),Scores(:),'1');
                    else
                        [~,~,~,OOBAUC{i,j}.(cl)(trial,k)] = ...
                            perfcurve(Ytrain,Scores(:,2),'1');
                    end                                    
                    
                    if ~strcmp(Forest{k}.Rescale,'off')
                        Scores = rerf_classprob(Forest{k},Xtest,'last',Xtrain);
                    else
                        Scores = rerf_classprob(Forest{k},Xtest,'last');
                    end
                    TestPredictions(:,trial,k) = predict_class(Scores,Labels);
                    TestError{i,j}.(cl)(trial,k) = ...
                        misclassification_rate(TestPredictions(:,trial,k),Ytest,false);
                end
                
                % select best model for test predictions
                BI = hp_optimize(OOBError{i,j}.(cl)(trial,:),...
                    OOBAUC{i,j}.(cl)(trial,:));
                BestIdx{i,j}.(cl)(trial) = BI(end);
                
                TP(:,trial) = TestPredictions(:,trial,BestIdx{i,j}.(cl)(trial));

                clear Forest
            end
            [Bias{i,j}.(cl),Variance{i,j}.(cl)] = bias_variance(TP,ClassPosteriors);
            fprintf('%s complete\n',cl)
            save([rerfPath 'RandomerForest/Results/2017.06.04/Gaussian/Gaussian_45_deg_2017_06_04.mat'],'ps',...
                'ns','Params','OOBError','OOBAUC','TestError',...
                'TrainTime','Bias','Variance','BestIdx')
        end
    end   
end