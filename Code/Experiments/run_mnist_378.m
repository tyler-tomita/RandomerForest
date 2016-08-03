close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

rng(1);

load MNIST_train
load MNIST_test
load MNIST_378_train_indices
load Random_matrix_adjustment_factor

[ntrain,p] = size(Xtrain);
ih = sqrt(p);
iw = ih;

Ytrain = cellstr(num2str(Ytrain));
Xtest = Xtest(Ytest==3|Ytest==7|Ytest==8,:);
Ytest = cellstr(num2str(Ytest(Ytest==3|Ytest==7|Ytest==8)));
Labels = unique(Ytest);

ntrees = 2000;
Stratified = true;
NWorkers = 2;

for k = 1:length(ns)
        nTrain = ns(k);
        fprintf('\nn = %d\n\n',nTrain);
        
        BaggedError.srerf{k} = NaN(ntrials,5);
        BaggedError.control{k} = NaN(ntrials,5);
        BaggedError.rerf{k} = NaN(ntrials,5);
        BaggedError.rf{k} = NaN(ntrials,4);
        AUC.srerf{k} = NaN(ntrials,5);
        AUC.control{k} = NaN(ntrials,5);
        AUC.rerf{k} = NaN(ntrials,5);
        AUC.rf{k} = NaN(ntrials,4);
    
    for trial = 1:ntrials
        fprintf('trial = %d\n',trial);

        %% Structured RerF %%

        fprintf('Structured RerF\n');

        ds = [ceil(p.^[1/4 1/2 3/4 1]) 10*p];
        for j = 1:length(ds)
            d = ds(j);
            
            fprintf('d = %d\n',d);

            srerf{j} = rpclassificationforest(ntrees,...
                Xtrain(TrainIdx{k}(trial,:),:),Ytrain(TrainIdx{k}(trial,:)),...
                'Image','on','ih',ih,'iw',iw,'nvartosample',d,'NWorkers',...
                NWorkers,'Stratified',Stratified);
            Scores = rerf_oob_classprob(srerf{j},Xtrain(TrainIdx{k}(trial,:),:),'last');
            Predictions = predict_class(Scores,srerf{j}.classname);
            BaggedError.srerf{k}(trial,j) = misclassification_rate(Predictions,...
                Ytrain(TrainIdx{k}(trial,:)),false);
            Yb = binarize_labels(Ytrain(TrainIdx{k}(trial,:)),srerf{j}.classname);
            [~,~,~,AUC.srerf{k}(trial,j)] = perfcurve(Yb(:),Scores(:),1);
        end
        
        BestIdx = hp_optimize(BaggedError.srerf{k}(trial,:),AUC.srerf{k}(trial,:));
        if length(BestIdx)>1
            BestIdx = BestIdx(end);
        end
        
        Scores = rerf_classprob(srerf{BestIdx},Xtest,'last');
        Predictions = predict_class(Scores,srerf{BestIdx}.classname);
        TestError.srerf{k}(trial) = misclassification_rate(Predictions,...
            Ytest,false);
        
        clear srerf
        
        %% RerF controlled for density%%
        fprintf('RerF control\n');

        ds = [ceil(p.^[1/4 1/2 3/4 1]) 10*p];

        for j = 1:length(ds)
            d = ds(j);

            fprintf('d = %d\n',d);

            control{j} = rpclassificationforest(ntrees,Xtrain(TrainIdx{k}(trial,:),:),Ytrain(TrainIdx{k}(trial,:)),...
                'Image','control','ih',ih,'iw',iw,'nvartosample',d,'NWorkers',...
                NWorkers,'Stratified',Stratified);
            Scores = rerf_oob_classprob(control{j},Xtrain(TrainIdx{k}(trial,:),:),'last');
            Predictions = predict_class(Scores,control{j}.classname);
            BaggedError.control{k}(trial,j) = misclassification_rate(Predictions,...
                Ytrain(TrainIdx{k}(trial,:)),false);
            Yb = binarize_labels(Ytrain(TrainIdx{k}(trial,:)),control{j}.classname);
            [~,~,~,AUC.control{k}(trial,j)] = perfcurve(Yb(:),Scores(:),1);
        end
        
        BestIdx = hp_optimize(BaggedError.control{k}(trial,:),AUC.control{k}(trial,:));
        if length(BestIdx)>1
            BestIdx = BestIdx(end);
        end
        
        Scores = rerf_classprob(control{BestIdx},Xtest,'last');
        Predictions = predict_class(Scores,control{BestIdx}.classname);
        TestError.control{k}(trial) = misclassification_rate(Predictions,...
            Ytest,false);
        
        clear control

        %% RerF %%
        fprintf('RerF\n');

        ds = [ceil(p.^[1/4 1/2 3/4 1]) 10*p];

        for j = 1:length(ds)
            d = ds(j);

            fprintf('d = %d\n',d);

            dprime = ceil(d^(1/interp1(ps,slope,p,'linear','extrap')));

            rerf{j} = rpclassificationforest(ntrees,Xtrain(TrainIdx{k}(trial,:),:),Ytrain(TrainIdx{k}(trial,:)),'sparsemethod',...
                'sparse-adjusted','nvartosample',d,'dprime',dprime,'NWorkers',NWorkers,...
                'Stratified',Stratified);
            Scores = rerf_oob_classprob(rerf{j},Xtrain(TrainIdx{k}(trial,:),:),'last');
            Predictions = predict_class(Scores,rerf{j}.classname);
            BaggedError.rerf{k}(trial,j) = misclassification_rate(Predictions,...
                Ytrain(TrainIdx{k}(trial,:)),false);
            Yb = binarize_labels(Ytrain(TrainIdx{k}(trial,:)),rerf{j}.classname);
            [~,~,~,AUC.rerf{k}(trial,j)] = perfcurve(Yb(:),Scores(:),1);
        end
        
        BestIdx = hp_optimize(BaggedError.rerf{k}(trial,:),AUC.rerf{k}(trial,:));
        if length(BestIdx)>1
            BestIdx = BestIdx(end);
        end
        
        Scores = rerf_classprob(rerf{BestIdx},Xtest,'last');
        Predictions = predict_class(Scores,rerf{BestIdx}.classname);
        TestError.rerf{k}(trial) = misclassification_rate(Predictions,...
            Ytest,false);
        
        clear rerf

        %% RF %%
        fprintf('Random Forest\n');

        ds = ceil(p.^[1/4 1/2 3/4 1]);

        for j = 1:length(ds)
            d = ds(j);

            fprintf('d = %d\n',d);

            rf{j} = rpclassificationforest(ntrees,Xtrain(TrainIdx{k}(trial,:),:),Ytrain(TrainIdx{k}(trial,:)),'RandomForest',true,...
                'nvartosample',d,'NWorkers',NWorkers,'Stratified',Stratified);
            Scores = rerf_oob_classprob(rf{j},Xtrain(TrainIdx{k}(trial,:),:),'last');
            Predictions = predict_class(Scores,rf{j}.classname);
            BaggedError.rf{k}(trial,j) = misclassification_rate(Predictions,...
                Ytrain(TrainIdx{k}(trial,:)),false);
            Yb = binarize_labels(Ytrain(TrainIdx{k}(trial,:)),rf{j}.classname);
            [~,~,~,AUC.rf{k}(trial,j)] = perfcurve(Yb(:),Scores(:),1);
        end
        
        BestIdx = hp_optimize(BaggedError.rf{k}(trial,:),AUC.rf{k}(trial,:));
        if length(BestIdx)>1
            BestIdx = BestIdx(end);
        end
        
        Scores = rerf_classprob(rf{BestIdx},Xtest,'last');
        Predictions = predict_class(Scores,rf{BestIdx}.classname);
        TestError.rf{k}(trial) = misclassification_rate(Predictions,...
            Ytest,false);
        
        clear rf
        
        save([rerfPath 'RandomerForest/Results/mnist_378.mat'],...
            'ntrees','BaggedError','AUC','TestError')
    end
end