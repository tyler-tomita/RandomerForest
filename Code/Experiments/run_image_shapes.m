close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

rng(1);

load image_shapes_data

[ih,iw,ntrain] = size(Xtrain_image);
[~,~,ntest] = size(Xtest_image);
p = ih*iw;
Xtrain = reshape(Xtrain_image,p,ntrain)';
Xtest = reshape(Xtest_image,p,ntest)';
Ytrain = cellstr(num2str(Ytrain));
Ytest = cellstr(num2str(Ytest));

clear Xtrain_image Xtest_image

load Random_matrix_adjustment_factor

ntrees = 2000;
Stratified = true;
NWorkers = 24;

FileID = fopen('~/shapes.out','w');

for k = 3:length(ns)
        nTrain = ns(k);
        fprintf(FileID,'\nn = %d\n\n',nTrain);
        
        BaggedError.srerf{k} = NaN(ntrials,5);
        BaggedError.control{k} = NaN(ntrials,5);
        BaggedError.rerf{k} = NaN(ntrials,5);
        BaggedError.rf{k} = NaN(ntrials,4);
        AUC.srerf{k} = NaN(ntrials,5);
        AUC.control{k} = NaN(ntrials,5);
        AUC.rerf{k} = NaN(ntrials,5);
        AUC.rf{k} = NaN(ntrials,4);
    
    for trial = 1:ntrials
        fprintf(FileID,'trial = %d\n',trial);

        %% Structured RerF %%

        fprintf(FileID,'Structured RerF\n');

        ds = [ceil(p.^[1/4 1/2 3/4 1]) 10*p];
        for j = 1:length(ds)
            d = ds(j);
            
            fprintf(FileID,'d = %d\n',d);

            srerf{j} = rpclassificationforest(ntrees,...
                Xtrain(TrainIdx{k}(trial,:),:),Ytrain(TrainIdx{k}(trial,:)),...
                'Image','on','ih',ih,'iw',iw,'nvartosample',d,'NWorkers',...
                NWorkers,'Stratified',Stratified);
            Scores = rerf_oob_classprob(srerf{j},Xtrain(TrainIdx{k}(trial,:),:),'last');
            Predictions = predict_class(Scores,srerf{j}.classname);
            BaggedError.srerf{k}(trial,j) = misclassification_rate(Predictions,...
                Ytrain(TrainIdx{k}(trial,:)),false);
            [~,~,~,AUC.srerf{k}(trial,j)] = perfcurve(Ytrain(TrainIdx{k}(trial,:)),Scores(:,2),'1');
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
        fprintf(FileID,'RerF control\n');

        ds = [ceil(p.^[1/4 1/2 3/4 1]) 10*p];

        for j = 1:length(ds)
            d = ds(j);

            fprintf(FileID,'d = %d\n',d);

            control{j} = rpclassificationforest(ntrees,Xtrain(TrainIdx{k}(trial,:),:),Ytrain(TrainIdx{k}(trial,:)),...
                'Image','control','ih',ih,'iw',iw,'nvartosample',d,'NWorkers',...
                NWorkers,'Stratified',Stratified);
            Scores = rerf_oob_classprob(control{j},Xtrain(TrainIdx{k}(trial,:),:),'last');
            Predictions = predict_class(Scores,control{j}.classname);
            BaggedError.control{k}(trial,j) = misclassification_rate(Predictions,...
                Ytrain(TrainIdx{k}(trial,:)),false);
            [~,~,~,AUC.control{k}(trial,j)] = perfcurve(Ytrain(TrainIdx{k}(trial,:)),Scores(:,2),'1');
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
        fprintf(FileID,'RerF\n');

        ds = [ceil(p.^[1/4 1/2 3/4 1]) 10*p];

        for j = 1:length(ds)
            d = ds(j);

            fprintf(FileID,'d = %d\n',d);

            dprime = ceil(d^(1/interp1(ps,slope,p,'linear','extrap')));

            rerf{j} = rpclassificationforest(ntrees,Xtrain(TrainIdx{k}(trial,:),:),Ytrain(TrainIdx{k}(trial,:)),'sparsemethod',...
                'sparse-adjusted','nvartosample',d,'dprime',dprime,'NWorkers',NWorkers,...
                'Stratified',Stratified);
            Scores = rerf_oob_classprob(rerf{j},Xtrain(TrainIdx{k}(trial,:),:),'last');
            Predictions = predict_class(Scores,rerf{j}.classname);
            BaggedError.rerf{k}(trial,j) = misclassification_rate(Predictions,...
                Ytrain(TrainIdx{k}(trial,:)),false);
            [~,~,~,AUC.rerf{k}(trial,j)] = perfcurve(Ytrain(TrainIdx{k}(trial,:)),Scores(:,2),'1');
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
        fprintf(FileID,'Random Forest\n');

        ds = ceil(p.^[1/4 1/2 3/4 1]);

        for j = 1:length(ds)
            d = ds(j);

            fprintf(FileID,'d = %d\n',d);

            rf{j} = rpclassificationforest(ntrees,Xtrain(TrainIdx{k}(trial,:),:),Ytrain(TrainIdx{k}(trial,:)),'RandomForest',true,...
                'nvartosample',d,'NWorkers',NWorkers,'Stratified',Stratified);
            Scores = rerf_oob_classprob(rf{j},Xtrain(TrainIdx{k}(trial,:),:),'last');
            Predictions = predict_class(Scores,rf{j}.classname);
            BaggedError.rf{k}(trial,j) = misclassification_rate(Predictions,...
                Ytrain(TrainIdx{k}(trial,:)),false);
            [~,~,~,AUC.rf{k}(trial,j)] = perfcurve(Ytrain(TrainIdx{k}(trial,:)),Scores(:,2),'1');
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
        
        save([rerfPath 'RandomerForest/Results/image_shapes_n200.mat'],...
            'ntrees','BaggedError','AUC','TestError')
    end
end