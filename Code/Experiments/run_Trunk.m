close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

rng(1);

load Trunk_data
load Random_matrix_adjustment_factor

ntrees = 1000;
NWorkers = 16;
Class = [0;1];
BaggedError.rf = NaN(length(dims),5,ntrials);
BaggedError.rerf = NaN(length(dims),7,ntrials);
BaggedError.rerfr = NaN(length(dims),7,ntrials);
BaggedError.rerfdn = NaN(length(dims),7,ntrials);
BaggedError.rf_rot = NaN(length(dims),5,ntrials);
BaggedError.frc = NaN(length(dims),7,ntrials);
AUC.rf = NaN(length(dims),5,ntrials);
AUC.rerf = NaN(length(dims),7,ntrials);
AUC.rerfr = NaN(length(dims),7,ntrials);
AUC.rerfdn = NaN(length(dims),7,ntrials);
AUC.rf_rot = NaN(length(dims),5,ntrials);
AUC.frc = NaN(length(dims),7,ntrials);
TestError.rf = NaN(ntrials,length(dims));
TestError.rerf = NaN(ntrials,length(dims));
TestError.rerfr = NaN(ntrials,length(dims));
TestError.rerfdn = NaN(ntrials,length(dims));
TestError.rf_rot = NaN(ntrials,length(dims));
TestError.frc = NaN(ntrials,length(dims));
trainTime.rf = NaN(length(dims),5,ntrials);
trainTime.rerf = NaN(length(dims),7,ntrials);
trainTime.rerfr = NaN(length(dims),7,ntrials);
trainTime.rerfdn = NaN(length(dims),7,ntrials);
trainTime.rf_rot = NaN(length(dims),5,ntrials);
trainTime.frc = NaN(length(dims),7,ntrials);

for j = 1:length(dims)
    
    d = dims(j);
    
    if d <= 5
        mtrys = [1:d ceil(d.^[1.5 2])];
    elseif d > 5 && d <= 100
        mtrys = ceil(d.^[0 1/4 1/2 3/4 1 1.5 2]);
    else
        mtrys = [ceil(d.^[0 1/4 1/2 3/4 1]) 5*d 10*d];
    end
    
    nmix = 2;

    for trial = 1:ntrials

        fprintf('trial %d\n',trial)

        i = 1;

        poolobj = gcp('nocreate');
        if isempty(poolobj)
            parpool('local',NWorkers);
        end
        
        fprintf('RF\n')
        for i = 1:length(mtrys)
            
            mtry = mtrys(i);

            if mtry <= d
                fprintf('mtry = %d\n',mtry)

                tic;
                rf{i} = rpclassificationforest(ntrees,Xtrain{j}(:,:,trial),...
                    Ytrain{j}(:,trial),'RandomForest',true,'nvartosample',...
                    mtry,'NWorkers',NWorkers,'Stratified',true);
                trainTime.rf(j,i,trial) = toc;
            Scores = rerf_oob_classprob(rf{i},Xtrain{j}(:,:,trial),'last');
            Predictions = predict_class(Scores,rf{i}.classname);
            BaggedError.rf(j,i,trial) = ...
                misclassification_rate(Predictions,Ytrain{j}(:,trial),false);
            [~,~,~,AUC.rf(j,i,trial)] = ...
                perfcurve(Ytrain{j}(:,trial),Scores(:,2),'1');
            end          
        end
        
        BestIdx = hp_optimize(BaggedError.rf(j,:,trial),AUC.rf(j,:,trial));
        if length(BestIdx)>1
            BestIdx = BestIdx(end);
        end
        
        Scores = rerf_classprob(rf{BestIdx},Xtest{j},'last');
        Predictions = predict_class(Scores,rf{BestIdx}.classname);
        TestError.rf(trial,j) = misclassification_rate(Predictions,...
            Ytest{j},false);
        
        clear rf
        
        fprintf('RerF\n')
        for i = 1:length(mtrys)
            
            mtry = mtrys(i);
            dprime = ceil(mtry^(1/interp1(ps,slope,d)));
            
            fprintf('mtry = %d\n',mtry)

            tic;
            rerf{i} = rpclassificationforest(ntrees,Xtrain{j}(:,:,trial),...
                Ytrain{j}(:,trial),'sparsemethod','sparse-adjusted',...
                'nvartosample',mtry,'NWorkers',NWorkers,...
                'Stratified',true,'dprime',dprime);
            trainTime.rerf(j,i,trial) = toc;
            Scores = rerf_oob_classprob(rerf{i},Xtrain{j}(:,:,trial),'last');
            Predictions = predict_class(Scores,rerf{i}.classname);
            BaggedError.rerf(j,i,trial) = misclassification_rate(Predictions,...
                Ytrain{j}(:,trial),false);
            [~,~,~,AUC.rerf(j,i,trial)] = ...
                perfcurve(Ytrain{j}(:,trial),Scores(:,2),'1');
        end
        
        BestIdx = hp_optimize(BaggedError.rerf(j,:,trial),AUC.rerf(j,:,trial));
        if length(BestIdx)>1
            BestIdx = BestIdx(end);
        end
        
        Scores = rerf_classprob(rerf{BestIdx},Xtest{j},'last');
        Predictions = predict_class(Scores,rerf{BestIdx}.classname);
        TestError.rerf(trial,j) = misclassification_rate(Predictions,...
            Ytest{j},false);
        
        clear rerf
        
        fprintf('RerF(r)\n')
        for i = 1:length(mtrys)
            
            mtry = mtrys(i);
            dprime = ceil(mtry^(1/interp1(ps,slope,d)));

            fprintf('mtry = %d\n',mtry)
            
            tic;
            rerfr{i} = rpclassificationforest(ntrees,Xtrain{j}(:,:,trial),...
                Ytrain{j}(:,trial),'sparsemethod','sparse-adjusted',...
                'Robust',true,'nvartosample',mtry,'NWorkers',NWorkers,...
                'Stratified',true,'dprime',dprime);
            trainTime.rerfr(j,i,trial) = toc;
            Scores = rerf_oob_classprob(rerfr{i},Xtrain{j}(:,:,trial),'last');
            Predictions = predict_class(Scores,rerfr{i}.classname);
            BaggedError.rerfr(j,i,trial) = misclassification_rate(Predictions,...
                Ytrain{j}(:,trial),false);
            [~,~,~,AUC.rerfr(j,i,trial)] = ...
                perfcurve(Ytrain{j}(:,trial),Scores(:,2),'1');
        end
        
        BestIdx = hp_optimize(BaggedError.rerfr(j,:,trial),AUC.rerfr(j,:,trial));
        if length(BestIdx)>1
            BestIdx = BestIdx(end);
        end
        
        Scores = rerf_classprob(rerfr{BestIdx},Xtest{j},'last',Xtrain{j}(:,:,trial));
        Predictions = predict_class(Scores,rerfr{BestIdx}.classname);
        TestError.rerfr(trial,j) = misclassification_rate(Predictions,...
            Ytest{j},false);
        
        clear rerfr
        
        fprintf('RerF(d)\n')
        for i = 1:length(mtrys)

            mtry = mtrys(i);
            dprime = ceil(mtry^(1/interp1(ps,slope,d)));

            fprintf('mtry = %d\n',mtry)

            tic;
            rerfdn{i} = rpclassificationforest(ntrees,Xtrain{j}(:,:,trial),...
                Ytrain{j}(:,trial),'sparsemethod','sparse-adjusted',...
                'mdiff','node','nvartosample',mtry,'NWorkers',NWorkers,...
                'Stratified',true,'dprime',dprime);
            trainTime.rerfdn(j,i,trial) = toc;
            Scores = rerf_oob_classprob(rerfdn{i},Xtrain{j}(:,:,trial),'last');
            Predictions = predict_class(Scores,rerfdn{i}.classname);
            BaggedError.rerfdn(j,i,trial) = misclassification_rate(Predictions,...
                Ytrain{j}(:,trial),false);
            [~,~,~,AUC.rerfdn(j,i,trial)] = ...
                perfcurve(Ytrain{j}(:,trial),Scores(:,2),'1');
        end
        
        BestIdx = hp_optimize(BaggedError.rerfdn(j,:,trial),AUC.rf(j,:,trial));
        if length(BestIdx)>1
            BestIdx = BestIdx(end);
        end
        
        Scores = rerf_classprob(rerfdn{BestIdx},Xtest{j},'last');
        Predictions = predict_class(Scores,rerfdn{BestIdx}.classname);
        TestError.rerfdn(trial,j) = misclassification_rate(Predictions,...
            Ytest{j},false);
        
        clear rerfdn
        
        fprintf('RR-RF\n')
        for i = 1:length(mtrys)

            mtry = mtrys(i);

            fprintf('mtry = %d\n',mtry)

            if mtry <= d && d <= 500
                tic;
                rf_rot{i} = rpclassificationforest(ntrees,Xtrain{j}(:,:,trial),...
                    Ytrain{j}(:,trial),'RandomForest',true,'rotate',true,...
                    'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
                trainTime.rf_rot(j,i,trial) = toc;
                Scores = rerf_oob_classprob(rf_rot{i},Xtrain{j}(:,:,trial),'last');
                Predictions = predict_class(Scores,rf_rot{i}.classname);
                BaggedError.rf_rot(j,i,trial) = misclassification_rate(Predictions,...
                    Ytrain{j}(:,trial),false);
                [~,~,~,AUC.rf_rot(j,i,trial)] = ...
                    perfcurve(Ytrain{j}(:,trial),Scores(:,2),'1');
            end
        end
        
        BestIdx = hp_optimize(BaggedError.rf_rot(j,:,trial),AUC.rf_rot(j,:,trial));
        if length(BestIdx)>1
            BestIdx = BestIdx(end);
        end
        
        Scores = rerf_classprob(rf_rot{BestIdx},Xtest{j},'last');
        Predictions = predict_class(Scores,rf_rot{BestIdx}.classname);
        TestError.rf_rot(trial,j) = misclassification_rate(Predictions,...
            Ytest{j},false);
        
        clear rf_rot
        
        fprintf('FR-C\n')
        for i = 1:length(mtrys)

            mtry = mtrys(i);

            fprintf('mtry = %d\n',mtry)

            tic;
            frc{i} = rpclassificationforest(ntrees,Xtrain{j}(:,:,trial),...
                Ytrain{j}(:,trial),'sparsemethod','frc','nmix',nmix,...
                'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            trainTime.frc(j,i,trial) = toc;
            Scores = rerf_oob_classprob(frc{i},Xtrain{j}(:,:,trial),'last');
            Predictions = predict_class(Scores,frc{i}.classname);
            BaggedError.frc(j,i,trial) = misclassification_rate(Predictions,...
                Ytrain{j}(:,trial),false);
            [~,~,~,AUC.frc(j,i,trial)] = ...
                perfcurve(Ytrain{j}(:,trial),Scores(:,2),'1');
        end
        
        BestIdx = hp_optimize(BaggedError.frc(j,:,trial),AUC.frc(j,:,trial));
        if length(BestIdx)>1
            BestIdx = BestIdx(end);
        end
        
        Scores = rerf_classprob(frc{BestIdx},Xtest{j},'last');
        Predictions = predict_class(Scores,frc{BestIdx}.classname);
        TestError.frc(trial,j) = misclassification_rate(Predictions,...
            Ytest{j},false);
        
        clear frc
    end
    save([rerfPath 'RandomerForest/Results/Trunk.mat'],'dims',...
        'BaggedError','AUC','TestError','trainTime')
end