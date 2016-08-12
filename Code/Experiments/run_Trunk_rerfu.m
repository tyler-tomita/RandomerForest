close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

rng(1);

load Trunk_data
load Random_matrix_adjustment_factor

for i = length(dims):length(dims)
    p = dims(i);
    fprintf('p = %d\n',p)
      
    if p <= 5
        mtrys = [1:p ceil(p.^[1.5 2])];
    elseif p > 5 && p <= 100
        mtrys = ceil(p.^[0 1/4 1/2 3/4 1 1.5 2]);
    else
        mtrys = [ceil(p.^[0 1/4 1/2 3/4 1]) 5*p 10*p];
    end

    Params{i}.rerfu.nTrees = 1000;
    Params{i}.rerfu.Stratified = true;
    Params{i}.rerfu.NWorkers = 16;
    Params{i}.rerfu.d = mtrys;
    Params{i}.rerfu.Method = 'uniform-nnzs';
    if p <5
        Params{i}.rerfu.nmix = 1:p;
    else
        Params{i}.rerfu.nmix = 1:5;
    end

    OOBError{i}.rerfu = NaN(ntrials,length(Params{i}.rerfu.d));
    OOBAUC{i}.rerfu = NaN(ntrials,length(Params{i}.rerfu.d));
    TrainTime{i}.rerfu = NaN(ntrials,length(Params{i}.rerfu.d));

    for trial = 1:ntrials
        
        fprintf('trial %d\n',trial)

        % train classifier
        poolobj = gcp('nocreate');
        if isempty(poolobj)
            parpool('local',Params{i}.rerfu.NWorkers,...
                'IdleTimeout',360);
        end

        [Forest,~,TrainTime{i}.rerfu(trial,:)] = ...
            RerF_train(Xtrain{i}(:,:,trial),Ytrain{i}(:,trial),...
            Params{i}.rerfu);

        % select best hyperparameter

        for j = 1:length(Params{i}.rerfu.d)
            Scores = rerf_oob_classprob(Forest{j},Xtrain{i}(:,:,trial),'last');
            Predictions = predict_class(Scores,Forest{j}.classname);
            OOBError{i}.rerfu(trial,j) = ...
                misclassification_rate(Predictions,Ytrain{i}(:,trial),...
                false);
            if size(Scores,2) > 2
                Yb = binarize_labels(Ytrain{i}(:,trial),Forest{j}.classname);
                [~,~,~,OOBAUC{i}.rerfu(trial,j)] = ...
                    perfcurve(Yb(:),Scores(:),'1');
            else
                [~,~,~,OOBAUC{i}.rerfu(trial,j)] = ...
                    perfcurve(Ytrain{i}(:,trial),Scores(:,2),'1');
            end
        end
        BestIdx = hp_optimize(OOBError{i}.rerfu(trial,:),...
            OOBAUC{i}.rerfu(trial,:));
        if length(BestIdx)>1
            BestIdx = BestIdx(end);
        end

        if ~Forest{BestIdx}.Robust
            Scores = rerf_classprob(Forest{BestIdx},Xtest{i},'last');
        else
            Scores = rerf_classprob(Forest{BestIdx},Xtest{i},...
                'last',Xtrain{i}(:,:,trial));
        end
        Predictions = predict_class(Scores,Forest{BestIdx}.classname);
        TestError{i}.rerfu(trial) = misclassification_rate(Predictions,...
            Ytest{i},false);

        clear Forest

        save([rerfPath 'RandomerForest/Results/Trunk_rerfu_p1000.mat'],'dims',...
            'OOBError','OOBAUC','TestError','TrainTime')
    end
end