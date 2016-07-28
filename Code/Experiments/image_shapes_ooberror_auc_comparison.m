close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

rng(1);

load image_shapes_data

[ih,iw,n] = size(X_image);
p = ih*iw;
X = reshape(X_image,p,n)';
Ystr = cellstr(num2str(Y));

Xtrain = X(TrainIdx{1}(1,:),:);
Ytrain = Y(TrainIdx{1}(1,:));
Ytrain_str = Ystr(TrainIdx{1}(1,:));

clear X_image

ntrees = 2000;
Stratified = true;
NWorkers = 24;
ntrials = 3;

nTrain = ns(1);

ds = ceil(p.^[0 1/4 1/2 3/4 1]);

BaggedError = NaN(ntrials,length(ds));
AUC = NaN(ntrials,length(ds));
for trial = 1:ntrials
    fprintf('trial %d\n',trial)
    for i = 1:length(ds)
        d = ds(i);

        fprintf('d = %d\n',d)

        rf = rpclassificationforest(ntrees,Xtrain,Ytrain_str,...
            'RandomForest',true,'nvartosample',d,'NWorkers',NWorkers,...
            'Stratified',Stratified);
        Scores = rerf_oob_classprob(rf,Xtrain,'last');
        Predictions = predict_class(Scores,rf.classname);
        BaggedError(trial,i) = misclassification_rate(Predictions,Ytrain_str);
        [~,~,~,AUC(trial,i)] = perfcurve(Ytrain_str,Scores(:,2),'1');
    end

    save([rerfPath 'RandomerForest/Results/image_shapes_ooberror_auc_comparison.mat'],...
        'BaggedError','AUC')
end