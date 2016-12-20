close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

rng(1);

load Trunk_transformations_data
ntrees = 1000;
NWorkers = 16;


for i = 1:length(dims)
    
    d = dims(i);
    fprintf('dimension %d\n',d)
    if d <= 5
        mtrys = 1:d;
    else
        mtrys = ceil(d.^[0 1/4 1/2 3/4 1]);
    end
    
    if d >= 6
        nmixs = 2:6;
    else
        nmixs = 2:d;
    end 

    err_rfr.Untransformed = zeros(ntrees,length(mtrys),ntrials);
    err_rf_rotr.Untransformed = zeros(ntrees,length(mtrys),ntrials);
    
    err_rfr.Rotated = zeros(ntrees,length(mtrys),ntrials);
    err_rf_rotr.Rotated = zeros(ntrees,length(mtrys),ntrials);
    
    err_rfr.Scaled = zeros(ntrees,length(mtrys),ntrials);
    err_rf_rotr.Scaled = zeros(ntrees,length(mtrys),ntrials);
    
    err_rfr.Affine = zeros(ntrees,length(mtrys),ntrials);
    err_rf_rotr.Affine = zeros(ntrees,length(mtrys),ntrials);
    
    err_rfr.Outlier = zeros(ntrees,length(mtrys),ntrials);
    err_rf_rotr.Outlier = zeros(ntrees,length(mtrys),ntrials);
    
    Class = [0;1];

    for trial = 1:ntrials

        fprintf('trial %d\n',trial)

        for j = 1:length(mtrys)
            mtry = mtrys(j);

            fprintf('mtry = %d\n',mtry)
            
            %Untransformed
            rfr = rpclassificationforest(ntrees,X{i}(:,:,trial),Y{i}(:,trial),'RandomForest',true,'Robust',true,'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            Predictions = oobpredict(rfr,X{i}(:,:,trial),Y{i}(:,trial));
            err_rfr.Untransformed(:,j,trial) = oob_error(Predictions,Y{i}(:,trial),'every');

            rf_rotr = rpclassificationforest(ntrees,X{i}(:,:,trial),Y{i}(:,trial),'RandomForest',true,'Robust',true,'rotate',true,'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            Predictions = oobpredict(rf_rotr,X{i}(:,:,trial),Y{i}(:,trial));
            err_rf_rotr.Untransformed(:,j,trial) = oob_error(Predictions,Y{i}(:,trial),'every');

            %Rotated
            rfr = rpclassificationforest(ntrees,X_rot{i}(:,:,trial),Y{i}(:,trial),'RandomForest',true,'Robust',true,'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            Predictions = oobpredict(rfr,X_rot{i}(:,:,trial),Y{i}(:,trial));
            err_rfr.Rotated(:,j,trial) = oob_error(Predictions,Y{i}(:,trial),'every');

            rf_rotr = rpclassificationforest(ntrees,X_rot{i}(:,:,trial),Y{i}(:,trial),'RandomForest',true,'Robust',true,'rotate',true,'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            Predictions = oobpredict(rf_rotr,X_rot{i}(:,:,trial),Y{i}(:,trial));
            err_rf_rotr.Rotated(:,j,trial) = oob_error(Predictions,Y{i}(:,trial),'every');

            %Scaled
            rfr = rpclassificationforest(ntrees,X_scale{i}(:,:,trial),Y{i}(:,trial),'RandomForest',true,'Robust',true,'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            Predictions = oobpredict(rfr,X_scale{i}(:,:,trial),Y{i}(:,trial));
            err_rfr.Scaled(:,j,trial) = oob_error(Predictions,Y{i}(:,trial),'every');

            rf_rotr = rpclassificationforest(ntrees,X_scale{i}(:,:,trial),Y{i}(:,trial),'RandomForest',true,'Robust',true,'rotate',true,'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            Predictions = oobpredict(rf_rotr,X_scale{i}(:,:,trial),Y{i}(:,trial));
            err_rf_rotr.Scaled(:,j,trial) = oob_error(Predictions,Y{i}(:,trial),'every');

            %Affine
            rfr = rpclassificationforest(ntrees,X_affine{i}(:,:,trial),Y{i}(:,trial),'RandomForest',true,'Robust',true,'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            Predictions = oobpredict(rfr,X_affine{i}(:,:,trial),Y{i}(:,trial));
            err_rfr.Affine(:,j,trial) = oob_error(Predictions,Y{i}(:,trial),'every');

            rf_rotr = rpclassificationforest(ntrees,X_affine{i}(:,:,trial),Y{i}(:,trial),'RandomForest',true,'Robust',true,'rotate',true,'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            Predictions = oobpredict(rf_rotr,X_affine{i}(:,:,trial),Y{i}(:,trial));
            err_rf_rotr.Affine(:,j,trial) = oob_error(Predictions,Y{i}(:,trial),'every');

            %Outlier
            rfr = rpclassificationforest(ntrees,X_out{i}(:,:,trial),Y_out{i}(:,trial),'RandomForest',true,'Robust',true,'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            Predictions = oobpredict(rfr,X_out{i}(:,:,trial),Y_out{i}(:,trial));
            err_rfr.Outlier(:,j,trial) = oob_error(Predictions,Y_out{i}(:,trial),'every');

            rf_rotr = rpclassificationforest(ntrees,X_out{i}(:,:,trial),Y_out{i}(:,trial),'RandomForest',true,'Robust',true,'rotate',true,'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            Predictions = oobpredict(rf_rotr,X_out{i}(:,:,trial),Y_out{i}(:,trial));
            err_rf_rotr.Outlier(:,j,trial) = oob_error(Predictions,Y_out{i}(:,trial),'every');
        end
    end

    if length(mtrys) < 5
        emptyCol = 5 - length(mtrys);
        
        %Untransformed
        sem_rfr.Untransformed(:,1:length(mtrys),i) = std(err_rfr.Untransformed,[],3)/sqrt(ntrials);
        sem_rf_rotr.Untransformed(:,1:length(mtrys),i) = std(err_rf_rotr.Untransformed,[],3)/sqrt(ntrials);
        sem_rfr.Untransformed(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        sem_rf_rotr.Untransformed(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);

        var_rfr.Untransformed(:,1:length(mtrys),i) = var(err_rfr.Untransformed,0,3);
        var_rf_rotr.Untransformed(:,1:length(mtrys),i) = var(err_rf_rotr.Untransformed,0,3);
        var_rfr.Untransformed(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        var_rf_rotr.Untransformed(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);

        mean_err_rfr.Untransformed(:,1:length(mtrys),i) = mean(err_rfr.Untransformed,3);
        mean_err_rf_rotr.Untransformed(:,1:length(mtrys),i) = mean(err_rf_rotr.Untransformed,3);
        mean_err_rfr.Untransformed(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        mean_err_rf_rotr.Untransformed(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        
        %Rotated
        sem_rfr.Rotated(:,1:length(mtrys),i) = std(err_rfr.Rotated,[],3)/sqrt(ntrials);
        sem_rf_rotr.Rotated(:,1:length(mtrys),i) = std(err_rf_rotr.Rotated,[],3)/sqrt(ntrials);
        sem_rfr.Rotated(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        sem_rf_rotr.Rotated(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);

        var_rfr.Rotated(:,1:length(mtrys),i) = var(err_rfr.Rotated,0,3);
        var_rf_rotr.Rotated(:,1:length(mtrys),i) = var(err_rf_rotr.Rotated,0,3);
        var_rfr.Rotated(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        var_rf_rotr.Rotated(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);

        mean_err_rfr.Rotated(:,1:length(mtrys),i) = mean(err_rfr.Rotated,3);
        mean_err_rf_rotr.Rotated(:,1:length(mtrys),i) = mean(err_rf_rotr.Rotated,3);
        mean_err_rfr.Rotated(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        mean_err_rf_rotr.Rotated(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        
        %Scaled
        sem_rfr.Scaled(:,1:length(mtrys),i) = std(err_rfr.Scaled,[],3)/sqrt(ntrials);
        sem_rf_rotr.Scaled(:,1:length(mtrys),i) = std(err_rf_rotr.Scaled,[],3)/sqrt(ntrials);
        sem_rfr.Scaled(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        sem_rf_rotr.Scaled(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);

        var_rfr.Scaled(:,1:length(mtrys),i) = var(err_rfr.Scaled,0,3);
        var_rf_rotr.Scaled(:,1:length(mtrys),i) = var(err_rf_rotr.Scaled,0,3);
        var_rfr.Scaled(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        var_rf_rotr.Scaled(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);

        mean_err_rfr.Scaled(:,1:length(mtrys),i) = mean(err_rfr.Scaled,3);
        mean_err_rf_rotr.Scaled(:,1:length(mtrys),i) = mean(err_rf_rotr.Scaled,3);
        mean_err_rfr.Scaled(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        mean_err_rf_rotr.Scaled(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        
        %Affine
        sem_rfr.Affine(:,1:length(mtrys),i) = std(err_rfr.Affine,[],3)/sqrt(ntrials);
        sem_rf_rotr.Affine(:,1:length(mtrys),i) = std(err_rf_rotr.Affine,[],3)/sqrt(ntrials);
        sem_rfr.Affine(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        sem_rf_rotr.Affine(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);

        var_rfr.Affine(:,1:length(mtrys),i) = var(err_rfr.Affine,0,3);
        var_rf_rotr.Affine(:,1:length(mtrys),i) = var(err_rf_rotr.Affine,0,3);
        var_rfr.Affine(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        var_rf_rotr.Affine(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);

        mean_err_rfr.Affine(:,1:length(mtrys),i) = mean(err_rfr.Affine,3);
        mean_err_rf_rotr.Affine(:,1:length(mtrys),i) = mean(err_rf_rotr.Affine,3);
        mean_err_rfr.Affine(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        mean_err_rf_rotr.Affine(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        
        %Outlier
        sem_rfr.Outlier(:,1:length(mtrys),i) = std(err_rfr.Outlier,[],3)/sqrt(ntrials);
        sem_rf_rotr.Outlier(:,1:length(mtrys),i) = std(err_rf_rotr.Outlier,[],3)/sqrt(ntrials);
        sem_rfr.Outlier(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        sem_rf_rotr.Outlier(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);

        var_rfr.Outlier(:,1:length(mtrys),i) = var(err_rfr.Outlier,0,3);
        var_rf_rotr.Outlier(:,1:length(mtrys),i) = var(err_rf_rotr.Outlier,0,3);
        var_rfr.Outlier(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        var_rf_rotr.Outlier(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);

        mean_err_rfr.Outlier(:,1:length(mtrys),i) = mean(err_rfr.Outlier,3);
        mean_err_rf_rotr.Outlier(:,1:length(mtrys),i) = mean(err_rf_rotr.Outlier,3);
        mean_err_rfr.Outlier(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        mean_err_rf_rotr.Outlier(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
    else
        %Untransformed
        sem_rfr.Untransformed(:,:,i) = std(err_rfr.Untransformed,[],3)/sqrt(ntrials);
        sem_rf_rotr.Untransformed(:,:,i) = std(err_rf_rotr.Untransformed,[],3)/sqrt(ntrials);

        var_rfr.Untransformed(:,:,i) = var(err_rfr.Untransformed,0,3);
        var_rf_rotr.Untransformed(:,:,i) = var(err_rf_rotr.Untransformed,0,3);

        mean_err_rfr.Untransformed(:,:,i) = mean(err_rfr.Untransformed,3);
        mean_err_rf_rotr.Untransformed(:,:,i) = mean(err_rf_rotr.Untransformed,3);
        
        %Rotated
        sem_rfr.Rotated(:,:,i) = std(err_rfr.Rotated,[],3)/sqrt(ntrials);
        sem_rf_rotr.Rotated(:,:,i) = std(err_rf_rotr.Rotated,[],3)/sqrt(ntrials);

        var_rfr.Rotated(:,:,i) = var(err_rfr.Rotated,0,3);
        var_rf_rotr.Rotated(:,:,i) = var(err_rf_rotr.Rotated,0,3);

        mean_err_rfr.Rotated(:,:,i) = mean(err_rfr.Rotated,3);
        mean_err_rf_rotr.Rotated(:,:,i) = mean(err_rf_rotr.Rotated,3);
        
        %Scaled
        sem_rfr.Scaled(:,:,i) = std(err_rfr.Scaled,[],3)/sqrt(ntrials);
        sem_rf_rotr.Scaled(:,:,i) = std(err_rf_rotr.Scaled,[],3)/sqrt(ntrials);

        var_rfr.Scaled(:,:,i) = var(err_rfr.Scaled,0,3);
        var_rf_rotr.Scaled(:,:,i) = var(err_rf_rotr.Scaled,0,3);

        mean_err_rfr.Scaled(:,:,i) = mean(err_rfr.Scaled,3);
        mean_err_rf_rotr.Scaled(:,:,i) = mean(err_rf_rotr.Scaled,3);
        
        %Affine
        sem_rfr.Affine(:,:,i) = std(err_rfr.Affine,[],3)/sqrt(ntrials);
        sem_rf_rotr.Affine(:,:,i) = std(err_rf_rotr.Affine,[],3)/sqrt(ntrials);

        var_rfr.Affine(:,:,i) = var(err_rfr.Affine,0,3);
        var_rf_rotr.Affine(:,:,i) = var(err_rf_rotr.Affine,0,3);

        mean_err_rfr.Affine(:,:,i) = mean(err_rfr.Affine,3);
        mean_err_rf_rotr.Affine(:,:,i) = mean(err_rf_rotr.Affine,3);
        
        %Outlier
        sem_rfr.Outlier(:,:,i) = std(err_rfr.Outlier,[],3)/sqrt(ntrials);
        sem_rf_rotr.Outlier(:,:,i) = std(err_rf_rotr.Outlier,[],3)/sqrt(ntrials);

        var_rfr.Outlier(:,:,i) = var(err_rfr.Outlier,0,3);
        var_rf_rotr.Outlier(:,:,i) = var(err_rf_rotr.Outlier,0,3);

        mean_err_rfr.Outlier(:,:,i) = mean(err_rfr.Outlier,3);
        mean_err_rf_rotr.Outlier(:,:,i) = mean(err_rf_rotr.Outlier,3);
    end
    
    save([rerfPath 'RandomerForest/Results/Trunk_transformations_rf_rotr.mat'],...
        'dims','sem_rfr','var_rfr','mean_err_rfr','sem_rf_rotr','var_rf_rotr','mean_err_rf_rotr')
end
