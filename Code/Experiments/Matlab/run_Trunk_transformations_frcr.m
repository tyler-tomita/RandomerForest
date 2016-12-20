close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

rng(1);

load Trunk_transformations_data
ntrees = 1000;
NWorkers = 24;
nmix = 2;

for i = 1:length(dims)
    
    d = dims(i);
    fprintf('dimension %d\n',d)
    if d <= 5
        mtrys = 1:d;
    else
        mtrys = ceil(d.^[0 1/4 1/2 3/4 1]);
    end

    err_frcr.Untransformed = zeros(ntrees,length(mtrys),ntrials);
    err_frcr.Rotated = zeros(ntrees,length(mtrys),ntrials);
    err_frcr.Scaled = zeros(ntrees,length(mtrys),ntrials);
    err_frcr.Affine = zeros(ntrees,length(mtrys),ntrials);
    err_frcr.Outlier = zeros(ntrees,length(mtrys),ntrials);
    
    Class = [0;1];

    for trial = 1:ntrials

        fprintf('trial %d\n',trial)

        for j = 1:length(mtrys)
            mtry = mtrys(j);

            fprintf('mtry = %d\n',mtry)
            
            %Untransformed
            frcr = rpclassificationforest(ntrees,X{i}(:,:,trial),Y{i}(:,trial),'sparsemethod','frc','nmix',nmix,'Robust',true,'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            Predictions = oobpredict(frcr,X{i}(:,:,trial),Y{i}(:,trial));
            err_frcr.Untransformed(:,j,trial) = oob_error(Predictions,Y{i}(:,trial),'every');

            %Rotated
            frcr = rpclassificationforest(ntrees,X_rot{i}(:,:,trial),Y{i}(:,trial),'sparsemethod','frc','nmix',nmix,'Robust',true,'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            Predictions = oobpredict(frcr,X_rot{i}(:,:,trial),Y{i}(:,trial));
            err_frcr.Rotated(:,j,trial) = oob_error(Predictions,Y{i}(:,trial),'every');

            %Scaled
            frcr = rpclassificationforest(ntrees,X_scale{i}(:,:,trial),Y{i}(:,trial),'sparsemethod','frc','nmix',nmix,'Robust',true,'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            Predictions = oobpredict(frcr,X_scale{i}(:,:,trial),Y{i}(:,trial));
            err_frcr.Scaled(:,j,trial) = oob_error(Predictions,Y{i}(:,trial),'every');

            %Affine
            frcr = rpclassificationforest(ntrees,X_affine{i}(:,:,trial),Y{i}(:,trial),'sparsemethod','frc','nmix',nmix,'Robust',true,'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            Predictions = oobpredict(frcr,X_affine{i}(:,:,trial),Y{i}(:,trial));
            err_frcr.Affine(:,j,trial) = oob_error(Predictions,Y{i}(:,trial),'every');

            %Outlier
            frcr = rpclassificationforest(ntrees,X_out{i}(:,:,trial),Y_out{i}(:,trial),'sparsemethod','frc','nmix',nmix,'Robust',true,'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            Predictions = oobpredict(frcr,X_out{i}(:,:,trial),Y_out{i}(:,trial));
            err_frcr.Outlier(:,j,trial) = oob_error(Predictions,Y_out{i}(:,trial),'every');
        end
    end

    if length(mtrys) < 5
        emptyCol = 5 - length(mtrys);
        
        %Untransformed
        sem_frcr.Untransformed(:,1:length(mtrys),i) = std(err_frcr.Untransformed,[],3)/sqrt(ntrials);
        sem_frcr.Untransformed(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);

        var_frcr.Untransformed(:,1:length(mtrys),i) = var(err_frcr.Untransformed,0,3);
        var_frcr.Untransformed(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);

        mean_err_frcr.Untransformed(:,1:length(mtrys),i) = mean(err_frcr.Untransformed,3);
        mean_err_frcr.Untransformed(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        
        %Rotated
        sem_frcr.Rotated(:,1:length(mtrys),i) = std(err_frcr.Rotated,[],3)/sqrt(ntrials);
        sem_frcr.Rotated(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);

        var_frcr.Rotated(:,1:length(mtrys),i) = var(err_frcr.Rotated,0,3);
        var_frcr.Rotated(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);

        mean_err_frcr.Rotated(:,1:length(mtrys),i) = mean(err_frcr.Rotated,3);
        mean_err_frcr.Rotated(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        
        %Scaled
        sem_frcr.Scaled(:,1:length(mtrys),i) = std(err_frcr.Scaled,[],3)/sqrt(ntrials);
        sem_frcr.Scaled(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);

        var_frcr.Scaled(:,1:length(mtrys),i) = var(err_frcr.Scaled,0,3);
        var_frcr.Scaled(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);

        mean_err_frcr.Scaled(:,1:length(mtrys),i) = mean(err_frcr.Scaled,3);
        mean_err_frcr.Scaled(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        
        %Affine
        sem_frcr.Affine(:,1:length(mtrys),i) = std(err_frcr.Affine,[],3)/sqrt(ntrials);
        sem_frcr.Affine(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);

        var_frcr.Affine(:,1:length(mtrys),i) = var(err_frcr.Affine,0,3);
        var_frcr.Affine(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);

        mean_err_frcr.Affine(:,1:length(mtrys),i) = mean(err_frcr.Affine,3);
        mean_err_frcr.Affine(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        
        %Outlier
        sem_frcr.Outlier(:,1:length(mtrys),i) = std(err_frcr.Outlier,[],3)/sqrt(ntrials);
        sem_frcr.Outlier(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);

        var_frcr.Outlier(:,1:length(mtrys),i) = var(err_frcr.Outlier,0,3);
        var_frcr.Outlier(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);

        mean_err_frcr.Outlier(:,1:length(mtrys),i) = mean(err_frcr.Outlier,3);
        mean_err_frcr.Outlier(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
    else
        %Untransformed
        sem_frcr.Untransformed(:,:,i) = std(err_frcr.Untransformed,[],3)/sqrt(ntrials);

        var_frcr.Untransformed(:,:,i) = var(err_frcr.Untransformed,0,3);

        mean_err_frcr.Untransformed(:,:,i) = mean(err_frcr.Untransformed,3);
        
        %Rotated
        sem_frcr.Rotated(:,:,i) = std(err_frcr.Rotated,[],3)/sqrt(ntrials);

        var_frcr.Rotated(:,:,i) = var(err_frcr.Rotated,0,3);

        mean_err_frcr.Rotated(:,:,i) = mean(err_frcr.Rotated,3);
        
        %Scaled
        sem_frcr.Scaled(:,:,i) = std(err_frcr.Scaled,[],3)/sqrt(ntrials);

        var_frcr.Scaled(:,:,i) = var(err_frcr.Scaled,0,3);

        mean_err_frcr.Scaled(:,:,i) = mean(err_frcr.Scaled,3);
        
        %Affine
        sem_frcr.Affine(:,:,i) = std(err_frcr.Affine,[],3)/sqrt(ntrials);

        var_frcr.Affine(:,:,i) = var(err_frcr.Affine,0,3);

        mean_err_frcr.Affine(:,:,i) = mean(err_frcr.Affine,3);
        
        %Outlier
        sem_frcr.Outlier(:,:,i) = std(err_frcr.Outlier,[],3)/sqrt(ntrials);

        var_frcr.Outlier(:,:,i) = var(err_frcr.Outlier,0,3);

        mean_err_frcr.Outlier(:,:,i) = mean(err_frcr.Outlier,3);
    end
    
    save([rerfPath 'RandomerForest/Results/Trunk_transformations_frcr.mat'],...
        'dims','sem_frcr','var_frcr','mean_err_frcr')
end
