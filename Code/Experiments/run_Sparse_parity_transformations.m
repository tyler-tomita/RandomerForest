close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

rng(1);

load Sparse_parity_transformations_data
ntrees = 500;
NWorkers = 24;

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

    err_rf.Untransformed = zeros(ntrees,length(mtrys),ntrials);
    err_rerf.Untransformed = zeros(ntrees,length(mtrys),ntrials);
    err_rerfdn.Untransformed = zeros(ntrees,length(mtrys),ntrials);
    err_rerfr.Untransformed = zeros(ntrees,length(mtrys),ntrials);
    err_rf_rot.Untransformed = zeros(ntrees,length(mtrys),ntrials);
    
    err_rf.Rotated = zeros(ntrees,length(mtrys),ntrials);
    err_rerf.Rotated = zeros(ntrees,length(mtrys),ntrials);
    err_rerfdn.Rotated = zeros(ntrees,length(mtrys),ntrials);
    err_rerfr.Rotated = zeros(ntrees,length(mtrys),ntrials);
    err_rf_rot.Rotated = zeros(ntrees,length(mtrys),ntrials);
    
    err_rf.Scaled = zeros(ntrees,length(mtrys),ntrials);
    err_rerf.Scaled = zeros(ntrees,length(mtrys),ntrials);
    err_rerfdn.Scaled = zeros(ntrees,length(mtrys),ntrials);
    err_rerfr.Scaled = zeros(ntrees,length(mtrys),ntrials);
    err_rf_rot.Scaled = zeros(ntrees,length(mtrys),ntrials);
    
    err_rf.Affine = zeros(ntrees,length(mtrys),ntrials);
    err_rerf.Affine = zeros(ntrees,length(mtrys),ntrials);
    err_rerfdn.Affine = zeros(ntrees,length(mtrys),ntrials);
    err_rerfr.Affine = zeros(ntrees,length(mtrys),ntrials);
    err_rf_rot.Affine = zeros(ntrees,length(mtrys),ntrials);
    
    err_rf.Outlier = zeros(ntrees,length(mtrys),ntrials);
    err_rerf.Outlier = zeros(ntrees,length(mtrys),ntrials);
    err_rerfdn.Outlier = zeros(ntrees,length(mtrys),ntrials);
    err_rerfr.Outlier = zeros(ntrees,length(mtrys),ntrials);
    err_rf_rot.Outlier = zeros(ntrees,length(mtrys),ntrials);
    
    err_frc.Untransformed = zeros(ntrees,length(mtrys)*length(nmixs),ntrials);
    
    err_frc.Rotated = zeros(ntrees,length(mtrys)*length(nmixs),ntrials);
    
    err_frc.Scaled = zeros(ntrees,length(mtrys)*length(nmixs),ntrials);
    
    err_frc.Affine = zeros(ntrees,length(mtrys)*length(nmixs),ntrials);
    
    err_frc.Outlier = zeros(ntrees,length(mtrys)*length(nmixs),ntrials);
    
    Class = [0;1];

    for trial = 1:ntrials

        fprintf('trial %d\n',trial)

        for j = 1:length(mtrys)
            mtry = mtrys(j);

            fprintf('mtry = %d\n',mtry)
            
            %Untransformed
            rf = rpclassificationforest(ntrees,X{i}(:,:,trial),Y{i}(:,trial),'RandomForest',true,'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            err_rf.Untransformed(:,j,trial) = oobpredict(rf,X{i}(:,:,trial),Y{i}(:,trial),'every');

            rerf = rpclassificationforest(ntrees,X{i}(:,:,trial),Y{i}(:,trial),'sparsemethod','sparse','nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            err_rerf.Untransformed(:,j,trial) = oobpredict(rerf,X{i}(:,:,trial),Y{i}(:,trial),'every');

            rerfdn = rpclassificationforest(ntrees,X{i}(:,:,trial),Y{i}(:,trial),'sparsemethod','sparse','mdiff','node','nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            err_rerfdn.Untransformed(:,j,trial) = oobpredict(rerfdn,X{i}(:,:,trial),Y{i}(:,trial),'every');

            rerfr = rpclassificationforest(ntrees,X{i}(:,:,trial),Y{i}(:,trial),'sparsemethod','sparse','Robust',true,'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            err_rerfr.Untransformed(:,j,trial) = oobpredict(rerfr,X{i}(:,:,trial),Y{i}(:,trial),'every');

            rf_rot = rpclassificationforest(ntrees,X{i}(:,:,trial),Y{i}(:,trial),'RandomForest',true,'rotate',true,'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            err_rf_rot.Untransformed(:,j,trial) = oobpredict(rf_rot,X{i}(:,:,trial),Y{i}(:,trial),'every');

            %Rotated
            rf = rpclassificationforest(ntrees,X_rot{i}(:,:,trial),Y{i}(:,trial),'RandomForest',true,'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            err_rf.Rotated(:,j,trial) = oobpredict(rf,X_rot{i}(:,:,trial),Y{i}(:,trial),'every');

            rerf = rpclassificationforest(ntrees,X_rot{i}(:,:,trial),Y{i}(:,trial),'sparsemethod','sparse','nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            err_rerf.Rotated(:,j,trial) = oobpredict(rerf,X_rot{i}(:,:,trial),Y{i}(:,trial),'every');

            rerfdn = rpclassificationforest(ntrees,X_rot{i}(:,:,trial),Y{i}(:,trial),'sparsemethod','sparse','mdiff','node','nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            err_rerfdn.Rotated(:,j,trial) = oobpredict(rerfdn,X_rot{i}(:,:,trial),Y{i}(:,trial),'every');

            rerfr = rpclassificationforest(ntrees,X_rot{i}(:,:,trial),Y{i}(:,trial),'sparsemethod','sparse','Robust',true,'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            err_rerfr.Rotated(:,j,trial) = oobpredict(rerfr,X_rot{i}(:,:,trial),Y{i}(:,trial),'every');

            rf_rot = rpclassificationforest(ntrees,X_rot{i}(:,:,trial),Y{i}(:,trial),'RandomForest',true,'rotate',true,'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            err_rf_rot.Rotated(:,j,trial) = oobpredict(rf_rot,X_rot{i}(:,:,trial),Y{i}(:,trial),'every');

            %Scaled
            rf = rpclassificationforest(ntrees,X_scale{i}(:,:,trial),Y{i}(:,trial),'RandomForest',true,'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            err_rf.Scaled(:,j,trial) = oobpredict(rf,X_scale{i}(:,:,trial),Y{i}(:,trial),'every');

            rerf = rpclassificationforest(ntrees,X_scale{i}(:,:,trial),Y{i}(:,trial),'sparsemethod','sparse','nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            err_rerf.Scaled(:,j,trial) = oobpredict(rerf,X_scale{i}(:,:,trial),Y{i}(:,trial),'every');

            rerfdn = rpclassificationforest(ntrees,X_scale{i}(:,:,trial),Y{i}(:,trial),'sparsemethod','sparse','mdiff','node','nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            err_rerfdn.Scaled(:,j,trial) = oobpredict(rerfdn,X_scale{i}(:,:,trial),Y{i}(:,trial),'every');

            rerfr = rpclassificationforest(ntrees,X_scale{i}(:,:,trial),Y{i}(:,trial),'sparsemethod','sparse','Robust',true,'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            err_rerfr.Scaled(:,j,trial) = oobpredict(rerfr,X_scale{i}(:,:,trial),Y{i}(:,trial),'every');

            rf_rot = rpclassificationforest(ntrees,X_scale{i}(:,:,trial),Y{i}(:,trial),'RandomForest',true,'rotate',true,'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            err_rf_rot.Scaled(:,j,trial) = oobpredict(rf_rot,X_scale{i}(:,:,trial),Y{i}(:,trial),'every');

            %Affine
            rf = rpclassificationforest(ntrees,X_affine{i}(:,:,trial),Y{i}(:,trial),'RandomForest',true,'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            err_rf.Affine(:,j,trial) = oobpredict(rf,X_affine{i}(:,:,trial),Y{i}(:,trial),'every');

            rerf = rpclassificationforest(ntrees,X_affine{i}(:,:,trial),Y{i}(:,trial),'sparsemethod','sparse','nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            err_rerf.Affine(:,j,trial) = oobpredict(rerf,X_affine{i}(:,:,trial),Y{i}(:,trial),'every');

            rerfdn = rpclassificationforest(ntrees,X_affine{i}(:,:,trial),Y{i}(:,trial),'sparsemethod','sparse','mdiff','node','nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            err_rerfdn.Affine(:,j,trial) = oobpredict(rerfdn,X_affine{i}(:,:,trial),Y{i}(:,trial),'every');

            rerfr = rpclassificationforest(ntrees,X_affine{i}(:,:,trial),Y{i}(:,trial),'sparsemethod','sparse','Robust',true,'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            err_rerfr.Affine(:,j,trial) = oobpredict(rerfr,X_affine{i}(:,:,trial),Y{i}(:,trial),'every');

            rf_rot = rpclassificationforest(ntrees,X_affine{i}(:,:,trial),Y{i}(:,trial),'RandomForest',true,'rotate',true,'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            err_rf_rot.Affine(:,j,trial) = oobpredict(rf_rot,X_affine{i}(:,:,trial),Y{i}(:,trial),'every');

            %Outlier
            rf = rpclassificationforest(ntrees,X_out{i}(:,:,trial),Y_out{i}(:,trial),'RandomForest',true,'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            err_rf.Outlier(:,j,trial) = oobpredict(rf,X_out{i}(:,:,trial),Y_out{i}(:,trial),'every');

            rerf = rpclassificationforest(ntrees,X_out{i}(:,:,trial),Y_out{i}(:,trial),'sparsemethod','sparse','nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            err_rerf.Outlier(:,j,trial) = oobpredict(rerf,X_out{i}(:,:,trial),Y_out{i}(:,trial),'every');

            rerfdn = rpclassificationforest(ntrees,X_out{i}(:,:,trial),Y_out{i}(:,trial),'sparsemethod','sparse','mdiff','node','nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            err_rerfdn.Outlier(:,j,trial) = oobpredict(rerfdn,X_out{i}(:,:,trial),Y_out{i}(:,trial),'every');

            rerfr = rpclassificationforest(ntrees,X_out{i}(:,:,trial),Y_out{i}(:,trial),'sparsemethod','sparse','Robust',true,'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            err_rerfr.Outlier(:,j,trial) = oobpredict(rerfr,X_out{i}(:,:,trial),Y_out{i}(:,trial),'every');

            rf_rot = rpclassificationforest(ntrees,X_out{i}(:,:,trial),Y_out{i}(:,trial),'RandomForest',true,'rotate',true,'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            err_rf_rot.Outlier(:,j,trial) = oobpredict(rf_rot,X_out{i}(:,:,trial),Y_out{i}(:,trial),'every');

            for k = 1:length(nmixs)
                nmix = nmixs(k);
                fprintf('nmix = %d\n',nmix)

		%Untransformed
                frc = rpclassificationforest(ntrees,X{i}(:,:,trial),Y{i}(:,trial),'sparsemethod','frc','nmix',nmix,'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
                err_frc.Untransformed(:,length(nmixs)*(j-1)+k,trial) = oobpredict(frc,X{i}(:,:,trial),Y{i}(:,trial),'every');

                %Rotated
                frc = rpclassificationforest(ntrees,X_rot{i}(:,:,trial),Y{i}(:,trial),'sparsemethod','frc','nmix',nmix,'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
                err_frc.Rotated(:,length(nmixs)*(j-1)+k,trial) = oobpredict(frc,X_rot{i}(:,:,trial),Y{i}(:,trial),'every');

                %Scaled
                frc = rpclassificationforest(ntrees,X_scale{i}(:,:,trial),Y{i}(:,trial),'sparsemethod','frc','nmix',nmix,'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
                err_frc.Scaled(:,length(nmixs)*(j-1)+k,trial) = oobpredict(frc,X_scale{i}(:,:,trial),Y{i}(:,trial),'every');

                %Affine
                frc = rpclassificationforest(ntrees,X_affine{i}(:,:,trial),Y{i}(:,trial),'sparsemethod','frc','nmix',nmix,'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
                err_frc.Affine(:,length(nmixs)*(j-1)+k,trial) = oobpredict(frc,X_affine{i}(:,:,trial),Y{i}(:,trial),'every');

                %Outlier
                frc = rpclassificationforest(ntrees,X_out{i}(:,:,trial),Y_out{i}(:,trial),'sparsemethod','frc','nmix',nmix,'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
                err_frc.Outlier(:,length(nmixs)*(j-1)+k,trial) = oobpredict(frc,X_out{i}(:,:,trial),Y_out{i}(:,trial),'every');
            end
        end
    end

    if length(mtrys) < 5
        emptyCol = 5 - length(mtrys);
        
        %Untransformed
        sem_rf.Untransformed(:,1:length(mtrys),i) = std(err_rf.Untransformed,[],3)/sqrt(ntrials);
        sem_rerf.Untransformed(:,1:length(mtrys),i) = std(err_rerf.Untransformed,[],3)/sqrt(ntrials);
        sem_rerfr.Untransformed(:,1:length(mtrys),i) = std(err_rerfr.Untransformed,[],3)/sqrt(ntrials);
        sem_rerfdn.Untransformed(:,1:length(mtrys),i) = std(err_rerfdn.Untransformed,[],3)/sqrt(ntrials);
        sem_rf_rot.Untransformed(:,1:length(mtrys),i) = std(err_rf_rot.Untransformed,[],3)/sqrt(ntrials);
        sem_rf.Untransformed(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        sem_rerf.Untransformed(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        sem_rerfr.Untransformed(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        sem_rerfdn.Untransformed(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        sem_rf_rot.Untransformed(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);

        var_rf.Untransformed(:,1:length(mtrys),i) = var(err_rf.Untransformed,0,3);
        var_rerf.Untransformed(:,1:length(mtrys),i) = var(err_rerf.Untransformed,0,3);
        var_rerfr.Untransformed(:,1:length(mtrys),i) = var(err_rerfr.Untransformed,0,3);
        var_rerfdn.Untransformed(:,1:length(mtrys),i) = var(err_rerfdn.Untransformed,0,3);
        var_rf_rot.Untransformed(:,1:length(mtrys),i) = var(err_rf_rot.Untransformed,0,3);
        var_rf.Untransformed(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        var_rerf.Untransformed(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        var_rerfr.Untransformed(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        var_rerfdn.Untransformed(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        var_rf_rot.Untransformed(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);

        mean_err_rf.Untransformed(:,1:length(mtrys),i) = mean(err_rf.Untransformed,3);
        mean_err_rerf.Untransformed(:,1:length(mtrys),i) = mean(err_rerf.Untransformed,3);
        mean_err_rerfr.Untransformed(:,1:length(mtrys),i) = mean(err_rerfr.Untransformed,3);
        mean_err_rerfdn.Untransformed(:,1:length(mtrys),i) = mean(err_rerfdn.Untransformed,3);
        mean_err_rf_rot.Untransformed(:,1:length(mtrys),i) = mean(err_rf_rot.Untransformed,3);
        mean_err_rf.Untransformed(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        mean_err_rerf.Untransformed(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        mean_err_rerfr.Untransformed(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        mean_err_rerfdn.Untransformed(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        mean_err_rf_rot.Untransformed(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        
        %Rotated
        sem_rf.Rotated(:,1:length(mtrys),i) = std(err_rf.Rotated,[],3)/sqrt(ntrials);
        sem_rerf.Rotated(:,1:length(mtrys),i) = std(err_rerf.Rotated,[],3)/sqrt(ntrials);
        sem_rerfr.Rotated(:,1:length(mtrys),i) = std(err_rerfr.Rotated,[],3)/sqrt(ntrials);
        sem_rerfdn.Rotated(:,1:length(mtrys),i) = std(err_rerfdn.Rotated,[],3)/sqrt(ntrials);
        sem_rf_rot.Rotated(:,1:length(mtrys),i) = std(err_rf_rot.Rotated,[],3)/sqrt(ntrials);
        sem_rf.Rotated(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        sem_rerf.Rotated(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        sem_rerfr.Rotated(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        sem_rerfdn.Rotated(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        sem_rf_rot.Rotated(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);

        var_rf.Rotated(:,1:length(mtrys),i) = var(err_rf.Rotated,0,3);
        var_rerf.Rotated(:,1:length(mtrys),i) = var(err_rerf.Rotated,0,3);
        var_rerfr.Rotated(:,1:length(mtrys),i) = var(err_rerfr.Rotated,0,3);
        var_rerfdn.Rotated(:,1:length(mtrys),i) = var(err_rerfdn.Rotated,0,3);
        var_rf_rot.Rotated(:,1:length(mtrys),i) = var(err_rf_rot.Rotated,0,3);
        var_rf.Rotated(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        var_rerf.Rotated(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        var_rerfr.Rotated(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        var_rerfdn.Rotated(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        var_rf_rot.Rotated(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);

        mean_err_rf.Rotated(:,1:length(mtrys),i) = mean(err_rf.Rotated,3);
        mean_err_rerf.Rotated(:,1:length(mtrys),i) = mean(err_rerf.Rotated,3);
        mean_err_rerfr.Rotated(:,1:length(mtrys),i) = mean(err_rerfr.Rotated,3);
        mean_err_rerfdn.Rotated(:,1:length(mtrys),i) = mean(err_rerfdn.Rotated,3);
        mean_err_rf_rot.Rotated(:,1:length(mtrys),i) = mean(err_rf_rot.Rotated,3);
        mean_err_rf.Rotated(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        mean_err_rerf.Rotated(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        mean_err_rerfr.Rotated(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        mean_err_rerfdn.Rotated(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        mean_err_rf_rot.Rotated(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        
        %Scaled
        sem_rf.Scaled(:,1:length(mtrys),i) = std(err_rf.Scaled,[],3)/sqrt(ntrials);
        sem_rerf.Scaled(:,1:length(mtrys),i) = std(err_rerf.Scaled,[],3)/sqrt(ntrials);
        sem_rerfr.Scaled(:,1:length(mtrys),i) = std(err_rerfr.Scaled,[],3)/sqrt(ntrials);
        sem_rerfdn.Scaled(:,1:length(mtrys),i) = std(err_rerfdn.Scaled,[],3)/sqrt(ntrials);
        sem_rf_rot.Scaled(:,1:length(mtrys),i) = std(err_rf_rot.Scaled,[],3)/sqrt(ntrials);
        sem_rf.Scaled(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        sem_rerf.Scaled(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        sem_rerfr.Scaled(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        sem_rerfdn.Scaled(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        sem_rf_rot.Scaled(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);

        var_rf.Scaled(:,1:length(mtrys),i) = var(err_rf.Scaled,0,3);
        var_rerf.Scaled(:,1:length(mtrys),i) = var(err_rerf.Scaled,0,3);
        var_rerfr.Scaled(:,1:length(mtrys),i) = var(err_rerfr.Scaled,0,3);
        var_rerfdn.Scaled(:,1:length(mtrys),i) = var(err_rerfdn.Scaled,0,3);
        var_rf_rot.Scaled(:,1:length(mtrys),i) = var(err_rf_rot.Scaled,0,3);
        var_rf.Scaled(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        var_rerf.Scaled(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        var_rerfr.Scaled(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        var_rerfdn.Scaled(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        var_rf_rot.Scaled(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);

        mean_err_rf.Scaled(:,1:length(mtrys),i) = mean(err_rf.Scaled,3);
        mean_err_rerf.Scaled(:,1:length(mtrys),i) = mean(err_rerf.Scaled,3);
        mean_err_rerfr.Scaled(:,1:length(mtrys),i) = mean(err_rerfr.Scaled,3);
        mean_err_rerfdn.Scaled(:,1:length(mtrys),i) = mean(err_rerfdn.Scaled,3);
        mean_err_rf_rot.Scaled(:,1:length(mtrys),i) = mean(err_rf_rot.Scaled,3);
        mean_err_rf.Scaled(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        mean_err_rerf.Scaled(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        mean_err_rerfr.Scaled(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        mean_err_rerfdn.Scaled(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        mean_err_rf_rot.Scaled(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        
        %Affine
        sem_rf.Affine(:,1:length(mtrys),i) = std(err_rf.Affine,[],3)/sqrt(ntrials);
        sem_rerf.Affine(:,1:length(mtrys),i) = std(err_rerf.Affine,[],3)/sqrt(ntrials);
        sem_rerfr.Affine(:,1:length(mtrys),i) = std(err_rerfr.Affine,[],3)/sqrt(ntrials);
        sem_rerfdn.Affine(:,1:length(mtrys),i) = std(err_rerfdn.Affine,[],3)/sqrt(ntrials);
        sem_rf_rot.Affine(:,1:length(mtrys),i) = std(err_rf_rot.Affine,[],3)/sqrt(ntrials);
        sem_rf.Affine(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        sem_rerf.Affine(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        sem_rerfr.Affine(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        sem_rerfdn.Affine(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        sem_rf_rot.Affine(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);

        var_rf.Affine(:,1:length(mtrys),i) = var(err_rf.Affine,0,3);
        var_rerf.Affine(:,1:length(mtrys),i) = var(err_rerf.Affine,0,3);
        var_rerfr.Affine(:,1:length(mtrys),i) = var(err_rerfr.Affine,0,3);
        var_rerfdn.Affine(:,1:length(mtrys),i) = var(err_rerfdn.Affine,0,3);
        var_rf_rot.Affine(:,1:length(mtrys),i) = var(err_rf_rot.Affine,0,3);
        var_rf.Affine(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        var_rerf.Affine(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        var_rerfr.Affine(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        var_rerfdn.Affine(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        var_rf_rot.Affine(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);

        mean_err_rf.Affine(:,1:length(mtrys),i) = mean(err_rf.Affine,3);
        mean_err_rerf.Affine(:,1:length(mtrys),i) = mean(err_rerf.Affine,3);
        mean_err_rerfr.Affine(:,1:length(mtrys),i) = mean(err_rerfr.Affine,3);
        mean_err_rerfdn.Affine(:,1:length(mtrys),i) = mean(err_rerfdn.Affine,3);
        mean_err_rf_rot.Affine(:,1:length(mtrys),i) = mean(err_rf_rot.Affine,3);
        mean_err_rf.Affine(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        mean_err_rerf.Affine(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        mean_err_rerfr.Affine(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        mean_err_rerfdn.Affine(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        mean_err_rf_rot.Affine(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        
        %Outlier
        sem_rf.Outlier(:,1:length(mtrys),i) = std(err_rf.Outlier,[],3)/sqrt(ntrials);
        sem_rerf.Outlier(:,1:length(mtrys),i) = std(err_rerf.Outlier,[],3)/sqrt(ntrials);
        sem_rerfr.Outlier(:,1:length(mtrys),i) = std(err_rerfr.Outlier,[],3)/sqrt(ntrials);
        sem_rerfdn.Outlier(:,1:length(mtrys),i) = std(err_rerfdn.Outlier,[],3)/sqrt(ntrials);
        sem_rf_rot.Outlier(:,1:length(mtrys),i) = std(err_rf_rot.Outlier,[],3)/sqrt(ntrials);
        sem_rf.Outlier(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        sem_rerf.Outlier(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        sem_rerfr.Outlier(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        sem_rerfdn.Outlier(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        sem_rf_rot.Outlier(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);

        var_rf.Outlier(:,1:length(mtrys),i) = var(err_rf.Outlier,0,3);
        var_rerf.Outlier(:,1:length(mtrys),i) = var(err_rerf.Outlier,0,3);
        var_rerfr.Outlier(:,1:length(mtrys),i) = var(err_rerfr.Outlier,0,3);
        var_rerfdn.Outlier(:,1:length(mtrys),i) = var(err_rerfdn.Outlier,0,3);
        var_rf_rot.Outlier(:,1:length(mtrys),i) = var(err_rf_rot.Outlier,0,3);
        var_rf.Outlier(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        var_rerf.Outlier(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        var_rerfr.Outlier(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        var_rerfdn.Outlier(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        var_rf_rot.Outlier(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);

        mean_err_rf.Outlier(:,1:length(mtrys),i) = mean(err_rf.Outlier,3);
        mean_err_rerf.Outlier(:,1:length(mtrys),i) = mean(err_rerf.Outlier,3);
        mean_err_rerfr.Outlier(:,1:length(mtrys),i) = mean(err_rerfr.Outlier,3);
        mean_err_rerfdn.Outlier(:,1:length(mtrys),i) = mean(err_rerfdn.Outlier,3);
        mean_err_rf_rot.Outlier(:,1:length(mtrys),i) = mean(err_rf_rot.Outlier,3);
        mean_err_rf.Outlier(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        mean_err_rerf.Outlier(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        mean_err_rerfr.Outlier(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        mean_err_rerfdn.Outlier(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
        mean_err_rf_rot.Outlier(:,length(mtrys)+1:5,i) = NaN(ntrees,emptyCol);
    else
        %Untransformed
        sem_rf.Untransformed(:,:,i) = std(err_rf.Untransformed,[],3)/sqrt(ntrials);
        sem_rerf.Untransformed(:,:,i) = std(err_rerf.Untransformed,[],3)/sqrt(ntrials);
        sem_rerfr.Untransformed(:,:,i) = std(err_rerfr.Untransformed,[],3)/sqrt(ntrials);
        sem_rerfdn.Untransformed(:,:,i) = std(err_rerfdn.Untransformed,[],3)/sqrt(ntrials);
        sem_rf_rot.Untransformed(:,:,i) = std(err_rf_rot.Untransformed,[],3)/sqrt(ntrials);

        var_rf.Untransformed(:,:,i) = var(err_rf.Untransformed,0,3);
        var_rerf.Untransformed(:,:,i) = var(err_rerf.Untransformed,0,3);
        var_rerfr.Untransformed(:,:,i) = var(err_rerfr.Untransformed,0,3);
        var_rerfdn.Untransformed(:,:,i) = var(err_rerfdn.Untransformed,0,3);
        var_rf_rot.Untransformed(:,:,i) = var(err_rf_rot.Untransformed,0,3);

        mean_err_rf.Untransformed(:,:,i) = mean(err_rf.Untransformed,3);
        mean_err_rerf.Untransformed(:,:,i) = mean(err_rerf.Untransformed,3);
        mean_err_rerfr.Untransformed(:,:,i) = mean(err_rerfr.Untransformed,3);
        mean_err_rerfdn.Untransformed(:,:,i) = mean(err_rerfdn.Untransformed,3);
        mean_err_rf_rot.Untransformed(:,:,i) = mean(err_rf_rot.Untransformed,3);
        
        %Rotated
        sem_rf.Rotated(:,:,i) = std(err_rf.Rotated,[],3)/sqrt(ntrials);
        sem_rerf.Rotated(:,:,i) = std(err_rerf.Rotated,[],3)/sqrt(ntrials);
        sem_rerfr.Rotated(:,:,i) = std(err_rerfr.Rotated,[],3)/sqrt(ntrials);
        sem_rerfdn.Rotated(:,:,i) = std(err_rerfdn.Rotated,[],3)/sqrt(ntrials);
        sem_rf_rot.Rotated(:,:,i) = std(err_rf_rot.Rotated,[],3)/sqrt(ntrials);

        var_rf.Rotated(:,:,i) = var(err_rf.Rotated,0,3);
        var_rerf.Rotated(:,:,i) = var(err_rerf.Rotated,0,3);
        var_rerfr.Rotated(:,:,i) = var(err_rerfr.Rotated,0,3);
        var_rerfdn.Rotated(:,:,i) = var(err_rerfdn.Rotated,0,3);
        var_rf_rot.Rotated(:,:,i) = var(err_rf_rot.Rotated,0,3);

        mean_err_rf.Rotated(:,:,i) = mean(err_rf.Rotated,3);
        mean_err_rerf.Rotated(:,:,i) = mean(err_rerf.Rotated,3);
        mean_err_rerfr.Rotated(:,:,i) = mean(err_rerfr.Rotated,3);
        mean_err_rerfdn.Rotated(:,:,i) = mean(err_rerfdn.Rotated,3);
        mean_err_rf_rot.Rotated(:,:,i) = mean(err_rf_rot.Rotated,3);
        
        %Scaled
        sem_rf.Scaled(:,:,i) = std(err_rf.Scaled,[],3)/sqrt(ntrials);
        sem_rerf.Scaled(:,:,i) = std(err_rerf.Scaled,[],3)/sqrt(ntrials);
        sem_rerfr.Scaled(:,:,i) = std(err_rerfr.Scaled,[],3)/sqrt(ntrials);
        sem_rerfdn.Scaled(:,:,i) = std(err_rerfdn.Scaled,[],3)/sqrt(ntrials);
        sem_rf_rot.Scaled(:,:,i) = std(err_rf_rot.Scaled,[],3)/sqrt(ntrials);

        var_rf.Scaled(:,:,i) = var(err_rf.Scaled,0,3);
        var_rerf.Scaled(:,:,i) = var(err_rerf.Scaled,0,3);
        var_rerfr.Scaled(:,:,i) = var(err_rerfr.Scaled,0,3);
        var_rerfdn.Scaled(:,:,i) = var(err_rerfdn.Scaled,0,3);
        var_rf_rot.Scaled(:,:,i) = var(err_rf_rot.Scaled,0,3);

        mean_err_rf.Scaled(:,:,i) = mean(err_rf.Scaled,3);
        mean_err_rerf.Scaled(:,:,i) = mean(err_rerf.Scaled,3);
        mean_err_rerfr.Scaled(:,:,i) = mean(err_rerfr.Scaled,3);
        mean_err_rerfdn.Scaled(:,:,i) = mean(err_rerfdn.Scaled,3);
        mean_err_rf_rot.Scaled(:,:,i) = mean(err_rf_rot.Scaled,3);
        
        %Affine
        sem_rf.Affine(:,:,i) = std(err_rf.Affine,[],3)/sqrt(ntrials);
        sem_rerf.Affine(:,:,i) = std(err_rerf.Affine,[],3)/sqrt(ntrials);
        sem_rerfr.Affine(:,:,i) = std(err_rerfr.Affine,[],3)/sqrt(ntrials);
        sem_rerfdn.Affine(:,:,i) = std(err_rerfdn.Affine,[],3)/sqrt(ntrials);
        sem_rf_rot.Affine(:,:,i) = std(err_rf_rot.Affine,[],3)/sqrt(ntrials);

        var_rf.Affine(:,:,i) = var(err_rf.Affine,0,3);
        var_rerf.Affine(:,:,i) = var(err_rerf.Affine,0,3);
        var_rerfr.Affine(:,:,i) = var(err_rerfr.Affine,0,3);
        var_rerfdn.Affine(:,:,i) = var(err_rerfdn.Affine,0,3);
        var_rf_rot.Affine(:,:,i) = var(err_rf_rot.Affine,0,3);

        mean_err_rf.Affine(:,:,i) = mean(err_rf.Affine,3);
        mean_err_rerf.Affine(:,:,i) = mean(err_rerf.Affine,3);
        mean_err_rerfr.Affine(:,:,i) = mean(err_rerfr.Affine,3);
        mean_err_rerfdn.Affine(:,:,i) = mean(err_rerfdn.Affine,3);
        mean_err_rf_rot.Affine(:,:,i) = mean(err_rf_rot.Affine,3);
        
        %Outlier
        sem_rf.Outlier(:,:,i) = std(err_rf.Outlier,[],3)/sqrt(ntrials);
        sem_rerf.Outlier(:,:,i) = std(err_rerf.Outlier,[],3)/sqrt(ntrials);
        sem_rerfr.Outlier(:,:,i) = std(err_rerfr.Outlier,[],3)/sqrt(ntrials);
        sem_rerfdn.Outlier(:,:,i) = std(err_rerfdn.Outlier,[],3)/sqrt(ntrials);
        sem_rf_rot.Outlier(:,:,i) = std(err_rf_rot.Outlier,[],3)/sqrt(ntrials);

        var_rf.Outlier(:,:,i) = var(err_rf.Outlier,0,3);
        var_rerf.Outlier(:,:,i) = var(err_rerf.Outlier,0,3);
        var_rerfr.Outlier(:,:,i) = var(err_rerfr.Outlier,0,3);
        var_rerfdn.Outlier(:,:,i) = var(err_rerfdn.Outlier,0,3);
        var_rf_rot.Outlier(:,:,i) = var(err_rf_rot.Outlier,0,3);

        mean_err_rf.Outlier(:,:,i) = mean(err_rf.Outlier,3);
        mean_err_rerf.Outlier(:,:,i) = mean(err_rerf.Outlier,3);
        mean_err_rerfr.Outlier(:,:,i) = mean(err_rerfr.Outlier,3);
        mean_err_rerfdn.Outlier(:,:,i) = mean(err_rerfdn.Outlier,3);
        mean_err_rf_rot.Outlier(:,:,i) = mean(err_rf_rot.Outlier,3);
    end

    if d < 6
        emptyCol = 25 - length(mtrys)*length(nmixs);
        
        %Untransformed
        sem_frc.Untransformed(:,1:length(mtrys)*length(nmixs),i) = std(err_frc.Untransformed,[],3)/sqrt(ntrials);
        sem_frc.Untransformed(:,length(mtrys)*length(nmixs)+1:25,i) = NaN(ntrees,emptyCol);

        var_frc.Untransformed(:,1:length(mtrys)*length(nmixs),i) = var(err_frc.Untransformed,0,3);
        var_frc.Untransformed(:,length(mtrys)*length(nmixs)+1:25,i) = NaN(ntrees,emptyCol);

        mean_err_frc.Untransformed(:,1:length(mtrys)*length(nmixs),i) = mean(err_frc.Untransformed,3);
        mean_err_frc.Untransformed(:,length(mtrys)*length(nmixs)+1:25,i) = NaN(ntrees,emptyCol);
        
        %Rotated
        sem_frc.Rotated(:,1:length(mtrys)*length(nmixs),i) = std(err_frc.Rotated,[],3)/sqrt(ntrials);
        sem_frc.Rotated(:,length(mtrys)*length(nmixs)+1:25,i) = NaN(ntrees,emptyCol);

        var_frc.Rotated(:,1:length(mtrys)*length(nmixs),i) = var(err_frc.Rotated,0,3);
        var_frc.Rotated(:,length(mtrys)*length(nmixs)+1:25,i) = NaN(ntrees,emptyCol);

        mean_err_frc.Rotated(:,1:length(mtrys)*length(nmixs),i) = mean(err_frc.Rotated,3);
        mean_err_frc.Rotated(:,length(mtrys)*length(nmixs)+1:25,i) = NaN(ntrees,emptyCol);
        
        %Scaled
        sem_frc.Scaled(:,1:length(mtrys)*length(nmixs),i) = std(err_frc.Scaled,[],3)/sqrt(ntrials);
        sem_frc.Scaled(:,length(mtrys)*length(nmixs)+1:25,i) = NaN(ntrees,emptyCol);

        var_frc.Scaled(:,1:length(mtrys)*length(nmixs),i) = var(err_frc.Scaled,0,3);
        var_frc.Scaled(:,length(mtrys)*length(nmixs)+1:25,i) = NaN(ntrees,emptyCol);

        mean_err_frc.Scaled(:,1:length(mtrys)*length(nmixs),i) = mean(err_frc.Scaled,3);
        mean_err_frc.Scaled(:,length(mtrys)*length(nmixs)+1:25,i) = NaN(ntrees,emptyCol);
        
        %Affine
        sem_frc.Affine(:,1:length(mtrys)*length(nmixs),i) = std(err_frc.Affine,[],3)/sqrt(ntrials);
        sem_frc.Affine(:,length(mtrys)*length(nmixs)+1:25,i) = NaN(ntrees,emptyCol);

        var_frc.Affine(:,1:length(mtrys)*length(nmixs),i) = var(err_frc.Affine,0,3);
        var_frc.Affine(:,length(mtrys)*length(nmixs)+1:25,i) = NaN(ntrees,emptyCol);

        mean_err_frc.Affine(:,1:length(mtrys)*length(nmixs),i) = mean(err_frc.Affine,3);
        mean_err_frc.Affine(:,length(mtrys)*length(nmixs)+1:25,i) = NaN(ntrees,emptyCol);
        
        %Outlier
        sem_frc.Outlier(:,1:length(mtrys)*length(nmixs),i) = std(err_frc.Outlier,[],3)/sqrt(ntrials);
        sem_frc.Outlier(:,length(mtrys)*length(nmixs)+1:25,i) = NaN(ntrees,emptyCol);

        var_frc.Outlier(:,1:length(mtrys)*length(nmixs),i) = var(err_frc.Outlier,0,3);
        var_frc.Outlier(:,length(mtrys)*length(nmixs)+1:25,i) = NaN(ntrees,emptyCol);

        mean_err_frc.Outlier(:,1:length(mtrys)*length(nmixs),i) = mean(err_frc.Outlier,3);
        mean_err_frc.Outlier(:,length(mtrys)*length(nmixs)+1:25,i) = NaN(ntrees,emptyCol);
    else
        %Untransformed
        sem_frc.Untransformed(:,:,i) = std(err_frc.Untransformed,[],3)/sqrt(ntrials);

        var_frc.Untransformed(:,:,i) = var(err_frc.Untransformed,0,3);

        mean_err_frc.Untransformed(:,:,i) = mean(err_frc.Untransformed,3);
        
        %Rotated
        sem_frc.Rotated(:,:,i) = std(err_frc.Rotated,[],3)/sqrt(ntrials);

        var_frc.Rotated(:,:,i) = var(err_frc.Rotated,0,3);

        mean_err_frc.Rotated(:,:,i) = mean(err_frc.Rotated,3);
        
        %Scaled
        sem_frc.Scaled(:,:,i) = std(err_frc.Scaled,[],3)/sqrt(ntrials);

        var_frc.Scaled(:,:,i) = var(err_frc.Scaled,0,3);

        mean_err_frc.Scaled(:,:,i) = mean(err_frc.Scaled,3);
        
        %Affine
        sem_frc.Affine(:,:,i) = std(err_frc.Affine,[],3)/sqrt(ntrials);

        var_frc.Affine(:,:,i) = var(err_frc.Affine,0,3);

        mean_err_frc.Affine(:,:,i) = mean(err_frc.Affine,3);
        
        %Outlier
        sem_frc.Outlier(:,:,i) = std(err_frc.Outlier,[],3)/sqrt(ntrials);

        var_frc.Outlier(:,:,i) = var(err_frc.Outlier,0,3);

        mean_err_frc.Outlier(:,:,i) = mean(err_frc.Outlier,3);
    end

end

save([rerfPath 'RandomerForest/Results/Sparse_parity_transformations2.mat'],...
    'dims','sem_rf','var_rf','mean_err_rf','sem_rerf','var_rerf','mean_err_rerf',...
    'sem_rerfdn','var_rerfdn','mean_err_rerfdn','sem_rerfr','var_rerfr','mean_err_rerfr',...
    'sem_rf_rot','var_rf_rot','mean_err_rf_rot','sem_frc','var_frc','mean_err_frc')
