close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

rng(1);

n = 100;
ntrees = 1000;
ntrials = 10;
NWorkers = 2;
dims = [2 10 50 100 500];

for j = 1:length(dims)
    
    d = dims(j);
    fprintf('dimension %d\n',d)
    if d <= 5
        mtrys = 1:d;
    else
        mtrys = ceil(d.^[0 1/4 1/2 3/4 1]);
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
    
    Class = [0;1];

    for trial = 1:ntrials

        fprintf('trial %d\n',trial)

        d_idx = 1:d;
        mu1 = 1./sqrt(d_idx);
        mu0 = -1*mu1;
        Mu = cat(1,mu0,mu1);
        Sigma = ones(1,d);
        obj = gmdistribution(Mu,Sigma);
        [X,idx] = random(obj,n);
        Ystr = cellstr(num2str(Class(idx)));
        
        R = random_rotation(d);
        S = random_scaling(n,d,0,10);
        Sigma_outlier = 16*Sigma;
        X_rot = X*R;
        X_scale = X.*S;
        X_affine = (X*R).*S;
        outlier_model = gmdistribution(Mu,Sigma_outlier);
        [X_out,idx_out] = random(outlier_model,0.05*n);
        X_out = cat(1,X,X_out);
        Y_out = cellstr(num2str(Class(idx_out)));
        Y_out = cat(1,Ystr,Y_out);

        i = 1;

        for mtry = mtrys

            fprintf('mtry = %d\n',mtry)
            
            %Untransformed
            rf = rpclassificationforest(ntrees,X,Ystr,'RandomForest',true,'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            err_rf.Untransformed(:,i,trial) = oobpredict(rf,X,Ystr,'every');

            rerf = rpclassificationforest(ntrees,X,Ystr,'sparsemethod','sparse','nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            err_rerf.Untransformed(:,i,trial) = oobpredict(rerf,X,Ystr,'every');
            
            rerfr = rpclassificationforest(ntrees,X,Ystr,'sparsemethod','sparse','nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true,'Robust',true);
            err_rerfr.Untransformed(:,i,trial) = oobpredict(rerfr,X,Ystr,'every');

            rerfdn = rpclassificationforest(ntrees,X,Ystr,'sparsemethod','sparse','mdiff','node','nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            err_rerfdn.Untransformed(:,i,trial) = oobpredict(rerfdn,X,Ystr,'every');
            
            rf_rot = rpclassificationforest(ntrees,X,Ystr,'RandomForest',true,'rotate',true,'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            err_rf_rot.Untransformed(:,i,trial) = oobpredict(rf_rot,X,Ystr,'every');
            
            %Rotated
            rf = rpclassificationforest(ntrees,X_rot,Ystr,'RandomForest',true,'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            err_rf.Rotated(:,i,trial) = oobpredict(rf,X_rot,Ystr,'every');

            rerf = rpclassificationforest(ntrees,X_rot,Ystr,'sparsemethod','sparse','nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            err_rerf.Rotated(:,i,trial) = oobpredict(rerf,X_rot,Ystr,'every');
            
            rerfr = rpclassificationforest(ntrees,X_rot,Ystr,'sparsemethod','sparse','nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true,'Robust',true);
            err_rerfr.Rotated(:,i,trial) = oobpredict(rerfr,X_rot,Ystr,'every');

            rerfdn = rpclassificationforest(ntrees,X_rot,Ystr,'sparsemethod','sparse','mdiff','node','nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            err_rerfdn.Rotated(:,i,trial) = oobpredict(rerfdn,X_rot,Ystr,'every');
            
            rf_rot = rpclassificationforest(ntrees,X_rot,Ystr,'RandomForest',true,'rotate',true,'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            err_rf_rot.Rotated(:,i,trial) = oobpredict(rf_rot,X_rot,Ystr,'every');
            
            %Scaled
            rf = rpclassificationforest(ntrees,X_scale,Ystr,'RandomForest',true,'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            err_rf.Scaled(:,i,trial) = oobpredict(rf,X_scale,Ystr,'every');

            rerf = rpclassificationforest(ntrees,X_scale,Ystr,'sparsemethod','sparse','nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            err_rerf.Scaled(:,i,trial) = oobpredict(rerf,X_scale,Ystr,'every');
            
            rerfr = rpclassificationforest(ntrees,X_scale,Ystr,'sparsemethod','sparse','nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true,'Robust',true);
            err_rerfr.Scaled(:,i,trial) = oobpredict(rerfr,X_scale,Ystr,'every');

            rerfdn = rpclassificationforest(ntrees,X_scale,Ystr,'sparsemethod','sparse','mdiff','node','nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            err_rerfdn.Scaled(:,i,trial) = oobpredict(rerfdn,X_scale,Ystr,'every');
            
            rf_rot = rpclassificationforest(ntrees,X_scale,Ystr,'RandomForest',true,'rotate',true,'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            err_rf_rot.Scaled(:,i,trial) = oobpredict(rf_rot,X_scale,Ystr,'every');
            
            %Affine
            rf = rpclassificationforest(ntrees,X_affine,Ystr,'RandomForest',true,'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            err_rf.Affine(:,i,trial) = oobpredict(rf,X_affine,Ystr,'every');

            rerf = rpclassificationforest(ntrees,X_affine,Ystr,'sparsemethod','sparse','nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            err_rerf.Affine(:,i,trial) = oobpredict(rerf,X_affine,Ystr,'every');
            
            rerfr = rpclassificationforest(ntrees,X_affine,Ystr,'sparsemethod','sparse','nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true,'Robust',true);
            err_rerfr.Affine(:,i,trial) = oobpredict(rerfr,X_affine,Ystr,'every');

            rerfdn = rpclassificationforest(ntrees,X_affine,Ystr,'sparsemethod','sparse','mdiff','node','nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            err_rerfdn.Affine(:,i,trial) = oobpredict(rerfdn,X_affine,Ystr,'every');
            
            rf_rot = rpclassificationforest(ntrees,X_affine,Ystr,'RandomForest',true,'rotate',true,'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            err_rf_rot.Affine(:,i,trial) = oobpredict(rf_rot,X_affine,Ystr,'every');
            
            %Outlier
            rf = rpclassificationforest(ntrees,X_out,Y_out,'RandomForest',true,'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            err_rf.Outlier(:,i,trial) = oobpredict(rf,X_out,Y_out,'every');

            rerf = rpclassificationforest(ntrees,X_out,Y_out,'sparsemethod','sparse','nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            err_rerf.Outlier(:,i,trial) = oobpredict(rerf,X_out,Y_out,'every');
            
            rerfr = rpclassificationforest(ntrees,X_out,Y_out,'sparsemethod','sparse','nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true,'Robust',true);
            err_rerfr.Outlier(:,i,trial) = oobpredict(rerfr,X_out,Y_out,'every');

            rerfdn = rpclassificationforest(ntrees,X_out,Y_out,'sparsemethod','sparse','mdiff','node','nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            err_rerfdn.Outlier(:,i,trial) = oobpredict(rerfdn,X_out,Y_out,'every');
            
            rf_rot = rpclassificationforest(ntrees,X_out,Y_out,'RandomForest',true,'rotate',true,'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            err_rf_rot.Outlier(:,i,trial) = oobpredict(rf_rot,X_out,Y_out,'every');

            i = i + 1;
        end
    end

    save(sprintf([rerfPath 'RandomerForest/Results/Trunk_transformations_n%d_d%d.mat'],n,d),'err_rf','err_rerf','err_rerfr','err_rerfdn','err_rf_rot')

    if length(mtrys) < 5
        emptyCol = 5 - length(mtrys);
        
        %Untransformed
        sem_rf.Untransformed(:,1:length(mtrys),j) = std(err_rf.Untransformed,[],3)/sqrt(ntrials);
        sem_rerf.Untransformed(:,1:length(mtrys),j) = std(err_rerf.Untransformed,[],3)/sqrt(ntrials);
        sem_rerfr.Untransformed(:,1:length(mtrys),j) = std(err_rerfr.Untransformed,[],3)/sqrt(ntrials);
        sem_rerfdn.Untransformed(:,1:length(mtrys),j) = std(err_rerfdn.Untransformed,[],3)/sqrt(ntrials);
        sem_rf_rot.Untransformed(:,1:length(mtrys),j) = std(err_rf_rot.Untransformed,[],3)/sqrt(ntrials);
        sem_rf.Untransformed(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);
        sem_rerf.Untransformed(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);
        sem_rerfr.Untransformed(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);
        sem_rerfdn.Untransformed(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);
        sem_rf_rot.Untransformed(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);

        var_rf.Untransformed(:,1:length(mtrys),j) = var(err_rf.Untransformed,0,3);
        var_rerf.Untransformed(:,1:length(mtrys),j) = var(err_rerf.Untransformed,0,3);
        var_rerfr.Untransformed(:,1:length(mtrys),j) = var(err_rerfr.Untransformed,0,3);
        var_rerfdn.Untransformed(:,1:length(mtrys),j) = var(err_rerfdn.Untransformed,0,3);
        var_rf_rot.Untransformed(:,1:length(mtrys),j) = var(err_rf_rot.Untransformed,0,3);
        var_rf.Untransformed(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);
        var_rerf.Untransformed(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);
        var_rerfr.Untransformed(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);
        var_rerfdn.Untransformed(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);
        var_rf_rot.Untransformed(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);

        mean_err_rf.Untransformed(:,1:length(mtrys),j) = mean(err_rf.Untransformed,3);
        mean_err_rerf.Untransformed(:,1:length(mtrys),j) = mean(err_rerf.Untransformed,3);
        mean_err_rerfr.Untransformed(:,1:length(mtrys),j) = mean(err_rerfr.Untransformed,3);
        mean_err_rerfdn.Untransformed(:,1:length(mtrys),j) = mean(err_rerfdn.Untransformed,3);
        mean_err_rf_rot.Untransformed(:,1:length(mtrys),j) = mean(err_rf_rot.Untransformed,3);
        mean_err_rf.Untransformed(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);
        mean_err_rerf.Untransformed(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);
        mean_err_rerfr.Untransformed(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);
        mean_err_rerfdn.Untransformed(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);
        mean_err_rf_rot.Untransformed(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);
        
        %Rotated
        sem_rf.Rotated(:,1:length(mtrys),j) = std(err_rf.Rotated,[],3)/sqrt(ntrials);
        sem_rerf.Rotated(:,1:length(mtrys),j) = std(err_rerf.Rotated,[],3)/sqrt(ntrials);
        sem_rerfr.Rotated(:,1:length(mtrys),j) = std(err_rerfr.Rotated,[],3)/sqrt(ntrials);
        sem_rerfdn.Rotated(:,1:length(mtrys),j) = std(err_rerfdn.Rotated,[],3)/sqrt(ntrials);
        sem_rf_rot.Rotated(:,1:length(mtrys),j) = std(err_rf_rot.Rotated,[],3)/sqrt(ntrials);
        sem_rf.Rotated(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);
        sem_rerf.Rotated(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);
        sem_rerfr.Rotated(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);
        sem_rerfdn.Rotated(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);
        sem_rf_rot.Rotated(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);

        var_rf.Rotated(:,1:length(mtrys),j) = var(err_rf.Rotated,0,3);
        var_rerf.Rotated(:,1:length(mtrys),j) = var(err_rerf.Rotated,0,3);
        var_rerfr.Rotated(:,1:length(mtrys),j) = var(err_rerfr.Rotated,0,3);
        var_rerfdn.Rotated(:,1:length(mtrys),j) = var(err_rerfdn.Rotated,0,3);
        var_rf_rot.Rotated(:,1:length(mtrys),j) = var(err_rf_rot.Rotated,0,3);
        var_rf.Rotated(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);
        var_rerf.Rotated(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);
        var_rerfr.Rotated(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);
        var_rerfdn.Rotated(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);
        var_rf_rot.Rotated(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);

        mean_err_rf.Rotated(:,1:length(mtrys),j) = mean(err_rf.Rotated,3);
        mean_err_rerf.Rotated(:,1:length(mtrys),j) = mean(err_rerf.Rotated,3);
        mean_err_rerfr.Rotated(:,1:length(mtrys),j) = mean(err_rerfr.Rotated,3);
        mean_err_rerfdn.Rotated(:,1:length(mtrys),j) = mean(err_rerfdn.Rotated,3);
        mean_err_rf_rot.Rotated(:,1:length(mtrys),j) = mean(err_rf_rot.Rotated,3);
        mean_err_rf.Rotated(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);
        mean_err_rerf.Rotated(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);
        mean_err_rerfr.Rotated(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);
        mean_err_rerfdn.Rotated(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);
        mean_err_rf_rot.Rotated(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);
        
        %Scaled
        sem_rf.Scaled(:,1:length(mtrys),j) = std(err_rf.Scaled,[],3)/sqrt(ntrials);
        sem_rerf.Scaled(:,1:length(mtrys),j) = std(err_rerf.Scaled,[],3)/sqrt(ntrials);
        sem_rerfr.Scaled(:,1:length(mtrys),j) = std(err_rerfr.Scaled,[],3)/sqrt(ntrials);
        sem_rerfdn.Scaled(:,1:length(mtrys),j) = std(err_rerfdn.Scaled,[],3)/sqrt(ntrials);
        sem_rf_rot.Scaled(:,1:length(mtrys),j) = std(err_rf_rot.Scaled,[],3)/sqrt(ntrials);
        sem_rf.Scaled(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);
        sem_rerf.Scaled(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);
        sem_rerfr.Scaled(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);
        sem_rerfdn.Scaled(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);
        sem_rf_rot.Scaled(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);

        var_rf.Scaled(:,1:length(mtrys),j) = var(err_rf.Scaled,0,3);
        var_rerf.Scaled(:,1:length(mtrys),j) = var(err_rerf.Scaled,0,3);
        var_rerfr.Scaled(:,1:length(mtrys),j) = var(err_rerfr.Scaled,0,3);
        var_rerfdn.Scaled(:,1:length(mtrys),j) = var(err_rerfdn.Scaled,0,3);
        var_rf_rot.Scaled(:,1:length(mtrys),j) = var(err_rf_rot.Scaled,0,3);
        var_rf.Scaled(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);
        var_rerf.Scaled(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);
        var_rerfr.Scaled(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);
        var_rerfdn.Scaled(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);
        var_rf_rot.Scaled(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);

        mean_err_rf.Scaled(:,1:length(mtrys),j) = mean(err_rf.Scaled,3);
        mean_err_rerf.Scaled(:,1:length(mtrys),j) = mean(err_rerf.Scaled,3);
        mean_err_rerfr.Scaled(:,1:length(mtrys),j) = mean(err_rerfr.Scaled,3);
        mean_err_rerfdn.Scaled(:,1:length(mtrys),j) = mean(err_rerfdn.Scaled,3);
        mean_err_rf_rot.Scaled(:,1:length(mtrys),j) = mean(err_rf_rot.Scaled,3);
        mean_err_rf.Scaled(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);
        mean_err_rerf.Scaled(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);
        mean_err_rerfr.Scaled(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);
        mean_err_rerfdn.Scaled(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);
        mean_err_rf_rot.Scaled(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);
        
        %Affine
        sem_rf.Affine(:,1:length(mtrys),j) = std(err_rf.Affine,[],3)/sqrt(ntrials);
        sem_rerf.Affine(:,1:length(mtrys),j) = std(err_rerf.Affine,[],3)/sqrt(ntrials);
        sem_rerfr.Affine(:,1:length(mtrys),j) = std(err_rerfr.Affine,[],3)/sqrt(ntrials);
        sem_rerfdn.Affine(:,1:length(mtrys),j) = std(err_rerfdn.Affine,[],3)/sqrt(ntrials);
        sem_rf_rot.Affine(:,1:length(mtrys),j) = std(err_rf_rot.Affine,[],3)/sqrt(ntrials);
        sem_rf.Affine(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);
        sem_rerf.Affine(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);
        sem_rerfr.Affine(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);
        sem_rerfdn.Affine(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);
        sem_rf_rot.Affine(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);

        var_rf.Affine(:,1:length(mtrys),j) = var(err_rf.Affine,0,3);
        var_rerf.Affine(:,1:length(mtrys),j) = var(err_rerf.Affine,0,3);
        var_rerfr.Affine(:,1:length(mtrys),j) = var(err_rerfr.Affine,0,3);
        var_rerfdn.Affine(:,1:length(mtrys),j) = var(err_rerfdn.Affine,0,3);
        var_rf_rot.Affine(:,1:length(mtrys),j) = var(err_rf_rot.Affine,0,3);
        var_rf.Affine(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);
        var_rerf.Affine(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);
        var_rerfr.Affine(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);
        var_rerfdn.Affine(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);
        var_rf_rot.Affine(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);

        mean_err_rf.Affine(:,1:length(mtrys),j) = mean(err_rf.Affine,3);
        mean_err_rerf.Affine(:,1:length(mtrys),j) = mean(err_rerf.Affine,3);
        mean_err_rerfr.Affine(:,1:length(mtrys),j) = mean(err_rerfr.Affine,3);
        mean_err_rerfdn.Affine(:,1:length(mtrys),j) = mean(err_rerfdn.Affine,3);
        mean_err_rf_rot.Affine(:,1:length(mtrys),j) = mean(err_rf_rot.Affine,3);
        mean_err_rf.Affine(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);
        mean_err_rerf.Affine(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);
        mean_err_rerfr.Affine(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);
        mean_err_rerfdn.Affine(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);
        mean_err_rf_rot.Affine(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);
        
        %Outlier
        sem_rf.Outlier(:,1:length(mtrys),j) = std(err_rf.Outlier,[],3)/sqrt(ntrials);
        sem_rerf.Outlier(:,1:length(mtrys),j) = std(err_rerf.Outlier,[],3)/sqrt(ntrials);
        sem_rerfr.Outlier(:,1:length(mtrys),j) = std(err_rerfr.Outlier,[],3)/sqrt(ntrials);
        sem_rerfdn.Outlier(:,1:length(mtrys),j) = std(err_rerfdn.Outlier,[],3)/sqrt(ntrials);
        sem_rf_rot.Outlier(:,1:length(mtrys),j) = std(err_rf_rot.Outlier,[],3)/sqrt(ntrials);
        sem_rf.Outlier(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);
        sem_rerf.Outlier(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);
        sem_rerfr.Outlier(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);
        sem_rerfdn.Outlier(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);
        sem_rf_rot.Outlier(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);

        var_rf.Outlier(:,1:length(mtrys),j) = var(err_rf.Outlier,0,3);
        var_rerf.Outlier(:,1:length(mtrys),j) = var(err_rerf.Outlier,0,3);
        var_rerfr.Outlier(:,1:length(mtrys),j) = var(err_rerfr.Outlier,0,3);
        var_rerfdn.Outlier(:,1:length(mtrys),j) = var(err_rerfdn.Outlier,0,3);
        var_rf_rot.Outlier(:,1:length(mtrys),j) = var(err_rf_rot.Outlier,0,3);
        var_rf.Outlier(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);
        var_rerf.Outlier(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);
        var_rerfr.Outlier(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);
        var_rerfdn.Outlier(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);
        var_rf_rot.Outlier(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);

        mean_err_rf.Outlier(:,1:length(mtrys),j) = mean(err_rf.Outlier,3);
        mean_err_rerf.Outlier(:,1:length(mtrys),j) = mean(err_rerf.Outlier,3);
        mean_err_rerfr.Outlier(:,1:length(mtrys),j) = mean(err_rerfr.Outlier,3);
        mean_err_rerfdn.Outlier(:,1:length(mtrys),j) = mean(err_rerfdn.Outlier,3);
        mean_err_rf_rot.Outlier(:,1:length(mtrys),j) = mean(err_rf_rot.Outlier,3);
        mean_err_rf.Outlier(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);
        mean_err_rerf.Outlier(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);
        mean_err_rerfr.Outlier(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);
        mean_err_rerfdn.Outlier(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);
        mean_err_rf_rot.Outlier(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);
    else
        %Untransformed
        sem_rf.Untransformed(:,:,j) = std(err_rf.Untransformed,[],3)/sqrt(ntrials);
        sem_rerf.Untransformed(:,:,j) = std(err_rerf.Untransformed,[],3)/sqrt(ntrials);
        sem_rerfr.Untransformed(:,:,j) = std(err_rerfr.Untransformed,[],3)/sqrt(ntrials);
        sem_rerfdn.Untransformed(:,:,j) = std(err_rerfdn.Untransformed,[],3)/sqrt(ntrials);
        sem_rf_rot.Untransformed(:,:,j) = std(err_rf_rot.Untransformed,[],3)/sqrt(ntrials);

        var_rf.Untransformed(:,:,j) = var(err_rf.Untransformed,0,3);
        var_rerf.Untransformed(:,:,j) = var(err_rerf.Untransformed,0,3);
        var_rerfr.Untransformed(:,:,j) = var(err_rerfr.Untransformed,0,3);
        var_rerfdn.Untransformed(:,:,j) = var(err_rerfdn.Untransformed,0,3);
        var_rf_rot.Untransformed(:,:,j) = var(err_rf_rot.Untransformed,0,3);

        mean_err_rf.Untransformed(:,:,j) = mean(err_rf.Untransformed,3);
        mean_err_rerf.Untransformed(:,:,j) = mean(err_rerf.Untransformed,3);
        mean_err_rerfr.Untransformed(:,:,j) = mean(err_rerfr.Untransformed,3);
        mean_err_rerfdn.Untransformed(:,:,j) = mean(err_rerfdn.Untransformed,3);
        mean_err_rf_rot.Untransformed(:,:,j) = mean(err_rf_rot.Untransformed,3);
        
        %Rotated
        sem_rf.Rotated(:,:,j) = std(err_rf.Rotated,[],3)/sqrt(ntrials);
        sem_rerf.Rotated(:,:,j) = std(err_rerf.Rotated,[],3)/sqrt(ntrials);
        sem_rerfr.Rotated(:,:,j) = std(err_rerfr.Rotated,[],3)/sqrt(ntrials);
        sem_rerfdn.Rotated(:,:,j) = std(err_rerfdn.Rotated,[],3)/sqrt(ntrials);
        sem_rf_rot.Rotated(:,:,j) = std(err_rf_rot.Rotated,[],3)/sqrt(ntrials);

        var_rf.Rotated(:,:,j) = var(err_rf.Rotated,0,3);
        var_rerf.Rotated(:,:,j) = var(err_rerf.Rotated,0,3);
        var_rerfr.Rotated(:,:,j) = var(err_rerfr.Rotated,0,3);
        var_rerfdn.Rotated(:,:,j) = var(err_rerfdn.Rotated,0,3);
        var_rf_rot.Rotated(:,:,j) = var(err_rf_rot.Rotated,0,3);

        mean_err_rf.Rotated(:,:,j) = mean(err_rf.Rotated,3);
        mean_err_rerf.Rotated(:,:,j) = mean(err_rerf.Rotated,3);
        mean_err_rerfr.Rotated(:,:,j) = mean(err_rerfr.Rotated,3);
        mean_err_rerfdn.Rotated(:,:,j) = mean(err_rerfdn.Rotated,3);
        mean_err_rf_rot.Rotated(:,:,j) = mean(err_rf_rot.Rotated,3);
        
        %Scaled
        sem_rf.Scaled(:,:,j) = std(err_rf.Scaled,[],3)/sqrt(ntrials);
        sem_rerf.Scaled(:,:,j) = std(err_rerf.Scaled,[],3)/sqrt(ntrials);
        sem_rerfr.Scaled(:,:,j) = std(err_rerfr.Scaled,[],3)/sqrt(ntrials);
        sem_rerfdn.Scaled(:,:,j) = std(err_rerfdn.Scaled,[],3)/sqrt(ntrials);
        sem_rf_rot.Scaled(:,:,j) = std(err_rf_rot.Scaled,[],3)/sqrt(ntrials);

        var_rf.Scaled(:,:,j) = var(err_rf.Scaled,0,3);
        var_rerf.Scaled(:,:,j) = var(err_rerf.Scaled,0,3);
        var_rerfr.Scaled(:,:,j) = var(err_rerfr.Scaled,0,3);
        var_rerfdn.Scaled(:,:,j) = var(err_rerfdn.Scaled,0,3);
        var_rf_rot.Scaled(:,:,j) = var(err_rf_rot.Scaled,0,3);

        mean_err_rf.Scaled(:,:,j) = mean(err_rf.Scaled,3);
        mean_err_rerf.Scaled(:,:,j) = mean(err_rerf.Scaled,3);
        mean_err_rerfr.Scaled(:,:,j) = mean(err_rerfr.Scaled,3);
        mean_err_rerfdn.Scaled(:,:,j) = mean(err_rerfdn.Scaled,3);
        mean_err_rf_rot.Scaled(:,:,j) = mean(err_rf_rot.Scaled,3);
        
        %Affine
        sem_rf.Affine(:,:,j) = std(err_rf.Affine,[],3)/sqrt(ntrials);
        sem_rerf.Affine(:,:,j) = std(err_rerf.Affine,[],3)/sqrt(ntrials);
        sem_rerfr.Affine(:,:,j) = std(err_rerfr.Affine,[],3)/sqrt(ntrials);
        sem_rerfdn.Affine(:,:,j) = std(err_rerfdn.Affine,[],3)/sqrt(ntrials);
        sem_rf_rot.Affine(:,:,j) = std(err_rf_rot.Affine,[],3)/sqrt(ntrials);

        var_rf.Affine(:,:,j) = var(err_rf.Affine,0,3);
        var_rerf.Affine(:,:,j) = var(err_rerf.Affine,0,3);
        var_rerfr.Affine(:,:,j) = var(err_rerfr.Affine,0,3);
        var_rerfdn.Affine(:,:,j) = var(err_rerfdn.Affine,0,3);
        var_rf_rot.Affine(:,:,j) = var(err_rf_rot.Affine,0,3);

        mean_err_rf.Affine(:,:,j) = mean(err_rf.Affine,3);
        mean_err_rerf.Affine(:,:,j) = mean(err_rerf.Affine,3);
        mean_err_rerfr.Affine(:,:,j) = mean(err_rerfr.Affine,3);
        mean_err_rerfdn.Affine(:,:,j) = mean(err_rerfdn.Affine,3);
        mean_err_rf_rot.Affine(:,:,j) = mean(err_rf_rot.Affine,3);
        
        %Outlier
        sem_rf.Outlier(:,:,j) = std(err_rf.Outlier,[],3)/sqrt(ntrials);
        sem_rerf.Outlier(:,:,j) = std(err_rerf.Outlier,[],3)/sqrt(ntrials);
        sem_rerfr.Outlier(:,:,j) = std(err_rerfr.Outlier,[],3)/sqrt(ntrials);
        sem_rerfdn.Outlier(:,:,j) = std(err_rerfdn.Outlier,[],3)/sqrt(ntrials);
        sem_rf_rot.Outlier(:,:,j) = std(err_rf_rot.Outlier,[],3)/sqrt(ntrials);

        var_rf.Outlier(:,:,j) = var(err_rf.Outlier,0,3);
        var_rerf.Outlier(:,:,j) = var(err_rerf.Outlier,0,3);
        var_rerfr.Outlier(:,:,j) = var(err_rerfr.Outlier,0,3);
        var_rerfdn.Outlier(:,:,j) = var(err_rerfdn.Outlier,0,3);
        var_rf_rot.Outlier(:,:,j) = var(err_rf_rot.Outlier,0,3);

        mean_err_rf.Outlier(:,:,j) = mean(err_rf.Outlier,3);
        mean_err_rerf.Outlier(:,:,j) = mean(err_rerf.Outlier,3);
        mean_err_rerfr.Outlier(:,:,j) = mean(err_rerfr.Outlier,3);
        mean_err_rerfdn.Outlier(:,:,j) = mean(err_rerfdn.Outlier,3);
        mean_err_rf_rot.Outlier(:,:,j) = mean(err_rf_rot.Outlier,3);
    end

end

save([rerfPath 'RandomerForest/Results/Trunk_transformations.mat'],'dims','sem_rf','sem_rerf','sem_rerfr','sem_rerfdn','sem_rf_rot',...
    'var_rf','var_rerf','var_rerfr','var_rerfdn','var_rf_rot',...
    'mean_err_rf','mean_err_rerf','mean_err_rerfr','mean_err_rerfdn','mean_err_rf_rot')