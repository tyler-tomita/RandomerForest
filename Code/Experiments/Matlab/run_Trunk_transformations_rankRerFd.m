close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

rng(1);

n = 100;
ntrees = 1000;
ntrials = 10;
NWorkers = 24;
dims = [2 10 50 100 500 1000];

for j = 1:length(dims)
    
    d = dims(j);
    fprintf('dimension %d\n',d)
    if d <= 5
        mtrys = 1:d;
    else
        mtrys = ceil(d.^[0 1/4 1/2 3/4 1]);
    end
    
    err_rerfdnr.Untransformed = zeros(ntrees,length(mtrys),ntrials);
    err_rerfdnr.Rotated = zeros(ntrees,length(mtrys),ntrials);
    err_rerfdnr.Scaled = zeros(ntrees,length(mtrys),ntrials);
    err_rerfdnr.Affine = zeros(ntrees,length(mtrys),ntrials);
    err_rerfdnr.Outlier = zeros(ntrees,length(mtrys),ntrials);
    
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
            rerfdnr = rpclassificationforest(ntrees,X,Ystr,'sparsemethod','sparse','mdiff','node','Robust',true,'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            err_rerfdnr.Untransformed(:,i,trial) = oobpredict(rerfdnr,X,Ystr,'every');
            
            %Rotated

            rerfdnr = rpclassificationforest(ntrees,X_rot,Ystr,'sparsemethod','sparse','mdiff','node','Robust',true,'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            err_rerfdnr.Rotated(:,i,trial) = oobpredict(rerfdnr,X_rot,Ystr,'every');
            
            %Scaled

            rerfdnr = rpclassificationforest(ntrees,X_scale,Ystr,'sparsemethod','sparse','mdiff','node','Robust',true,'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            err_rerfdnr.Scaled(:,i,trial) = oobpredict(rerfdnr,X_scale,Ystr,'every');
            
            %Affine

            rerfdnr = rpclassificationforest(ntrees,X_affine,Ystr,'sparsemethod','sparse','mdiff','node','Robust',true,'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            err_rerfdnr.Affine(:,i,trial) = oobpredict(rerfdnr,X_affine,Ystr,'every');
            
            %Outlier

            rerfdnr = rpclassificationforest(ntrees,X_out,Y_out,'sparsemethod','sparse','mdiff','node','Robust',true,'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            err_rerfdnr.Outlier(:,i,trial) = oobpredict(rerfdnr,X_out,Y_out,'every');

            i = i + 1;
        end
    end

    save(sprintf([rerfPath 'RandomerForest/Results/Trunk_transformations_rankRerFd_n%d_d%d.mat'],n,d),'err_rerfdnr')

    if length(mtrys) < 5
        emptyCol = 5 - length(mtrys);
        
        %Untransformed
        sem_rerfdnr.Untransformed(:,1:length(mtrys),j) = std(err_rerfdnr.Untransformed,[],3)/sqrt(ntrials);
        sem_rerfdnr.Untransformed(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);

        var_rerfdnr.Untransformed(:,1:length(mtrys),j) = var(err_rerfdnr.Untransformed,0,3);
        var_rerfdnr.Untransformed(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);

        mean_err_rerfdnr.Untransformed(:,1:length(mtrys),j) = mean(err_rerfdnr.Untransformed,3);
        mean_err_rerfdnr.Untransformed(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);
        
        %Rotated
        sem_rerfdnr.Rotated(:,1:length(mtrys),j) = std(err_rerfdnr.Rotated,[],3)/sqrt(ntrials);
        sem_rerfdnr.Rotated(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);

        var_rerfdnr.Rotated(:,1:length(mtrys),j) = var(err_rerfdnr.Rotated,0,3);
        var_rerfdnr.Rotated(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);
        
        mean_err_rerfdnr.Rotated(:,1:length(mtrys),j) = mean(err_rerfdnr.Rotated,3);
        mean_err_rerfdnr.Rotated(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);
        
        %Scaled
        sem_rerfdnr.Scaled(:,1:length(mtrys),j) = std(err_rerfdnr.Scaled,[],3)/sqrt(ntrials);
        sem_rerfdnr.Scaled(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);

        var_rerfdnr.Scaled(:,1:length(mtrys),j) = var(err_rerfdnr.Scaled,0,3);
        var_rerfdnr.Scaled(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);

        mean_err_rerfdnr.Scaled(:,1:length(mtrys),j) = mean(err_rerfdnr.Scaled,3);
        mean_err_rerfdnr.Scaled(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);
        
        %Affine
        sem_rerfdnr.Affine(:,1:length(mtrys),j) = std(err_rerfdnr.Affine,[],3)/sqrt(ntrials);
        sem_rerfdnr.Affine(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);

        var_rerfdnr.Affine(:,1:length(mtrys),j) = var(err_rerfdnr.Affine,0,3);
        var_rerfdnr.Affine(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);

        mean_err_rerfdnr.Affine(:,1:length(mtrys),j) = mean(err_rerfdnr.Affine,3);
        mean_err_rerfdnr.Affine(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);
        
        %Outlier
        sem_rerfdnr.Outlier(:,1:length(mtrys),j) = std(err_rerfdnr.Outlier,[],3)/sqrt(ntrials);
        sem_rerfdnr.Outlier(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);

        var_rerfdnr.Outlier(:,1:length(mtrys),j) = var(err_rerfdnr.Outlier,0,3);
        var_rerfdnr.Outlier(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);

        mean_err_rerfdnr.Outlier(:,1:length(mtrys),j) = mean(err_rerfdnr.Outlier,3);
        mean_err_rerfdnr.Outlier(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);
    else
        %Untransformed
        sem_rerfdnr.Untransformed(:,:,j) = std(err_rerfdnr.Untransformed,[],3)/sqrt(ntrials);
        var_rerfdnr.Untransformed(:,:,j) = var(err_rerfdnr.Untransformed,0,3);
        mean_err_rerfdnr.Untransformed(:,:,j) = mean(err_rerfdnr.Untransformed,3);
        
        %Rotated
        sem_rerfdnr.Rotated(:,:,j) = std(err_rerfdnr.Rotated,[],3)/sqrt(ntrials);
        var_rerfdnr.Rotated(:,:,j) = var(err_rerfdnr.Rotated,0,3);
        mean_err_rerfdnr.Rotated(:,:,j) = mean(err_rerfdnr.Rotated,3);
        
        %Scaled
        sem_rerfdnr.Scaled(:,:,j) = std(err_rerfdnr.Scaled,[],3)/sqrt(ntrials);
        var_rerfdnr.Scaled(:,:,j) = var(err_rerfdnr.Scaled,0,3);
        mean_err_rerfdnr.Scaled(:,:,j) = mean(err_rerfdnr.Scaled,3);
        
        %Affine
        sem_rerfdnr.Affine(:,:,j) = std(err_rerfdnr.Affine,[],3)/sqrt(ntrials);
        var_rerfdnr.Affine(:,:,j) = var(err_rerfdnr.Affine,0,3);
        mean_err_rerfdnr.Affine(:,:,j) = mean(err_rerfdnr.Affine,3);
        
        %Outlier
        sem_rerfdnr.Outlier(:,:,j) = std(err_rerfdnr.Outlier,[],3)/sqrt(ntrials);
        var_rerfdnr.Outlier(:,:,j) = var(err_rerfdnr.Outlier,0,3);
        mean_err_rerfdnr.Outlier(:,:,j) = mean(err_rerfdnr.Outlier,3);
    end

end

save([rerfPath 'RandomerForest/Results/Trunk_transformations_rankRerFd.mat'],'dims','sem_rerfdnr','var_rerfdnr','mean_err_rerfdnr')