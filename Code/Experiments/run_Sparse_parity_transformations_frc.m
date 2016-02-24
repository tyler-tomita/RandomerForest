close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

rng(1);

load Sparse_parity_transformations_data
ntrees = 500;
NWorkers = 2;

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
    
    err_frc.Untransformed = zeros(ntrees,length(mtrys),ntrials);
    
    err_frc.Rotated = zeros(ntrees,length(mtrys),ntrials);
    
    err_frc.Scaled = zeros(ntrees,length(mtrys),ntrials);
    
    err_frc.Affine = zeros(ntrees,length(mtrys),ntrials);
    
    err_frc.Outlier = zeros(ntrees,length(mtrys),ntrials);
    
    Class = [0;1];

    for trial = 1:ntrials

        fprintf('trial %d\n',trial)

        for j = 1:length(mtrys)
            mtry = mtrys(j);

            fprintf('mtry = %d\n',mtry)
            
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

save([rerfPath 'RandomerForest/Results/Sparse_parity_transformations_frc.mat'],...
    'dims','sem_frc','var_frc','mean_err_frc')