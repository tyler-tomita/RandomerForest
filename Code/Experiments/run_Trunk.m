close all
clear
clc

dims = [10 100 500];

for j = 1:length(dims)
    
    d = dims(j);
    n = 100;
    ntrees = 1000;
    ntrials = 10;
    NWorkers = 2;
    mtrys = ceil(d.^[0 1/4 1/2 3/4 1]);
    err_rf = zeros(ntrees,length(mtrys),ntrials);
    err_rerf = zeros(ntrees,length(mtrys),ntrials);
    err_rerfdn = zeros(ntrees,length(mtrys),ntrials);
    err_rf_rot = zeros(ntrees,length(mtrys),ntrials);
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

        i = 1;

        for mtry = mtrys

            fprintf('mtry = %d\n',mtry)
            
            rf = rpclassificationforest(ntrees,X,Ystr,'RandomForest',true,'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            err_rf(:,i,trial) = oobpredict(rf,X,Ystr,'every');

            rerf = rpclassificationforest(ntrees,X,Ystr,'sparsemethod','sparse','nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            err_rerf(:,i,trial) = oobpredict(rerf,X,Ystr,'every');

            rerfdn = rpclassificationforest(ntrees,X,Ystr,'sparsemethod','sparse','mdiff','node','nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            err_rerfdn(:,i,trial) = oobpredict(rerfdn,X,Ystr,'every');
            
            rf_rot = rpclassificationforest(ntrees,X,Ystr,'RandomForest',true,'rotate',true,'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            err_rf_rot(:,i,trial) = oobpredict(rf_rot,X,Ystr,'every');

            i = i + 1;
        end
    end

    save(sprintf('~/RandomerForest/Results/Trunk_parameter_selection_n%d_d%d.mat',n,d),'err_rf','err_rerf','err_rerfdn','err_rf_rot')

    sem_rf(:,:,j) = std(err_rf,[],3)/sqrt(ntrials);
    sem_rerf(:,:,j) = std(err_rerf,[],3)/sqrt(ntrials);
    sem_rerfdn(:,:,j) = std(err_rerfdn,[],3)/sqrt(ntrials);
    sem_rf_rot(:,:,j) = std(err_rf_rot,[],3)/sqrt(ntrials);

    var_rf(:,:,j) = var(err_rf,0,3);
    var_rerf(:,:,j) = var(err_rerf,0,3);
    var_rerfdn(:,:,j) = var(err_rerfdn,0,3);
    var_rf_rot(:,:,j) = var(err_rf_rot,0,3);

    mean_err_rf(:,:,j) = mean(err_rf,3);
    mean_err_rerf(:,:,j) = mean(err_rerf,3);
    mean_err_rerfdn(:,:,j) = mean(err_rerfdn,3);
    mean_err_rf_rot(:,:,j) = mean(err_rf_rot,3);

end

save('~/RandomerForest/Results/Trunk.mat','sem_rf','sem_rerf','sem_rerfdn','sem_rf_rot',...
    'var_rf','var_rerf','var_rerfdn','var_rf_rot',...
    'mean_err_rf','mean_err_rerf','mean_err_rerfdn','mean_err_rf_rot')