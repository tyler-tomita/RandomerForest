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
Class = [0;1];
dims = [2 10 50 100 500];

for j = 1:length(dims)
    
    d = dims(j);
    
    if d <= 5
        mtrys = 1:d;
    else
        mtrys = ceil(d.^[0 1/4 1/2 3/4 1]);
    end
    Lhat.rf = zeros(ntrees,length(mtrys),ntrials);
    Lhat.rerf = zeros(ntrees,length(mtrys),ntrials);
    Lhat.rerfdn = zeros(ntrees,length(mtrys),ntrials);
    Lhat.rf_rot = zeros(ntrees,length(mtrys),ntrials);
    trainTime.rf = zeros(ntrials,length(mtrys));
    trainTime.rerf = zeros(ntrials,length(mtrys));
    trainTime.rerfdn = zeros(ntrials,length(mtrys));
    trainTime.rf_rot = zeros(ntrials,length(mtrys));

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
            
            tic;
            cl.rf = rpclassificationforest(ntrees,X,Ystr,'RandomForest',true,'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            trainTime.rf(trial,i) = toc;
            Lhat.rf(:,i,trial) = oobpredict(cl.rf,X,Ystr,'every');

            tic;
            cl.rerf = rpclassificationforest(ntrees,X,Ystr,'sparsemethod','sparse','nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            trainTime.rerf(trial,i) = toc;
            Lhat.rerf(:,i,trial) = oobpredict(cl.rerf,X,Ystr,'every');

            tic;
            cl.rerfdn = rpclassificationforest(ntrees,X,Ystr,'sparsemethod','sparse','mdiff','node','nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            trainTime.rerfdn(trial,i) = toc;
            Lhat.rerfdn(:,i,trial) = oobpredict(cl.rerfdn,X,Ystr,'every');
            
            tic;
            cl.rf_rot = rpclassificationforest(ntrees,X,Ystr,'RandomForest',true,'rotate',true,'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            trainTime.rf_rot(trial,i) = toc;
            Lhat.rf_rot(:,i,trial) = oobpredict(cl.rf_rot,X,Ystr,'every');

            i = i + 1;
        end
    end

    save(sprintf([rerfPath 'RandomerForest/Results/Trunk_parameter_selection_n%d_d%d.mat'],n,d),'Lhat','trainTime')

    if length(mtrys) < 5
        emptyCol = 5 - length(mtrys);
        
        semLhat.rf(:,1:length(mtrys),j) = std(Lhat.rf,[],3)/sqrt(ntrials);
        semLhat.rerf(:,1:length(mtrys),j) = std(Lhat.rerf,[],3)/sqrt(ntrials);
        semLhat.rerfdn(:,1:length(mtrys),j) = std(Lhat.rerfdn,[],3)/sqrt(ntrials);
        semLhat.rf_rot(:,1:length(mtrys),j) = std(Lhat.rf_rot,[],3)/sqrt(ntrials);
        semLhat.rf(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);
        semLhat.rerf(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);
        semLhat.rerfdn(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);
        semLhat.rf_rot(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);

        meanLhat.rf(:,1:length(mtrys),j) = mean(Lhat.rf,3);
        meanLhat.rerf(:,1:length(mtrys),j) = mean(Lhat.rerf,3);
        meanLhat.rerfdn(:,1:length(mtrys),j) = mean(Lhat.rerfdn,3);
        meanLhat.rf_rot(:,1:length(mtrys),j) = mean(Lhat.rf_rot,3);
        meanLhat.rf(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);
        meanLhat.rerf(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);
        meanLhat.rerfdn(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);
        meanLhat.rf_rot(:,length(mtrys)+1:5,j) = NaN(ntrees,emptyCol);
        
        semTrainTime.rf(1,1:length(mtrys),j) = std(trainTime.rf)/sqrt(ntrials);
        semTrainTime.rerf(1,1:length(mtrys),j) = std(trainTime.rerf)/sqrt(ntrials);
        semTrainTime.rerfdn(1,1:length(mtrys),j) = std(trainTime.rerfdn)/sqrt(ntrials);
        semTrainTime.rf_rot(1,1:length(mtrys),j) = std(trainTime.rf_rot)/sqrt(ntrials);
        semTrainTime.rf(1,length(mtrys)+1:5,j) = NaN(1,emptyCol);
        semTrainTime.rerf(1,length(mtrys)+1:5,j) = NaN(1,emptyCol);
        semTrainTime.rerfdn(1,length(mtrys)+1:5,j) = NaN(1,emptyCol);
        semTrainTime.rf_rot(1,length(mtrys)+1:5,j) = NaN(1,emptyCol);
        
        meanTrainTime.rf(1,1:length(mtrys),j) = mean(trainTime.rf);
        meanTrainTime.rerf(1,1:length(mtrys),j) = mean(trainTime.rerf);
        meanTrainTime.rerfdn(1,1:length(mtrys),j) = mean(trainTime.rerfdn);
        meanTrainTime.rf_rot(1,1:length(mtrys),j) = mean(trainTime.rf_rot);
        meanTrainTime.rf(1,length(mtrys)+1:5,j) = NaN(1,emptyCol);
        meanTrainTime.rerf(1,length(mtrys)+1:5,j) = NaN(1,emptyCol);
        meanTrainTime.rerfdn(1,length(mtrys)+1:5,j) = NaN(1,emptyCol);
        meanTrainTime.rf_rot(1,length(mtrys)+1:5,j) = NaN(1,emptyCol);
    else
        semLhat.rf(:,:,j) = std(Lhat.rf,[],3)/sqrt(ntrials);
        semLhat.rerf(:,:,j) = std(Lhat.rerf,[],3)/sqrt(ntrials);
        semLhat.rerfdn(:,:,j) = std(Lhat.rerfdn,[],3)/sqrt(ntrials);
        semLhat.rf_rot(:,:,j) = std(Lhat.rf_rot,[],3)/sqrt(ntrials);

        meanLhat.rf(:,:,j) = mean(Lhat.rf,3);
        meanLhat.rerf(:,:,j) = mean(Lhat.rerf,3);
        meanLhat.rerfdn(:,:,j) = mean(Lhat.rerfdn,3);
        meanLhat.rf_rot(:,:,j) = mean(Lhat.rf_rot,3);
        
        semTrainTime.rf(1,:,j) = std(trainTime.rf)/sqrt(ntrials);
        semTrainTime.rerf(1,:,j) = std(trainTime.rerf)/sqrt(ntrials);
        semTrainTime.rerfdn(1,:,j) = std(trainTime.rerfdn)/sqrt(ntrials);
        semTrainTime.rf_rot(1,:,j) = std(trainTime.rf_rot)/sqrt(ntrials);

        meanTrainTime.rf(1,:,j) = mean(trainTime.rf);
        meanTrainTime.rerf(1,:,j) = mean(trainTime.rerf);
        meanTrainTime.rerfdn(1,:,j) = mean(trainTime.rerfdn);
        meanTrainTime.rf_rot(1,:,j) = mean(trainTime.rf_rot);
    end

end

save([rerfPath 'RandomerForest/Results/Trunk.mat'],'dims','meanLhat','semLhat','meanTrainTime','semTrainTime')