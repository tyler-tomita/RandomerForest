close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

rng(1);

n = 1000;
ntrees = 500;
ntrials = 25;
NWorkers = 2;
Class = [0;1];
dims = [2 5 10 25 50 100];
ndims = length(dims);
Lhat.rf = NaN(ndims,5,ntrials);
Lhat.rerf = NaN(ndims,5,ntrials);
Lhat.rerfdn = NaN(ndims,5,ntrials);
Lhat.rf_rot = NaN(ndims,5,ntrials);
trainTime.rf = NaN(ndims,5,ntrials);
trainTime.rerf = NaN(ndims,5,ntrials);
trainTime.rerfdn = NaN(ndims,5,ntrials);
trainTime.rf_rot = NaN(ndims,5,ntrials);

for j = 1:length(dims)
    
    d = dims(j);
    dgood = min(3,d);
    
    if d <= 5
        mtrys = 1:d;
    else
        mtrys = ceil(d.^[0 1/4 1/2 3/4 1]);
    end

    for trial = 1:ntrials

        fprintf('trial %d\n',trial)

        X = zeros(n,d);
        Sigma = 1/32*ones(1,d);
        Mu = sparse(n,d);
        for jj = 1:n
            Mu(jj,:) = binornd(1,0.5,1,d);
            X(jj,1:d) = mvnrnd(Mu(jj,:),Sigma);
        end

        nones = sum(Mu(:,1:dgood),2);
        Ystr = cellstr(num2str(mod(nones,2)));

        i = 1;

        for i = 1:length(mtrys)
            
            mtry = mtrys(i);

            fprintf('mtry = %d\n',mtry)

            poolobj = gcp('nocreate');
            if isempty(poolobj)
                parpool('local',NWorkers);
            end

            tic;
            cl.rf = rpclassificationforest(ntrees,X,Ystr,'RandomForest',true,'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            trainTime.rf(j,i,trial) = toc;
            Lhat.rf(j,i,trial) = oobpredict(cl.rf,X,Ystr,'last');

            tic;
            cl.rerf = rpclassificationforest(ntrees,X,Ystr,'sparsemethod','sparse','nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            trainTime.rerf(j,i,trial) = toc;
            Lhat.rerf(j,i,trial) = oobpredict(cl.rerf,X,Ystr,'last');

            tic;
            cl.rerfdn = rpclassificationforest(ntrees,X,Ystr,'sparsemethod','sparse','mdiff','node','nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            trainTime.rerfdn(j,i,trial) = toc;
            Lhat.rerfdn(j,i,trial) = oobpredict(cl.rerfdn,X,Ystr,'last');

            tic;
            cl.rf_rot = rpclassificationforest(ntrees,X,Ystr,'RandomForest',true,'rotate',true,'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            trainTime.rf_rot(j,i,trial) = toc;
            Lhat.rf_rot(j,i,trial) = oobpredict(cl.rf_rot,X,Ystr,'last');
        end
    end
end

save([rerfPath 'RandomerForest/Results/Sparse_parity.mat'],'dims','Lhat','trainTime')