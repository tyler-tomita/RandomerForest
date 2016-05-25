close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

rng(1);

load Trunk_data

ntrees = 1000;
ntrials = 25;
NWorkers = 16;
Class = [0;1];
Lhat.rf = NaN(length(dims),5,ntrials);
Lhat.rerf = NaN(length(dims),7,ntrials);
Lhat.rerfr = NaN(length(dims),7,ntrials);
Lhat.rerfdn = NaN(length(dims),7,ntrials);
Lhat.rf_rot = NaN(length(dims),5,ntrials);
Lhat.frc = NaN(length(dims),7,ntrials);
trainTime.rf = NaN(length(dims),5,ntrials);
trainTime.rerf = NaN(length(dims),7,ntrials);
trainTime.rerfr = NaN(length(dims),7,ntrials);
trainTime.rerfdn = NaN(length(dims),7,ntrials);
trainTime.rf_rot = NaN(length(dims),5,ntrials);
trainTime.frc = NaN(length(dims),7,ntrials);

for j = 1:length(dims)
    
    d = dims(j);
    dgood = min(3,d);
    
    if d <= 5
        mtrys = [1:d ceil(d.^[1.5 2])];
    elseif d > 5 && d <= 100
        mtrys = ceil(d.^[0 1/4 1/2 3/4 1 1.5 2]);
    else
        mtrys = [ceil(d.^[0 1/4 1/2 3/4 1]) 5*d 10*d];
    end
    
    nmix = 2;

    for trial = 1:ntrials

        fprintf('trial %d\n',trial)

        i = 1;

        for i = 1:length(mtrys)
            
            mtry = mtrys(i);

            fprintf('mtry = %d\n',mtry)

            poolobj = gcp('nocreate');
            if isempty(poolobj)
                parpool('local',NWorkers);
            end

            if mtry <= d
                tic;
                cl.rf = rpclassificationforest(ntrees,X{j}(:,:,trial),Y{j}(:,trial),'RandomForest',true,'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
                trainTime.rf(j,i,trial) = toc;
                Predictions = oobpredict(cl.rf,X{j}(:,:,trial),Y{j}(:,trial));
                Lhat.rf(j,i,trial) = oob_error(Predictions,Y{j}(:,trial),'last');
            end

            tic;
            cl.rerf = rpclassificationforest(ntrees,X{j}(:,:,trial),Y{j}(:,trial),'sparsemethod','sparse','nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            trainTime.rerf(j,i,trial) = toc;
            Predictions = oobpredict(cl.rerf,X{j}(:,:,trial),Y{j}(:,trial));
            Lhat.rerf(j,i,trial) = oob_error(Predictions,Y{j}(:,trial),'last');
            
            tic;
            cl.rerfr = rpclassificationforest(ntrees,X{j}(:,:,trial),Y{j}(:,trial),'sparsemethod','sparse','Robust',true,'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            trainTime.rerfr(j,i,trial) = toc;
            Predictions = oobpredict(cl.rerfr,X{j}(:,:,trial),Y{j}(:,trial));
            Lhat.rerfr(j,i,trial) = oob_error(Predictions,Y{j}(:,trial),'last');

            tic;
            cl.rerfdn = rpclassificationforest(ntrees,X{j}(:,:,trial),Y{j}(:,trial),'sparsemethod','sparse','mdiff','node','nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            trainTime.rerfdn(j,i,trial) = toc;
            Predictions = oobpredict(cl.rerfdn,X{j}(:,:,trial),Y{j}(:,trial));
            Lhat.rerfdn(j,i,trial) = oob_error(Predictions,Y{j}(:,trial),'last');

            if mtry <= d && d <= 500
                tic;
                cl.rf_rot = rpclassificationforest(ntrees,X{j}(:,:,trial),Y{j}(:,trial),'RandomForest',true,'rotate',true,'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
                trainTime.rf_rot(j,i,trial) = toc;
                Predictions = oobpredict(cl.rf_rot,X{j}(:,:,trial),Y{j}(:,trial));
                Lhat.rf_rot(j,i,trial) = oob_error(Predictions,Y{j}(:,trial),'last');
            end
            
            tic;
            cl.frc = rpclassificationforest(ntrees,X{j}(:,:,trial),Y{j}(:,trial),'sparsemethod','frc','nmix',nmix,'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            trainTime.frc(j,i,trial) = toc;
            Predictions = oobpredict(cl.frc,X{j}(:,:,trial),Y{j}(:,trial));
            Lhat.frc(j,i,trial) = oob_error(Predictions,Y{j}(:,trial),'last');
        end
    end
end

save([rerfPath 'RandomerForest/Results/Trunk.mat'],'dims','Lhat','trainTime')
