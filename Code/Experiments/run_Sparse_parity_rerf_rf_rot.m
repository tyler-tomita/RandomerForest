close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

rng(1);

load Sparse_parity_data
ndims = length(dims);
ntrees = 500;
NWorkers = 2;

Lhat.rerf = NaN(ndims,5,ntrials);
Lhat.rf_rot = NaN(ndims,5,ntrials);
trainTime.rerf = NaN(ndims,5,ntrials);
trainTime.rf_rot = NaN(ndims,5,ntrials);

for i = 1:ndims

    d = dims(i);
    
    if d <= 5
        mtrys = 1:d;
    else
        mtrys = ceil(d.^[0 1/4 1/2 3/4 1]);
    endw
    
    for trial = 1:ntrials

        fprintf('trial %d\n',trial)

        for j = 1:length(mtrys)
            
            mtry = mtrys(j);

            fprintf('mtry = %d\n',mtry)
            
            poolobj = gcp('nocreate');
            if isempty(poolobj)
                parpool('local',NWorkers);
            end
            
            tic;
            cl.rerf = rpclassificationforest(ntrees,X{i}(:,:,trial),Y{i}(:,trial),'sparsemethod','sparse','nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            trainTime.rerf(i,j,trial) = toc;
            Lhat.rerf(i,j,trial) = oobpredict(cl.rerf,X{i}(:,:,trial),Y{i}(:,trial),'last');

            tic;
            cl.rf_rot = rpclassificationforest(ntrees,X{i}(:,:,trial),Y{i}(:,trial),'RandomForest',true,'rotate',true,'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            trainTime.rf_rot(i,j,trial) = toc;
            Lhat.rf_rot(i,j,trial) = oobpredict(cl.rf_rot,X{i}(:,:,trial),Y{i}(:,trial),'last');
        end
    end
end

save([rerfPath 'RandomerForest/Results/Sparse_parity_rerf_rf_rot.mat'],'dims','Lhat','trainTime')