close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

rng(1);

load Trunk_data
ndims = length(dims);
ntrees = 1000;
nmixs = 2:6;
NWorkers = 2;

Lhat.rf = NaN(ndims,5,ntrials);
Lhat.rerfr = NaN(ndims,5,ntrials);
Lhat.frc = NaN(ndims,25,ntrials);
trainTime.rf = NaN(ndims,5,ntrials);
trainTime.rerfr = NaN(ndims,5,ntrials);
trainTime.frc = NaN(ndims,25,ntrials);

for i = 1:ndims

    d = dims(i);
    
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
            cl.rf = rpclassificationforest(ntrees,X{i}(:,:,trial),Y{i}(:,trial),'RandomForest',true,'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            trainTime.rf(i,j,trial) = toc;
            Lhat.rf(i,j,trial) = oobpredict(cl.rf,X{i}(:,:,trial),Y{i}(:,trial),'last');

            tic;
            cl.rerfr = rpclassificationforest(ntrees,X{i}(:,:,trial),Y{i}(:,trial),'sparsemethod','sparse','Robust',true,'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            trainTime.rerfr(i,j,trial) = toc;
            Lhat.rerfr(i,j,trial) = oobpredict(cl.rerfr,X{i}(:,:,trial),Y{i}(:,trial),'last');
    
            for k = 1:length(nmixs)
                nmix = nmixs(k);
                fprintf('nmix = %d\n',nmix)
                tic;
                cl.frc = rpclassificationforest(ntrees,X{i}(:,:,trial),Y{i}(:,trial),'sparsemethod','frc','nmix',nmix,'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
                trainTime.frc(i,length(nmixs)*(j-1)+k,trial) = toc;
                Lhat.frc(i,length(nmixs)*(j-1)+k,trial) = oobpredict(cl.frc,X{i}(:,:,trial),Y{i}(:,trial),'last');
            end
        end
    end
end

save([rerfPath 'RandomerForest/Results/Trunk_rerfr_frc.mat'],'dims','Lhat','trainTime')