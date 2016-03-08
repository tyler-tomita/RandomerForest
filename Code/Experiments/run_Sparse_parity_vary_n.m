close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

rng(1);

load Sparse_parity_data_vary_n
nn = length(ns);
ntrees = 500;
NWorkers = 2;

Lhat.rf = NaN(nn,5,ntrials);
Lhat.rerf = NaN(nn,5,ntrials);
trainTime.rf = NaN(nn,5,ntrials);
trainTime.rerf = NaN(nn,5,ntrials);

for i = 1:nn

    n = ns(i);
    fprintf('n = %d\n',n)
    
    if d <= 5
        mtrys = 1:d;
    else
        mtrys = ceil(d.^[0 1/4 1/2 3/4 1]);
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
            cl.rerf = rpclassificationforest(ntrees,X{i}(:,:,trial),Y{i}(:,trial),'sparsemethod','sparse','nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            trainTime.rerf(i,j,trial) = toc;
            Lhat.rerf(i,j,trial) = oobpredict(cl.rerf,X{i}(:,:,trial),Y{i}(:,trial),'last');
        end
    end
end

save([rerfPath 'RandomerForest/Results/Sparse_parity_vary_n.mat'],'ns','d','Lhat','trainTime','ntrees')