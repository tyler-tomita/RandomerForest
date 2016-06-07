close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

rng(1);

load Sparse_parity_data

ntrees = 500;
ntrials = 25;
NWorkers = 16;
Class = [0;1];

Lhat.rerfu = NaN(length(dims),7,ntrials);
trainTime.rerfu = NaN(length(dims),7,ntrials);

for j = 1:length(dims)
    
    d = dims(j);
    
    if d <= 5
        mtrys = [1:d ceil(d.^[1.5 2])];
    elseif d > 5 && d <= 25
        mtrys = ceil(d.^[0 1/4 1/2 3/4 1 1.5 2]);
    else
        mtrys = [ceil(d.^[0 1/4 1/2 3/4 1]) 5*d 10*d];
    end

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

            tic;
            cl.rerfu = rpclassificationforest(ntrees,X{j}(:,:,trial),Y{j}(:,trial),'sparsemethod','sparse','nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            trainTime.rerfu(j,i,trial) = toc;
            Predictions = oobpredict(cl.rerfu,X{j}(:,:,trial),Y{j}(:,trial));
            Lhat.rerfu(j,i,trial) = oob_error(Predictions,Y{j}(:,trial),'last');

            clear cl
        end
    end
end

save([rerfPath 'RandomerForest/Results/Trunk_rerfu.mat'],'dims','Lhat','trainTime')
