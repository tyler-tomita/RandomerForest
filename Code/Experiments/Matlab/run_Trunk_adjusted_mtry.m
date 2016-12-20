close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

rng(1);

load Trunk_data
load Random_matrix_adjustment_factor

ntrees = 1000;
ntrials = 25;
NWorkers = 24;
Class = [0;1];
Lhat.rerf = NaN(length(dims),7,ntrials);
Lhat.rerfr = NaN(length(dims),7,ntrials);
Lhat.rerfdn = NaN(length(dims),7,ntrials);
trainTime.rerf = NaN(length(dims),7,ntrials);
trainTime.rerfr = NaN(length(dims),7,ntrials);
trainTime.rerfdn = NaN(length(dims),7,ntrials);

for j = 1:length(dims)
    
    d = dims(j);
    
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
            dprime = ceil(mtry^(1/interp1(ps,slope,d)));

            fprintf('mtry = %d\n',mtry)

            poolobj = gcp('nocreate');
            if isempty(poolobj)
                parpool('local',NWorkers);
            end

            tic;
            cl.rerf = rpclassificationforest(ntrees,X{j}(:,:,trial),...
                Y{j}(:,trial),'sparsemethod','sparse-adjusted',...
                'nvartosample',mtry,'NWorkers',NWorkers,...
                'Stratified',true,'dprime',dprime);
            trainTime.rerf(j,i,trial) = toc;
            Predictions = oobpredict(cl.rerf,X{j}(:,:,trial),Y{j}(:,trial));
            Lhat.rerf(j,i,trial) = oob_error(Predictions,Y{j}(:,trial),'last');
            
            tic;
            cl.rerfr = rpclassificationforest(ntrees,X{j}(:,:,trial),...
                Y{j}(:,trial),'sparsemethod','sparse-adjusted',...
                'Robust',true,'nvartosample',mtry,'NWorkers',NWorkers,...
                'Stratified',true,'dprime',dprime);
            trainTime.rerfr(j,i,trial) = toc;
            Predictions = oobpredict(cl.rerfr,X{j}(:,:,trial),Y{j}(:,trial));
            Lhat.rerfr(j,i,trial) = oob_error(Predictions,Y{j}(:,trial),'last');

            tic;
            cl.rerfdn = rpclassificationforest(ntrees,X{j}(:,:,trial),...
                Y{j}(:,trial),'sparsemethod','sparse-adjusted',...
                'mdiff','node','nvartosample',mtry,'NWorkers',NWorkers,...
                'Stratified',true,'dprime',dprime);
            trainTime.rerfdn(j,i,trial) = toc;
            Predictions = oobpredict(cl.rerfdn,X{j}(:,:,trial),Y{j}(:,trial));
            Lhat.rerfdn(j,i,trial) = oob_error(Predictions,Y{j}(:,trial),'last');

            clear cl
        end
    end
end

save([rerfPath 'RandomerForest/Results/Trunk_adjusted_mtry.mat'],'dims','Lhat','trainTime')
