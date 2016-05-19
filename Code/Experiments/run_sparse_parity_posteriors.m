%% Plot Sparse Parity Posterior Heat Maps 
% Trains RF, RerF, RerFdn, and Rotation RF on sparse parity and plots
% posterior heat maps for each classifier

%% Initialize parameters
close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

rng(1);

load Sparse_parity_transformations_data

[n,d,ntrials] = size(X{4});
ntrees = 500;
NWorkers = 24;
mtrys = ceil(d.^[0 1/4 1/2 3/4 1 1.5 2]);
rf.Lhat = zeros(ntrees,sum(mtrys<=d),ntrials);
rerf.Lhat = zeros(ntrees,length(mtrys),ntrials);
rerfr.Lhat = zeros(ntrees,length(mtrys),ntrials);
rerfdn.Lhat = zeros(ntrees,length(mtrys),ntrials);
rf_rot.Lhat = zeros(ntrees,sum(mtrys<=d),ntrials);
Class = [0;1];

xmin = -.5;
xmax = 1.5;
ymin = xmin;
ymax = xmax;
npoints = 50;
[xgv,ygv] = meshgrid(linspace(xmin,xmax,npoints),linspace(ymin,ymax,npoints));
Xpost = xgv(:);
Ypost = ygv(:);
Zpost = zeros(npoints^2,d-2);

fprintf('Untransformed\n')
for trial = 1:ntrials
    fprintf('trial %d\n',trial)
    
    x = X{4}(:,:,trial);
    y = Y{4}(:,trial);

    %% Train classifiers
    i = 1;

    for mtry = mtrys

        fprintf('mtry = %d\n',mtry)

        if mtry <= d
            rf.cl{i} = rpclassificationforest(ntrees,x,y,'RandomForest',true,'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            Predictions = oobpredict(rf.cl{i},x,y);
            rf.Lhat(:,i,trial) = oob_error(Predictions,y,'every');
        end

        rerf.cl{i} = rpclassificationforest(ntrees,x,y,'sparsemethod','sparse','nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
        Predictions = oobpredict(rerf.cl{i},x,y);
        rerf.Lhat(:,i,trial) = oob_error(Predictions,y,'every');
        
        rerfr.cl{i} = rpclassificationforest(ntrees,x,y,'sparsemethod','sparse','Robust',true,'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
        Predictions = oobpredict(rerfr.cl{i},x,y);
        rerfr.Lhat(:,i,trial) = oob_error(Predictions,y,'every');

        rerfdn.cl{i} = rpclassificationforest(ntrees,x,y,'sparsemethod','sparse','mdiff','node','nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
        Predictions = oobpredict(rerfdn.cl{i},x,y);
        rerfdn.Lhat(:,i,trial) = oob_error(Predictions,y,'every');

        if mtry <= d
            rf_rot.cl{i} = rpclassificationforest(ntrees,x,y,'RandomForest',true,'rotate',true,'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            Predictions = oobpredict(rf_rot.cl{i},x,y);
            rf_rot.Lhat(:,i,trial) = oob_error(Predictions,y,'every');
        end

        i = i + 1;
    end

    %% Select best hyperparameters for each algorithm
    [~,rf.minIdx(trial)] = min(rf.Lhat(end,:,trial));
    [~,rerf.minIdx(trial)] = min(rerf.Lhat(end,:,trial));
    [~,rerfr.minIdx(trial)] = min(rerfr.Lhat(end,:,trial));
    [~,rerfdn.minIdx(trial)] = min(rerfdn.Lhat(end,:,trial));
    [~,rf_rot.minIdx(trial)] = min(rf_rot.Lhat(end,:,trial));

    rf.posteriors(:,:,trial) = rerf_classprob(rf.cl{rf.minIdx(trial)},[Xpost Ypost Zpost]);
    rerf.posteriors(:,:,trial) = rerf_classprob(rerf.cl{rerf.minIdx(trial)},[Xpost Ypost Zpost]);
    rerfr.posteriors(:,:,trial) = rerf_classprob(rerfr.cl{rerf.minIdx(trial)},[Xpost Ypost Zpost],'last',x);
    rerfdn.posteriors(:,:,trial) = rerf_classprob(rerfdn.cl{rerfdn.minIdx(trial)},[Xpost Ypost Zpost]);
    rf_rot.posteriors(:,:,trial) = rerf_classprob(rf_rot.cl{rf_rot.minIdx(trial)},[Xpost Ypost Zpost]);
    
    rf = rmfield(rf,'cl');
    rerf = rmfield(rerf,'cl');
    rerfr = rmfield(rerfr,'cl');
    rerfdn = rmfield(rerfdn,'cl');
    rf_rot = rmfield(rf_rot,'cl');
end

save([rerfPath 'RandomerForest/Results/Sparse_parity_posteriors_npoints_50.mat'],...
    'n','d','dgood','ntrees','mtrys','rf','rerf','rerfr','rerfdn','rf_rot',...
    'Xpost','Ypost','Zpost')

%% Scaled
clear rf rerf rerfr rerfdn rf_rot

rf.Lhat = zeros(ntrees,sum(mtrys<=d),ntrials);
rerf.Lhat = zeros(ntrees,length(mtrys),ntrials);
rerfr.Lhat = zeros(ntrees,length(mtrys),ntrials);
rerfdn.Lhat = zeros(ntrees,length(mtrys),ntrials);
rf_rot.Lhat = zeros(ntrees,sum(mtrys<=d),ntrials);

ScaleFactors = [-5,0,5];
Xspost = Xpost*10^ScaleFactors(1);
Yspost = Ypost*10^ScaleFactors(2);
Zspost = Zpost;
fprintf('Scaled\n')
for trial = 1:ntrials
    
    fprintf('trial %d\n',trial)

    x = X_scale{4}(:,:,trial);
    y = Y{4}(:,trial);
    %% Train classifiers
    i = 1;

    for mtry = mtrys

        fprintf('mtry = %d\n',mtry)

        if mtry <= d
            rf.cl{i} = rpclassificationforest(ntrees,x,y,'RandomForest',true,'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            Predictions = oobpredict(rf.cl{i},x,y);
            rf.Lhat(:,i,trial) = oob_error(Predictions,y,'every');
        end

        rerf.cl{i} = rpclassificationforest(ntrees,x,y,'sparsemethod','sparse','nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
        Predictions = oobpredict(rerf.cl{i},x,y);
        rerf.Lhat(:,i,trial) = oob_error(Predictions,y,'every');
        
        rerfr.cl{i} = rpclassificationforest(ntrees,x,y,'sparsemethod','sparse','Robust',true,'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
        Predictions = oobpredict(rerfr.cl{i},x,y);
        rerfr.Lhat(:,i,trial) = oob_error(Predictions,y,'every');

        rerfdn.cl{i} = rpclassificationforest(ntrees,x,y,'sparsemethod','sparse','mdiff','node','nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
        Predictions = oobpredict(rerfdn.cl{i},x,y);
        rerfdn.Lhat(:,i,trial) = oob_error(Predictions,y,'every');

        if mtry <= d
            rf_rot.cl{i} = rpclassificationforest(ntrees,x,y,'RandomForest',true,'rotate',true,'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            Predictions = oobpredict(rf_rot.cl{i},x,y);
            rf_rot.Lhat(:,i,trial) = oob_error(Predictions,y,'every');
        end

        i = i + 1;
    end

    %% Select best hyperparameters for each algorithm
    [~,rf.minIdx(trial)] = min(rf.Lhat(end,:,trial));
    [~,rerf.minIdx(trial)] = min(rerf.Lhat(end,:,trial));
    [~,rerfr.minIdx(trial)] = min(rerfr.Lhat(end,:,trial));
    [~,rerfdn.minIdx(trial)] = min(rerfdn.Lhat(end,:,trial));
    [~,rf_rot.minIdx(trial)] = min(rf_rot.Lhat(end,:,trial));

    rf.posteriors(:,:,trial) = rerf_classprob(rf.cl{rf.minIdx(trial)},[Xspost Yspost Zspost]);
    rerf.posteriors(:,:,trial) = rerf_classprob(rerf.cl{rerf.minIdx(trial)},[Xspost Yspost Zspost]);
    rerfr.posteriors(:,:,trial) = rerf_classprob(rerfr.cl{rerf.minIdx(trial)},[Xspost Yspost Zspost],'last',x);
    rerfdn.posteriors(:,:,trial) = rerf_classprob(rerfdn.cl{rerfdn.minIdx(trial)},[Xspost Yspost Zspost]);
    rf_rot.posteriors(:,:,trial) = rerf_classprob(rf_rot.cl{rf_rot.minIdx(trial)},[Xspost Yspost Zspost]);
    
    rf = rmfield(rf,'cl');
    rerf = rmfield(rerf,'cl');
    rerfr = rmfield(rerfr,'cl');
    rerfdn = rmfield(rerfdn,'cl');
    rf_rot = rmfield(rf_rot,'cl');
end

save([rerfPath 'RandomerForest/Results/Sparse_parity_posteriors_scaled_npoints_50.mat'],...
    'n','d','dgood','ntrees','mtrys','rf','rerf','rerfr','rerfdn','rf_rot',...
    'Xspost','Yspost','Zspost')

%% Rotated
clear rf rerf rerfr rerfdn rf_rot

rf.Lhat = zeros(ntrees,sum(mtrys<=d),ntrials);
rerf.Lhat = zeros(ntrees,length(mtrys),ntrials);
rerfr.Lhat = zeros(ntrees,length(mtrys),ntrials);
rerfdn.Lhat = zeros(ntrees,length(mtrys),ntrials);
rf_rot.Lhat = zeros(ntrees,sum(mtrys<=d),ntrials);

theta = pi/4;
R = [cos(theta) sin(theta);-sin(theta) cos(theta)];
Xr = [Xpost Ypost]*R;
Xrpost = Xr(:,1);
Yrpost = Xr(:,2);
Zrpost = Zpost;
fprintf('Scaled\n')
for trial = 1:ntrials
    
    fprintf('trial %d\n',trial)

    x = [X{4}(:,1:2,trial)*R X{4}(:,3:end,trial)];
    y = Y{4}(:,trial);
    %% Train classifiers
    i = 1;

    for mtry = mtrys

        fprintf('mtry = %d\n',mtry)

        if mtry <= d
            rf.cl{i} = rpclassificationforest(ntrees,x,y,'RandomForest',true,'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            Predictions = oobpredict(rf.cl{i},x,y);
            rf.Lhat(:,i,trial) = oob_error(Predictions,y,'every');
        end

        rerf.cl{i} = rpclassificationforest(ntrees,x,y,'sparsemethod','sparse','nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
        Predictions = oobpredict(rerf.cl{i},x,y);
        rerf.Lhat(:,i,trial) = oob_error(Predictions,y,'every');
        
        rerfr.cl{i} = rpclassificationforest(ntrees,x,y,'sparsemethod','sparse','Robust',true,'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
        Predictions = oobpredict(rerfr.cl{i},x,y);
        rerfr.Lhat(:,i,trial) = oob_error(Predictions,y,'every');

        rerfdn.cl{i} = rpclassificationforest(ntrees,x,y,'sparsemethod','sparse','mdiff','node','nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
        Predictions = oobpredict(rerfdn.cl{i},x,y);
        rerfdn.Lhat(:,i,trial) = oob_error(Predictions,y,'every');

        if mtry <= d
            rf_rot.cl{i} = rpclassificationforest(ntrees,x,y,'RandomForest',true,'rotate',true,'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
            Predictions = oobpredict(rf_rot.cl{i},x,y);
            rf_rot.Lhat(:,i,trial) = oob_error(Predictions,y,'every');
        end

        i = i + 1;
    end

    %% Select best hyperparameters for each algorithm
    [~,rf.minIdx(trial)] = min(rf.Lhat(end,:,trial));
    [~,rerf.minIdx(trial)] = min(rerf.Lhat(end,:,trial));
    [~,rerfr.minIdx(trial)] = min(rerfr.Lhat(end,:,trial));
    [~,rerfdn.minIdx(trial)] = min(rerfdn.Lhat(end,:,trial));
    [~,rf_rot.minIdx(trial)] = min(rf_rot.Lhat(end,:,trial));

    rf.posteriors(:,:,trial) = rerf_classprob(rf.cl{rf.minIdx(trial)},[Xrpost Yrpost Zrpost]);
    rerf.posteriors(:,:,trial) = rerf_classprob(rerf.cl{rerf.minIdx(trial)},[Xrpost Yrpost Zrpost]);
    rerfr.posteriors(:,:,trial) = rerf_classprob(rerfr.cl{rerf.minIdx(trial)},[Xrpost Yrpost Zrpost],'last',x);
    rerfdn.posteriors(:,:,trial) = rerf_classprob(rerfdn.cl{rerfdn.minIdx(trial)},[Xrpost Yrpost Zrpost]);
    rf_rot.posteriors(:,:,trial) = rerf_classprob(rf_rot.cl{rf_rot.minIdx(trial)},[Xrpost Yrpost Zrpost]);
    
    rf = rmfield(rf,'cl');
    rerf = rmfield(rerf,'cl');
    rerfr = rmfield(rerfr,'cl');
    rerfdn = rmfield(rerfdn,'cl');
    rf_rot = rmfield(rf_rot,'cl');
end

save([rerfPath 'RandomerForest/Results/Sparse_parity_posteriors_rotated_npoints_50.mat'],...
    'n','d','dgood','ntrees','mtrys','rf','rerf','rerfr','rerfdn','rf_rot',...
    'Xrpost','Yrpost','Zrpost')