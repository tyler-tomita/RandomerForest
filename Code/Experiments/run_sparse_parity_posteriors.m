%% Plot Sparse Parity Posterior Heat Maps 
% Trains RF, RerF, RerFdn, and Rotation RF on sparse parity and plots
% posterior heat maps for each classifier

%% Initialize parameters
close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

n = 1000;
d = 20;
dgood = 3;
ntrees = 1000;
ntrials = 10;
NWorkers = 2;
if d <= 5
    mtrys = 1:d;
else
    mtrys = ceil(d.^[0 1/4 1/2 3/4 1]);
end
rf.Lhat = zeros(ntrees,length(mtrys),ntrials);
rerf.Lhat = zeros(ntrees,length(mtrys),ntrials);
rerfdn.Lhat = zeros(ntrees,length(mtrys),ntrials);
rf_rot.Lhat = zeros(ntrees,length(mtrys),ntrials);
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

for trial = 1:ntrials
    
    %% Generate datapoints
    X = zeros(n,d);
    %Sigma = 1/8*ones(1,d);
    Sigma = 1/32*ones(1,d);
    %nones = randi(d+1,n,1)-1;
    %Y = mod(nones,2);
    %Ystr = cellstr(num2str(Y));
    Mu = sparse(n,d);
    for jj = 1:n
        %onesidx = randsample(1:d,nones(j),false);
        %Mu(j,onesidx) = 1;
        Mu(jj,:) = binornd(1,0.5,1,d);
        X(jj,1:d) = mvnrnd(Mu(jj,:),Sigma);
    end

    nones = sum(Mu(:,1:dgood),2);
    Ystr = cellstr(num2str(mod(nones,2)));

    %% Train classifiers
    i = 1;

    for mtry = mtrys

        fprintf('mtry = %d\n',mtry)

        rf.cl{trial,i} = rpclassificationforest(ntrees,X,Ystr,'RandomForest',true,'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
        rf.Lhat(:,i,trial) = oobpredict(rf.cl{trial,i},X,Ystr,'every');

        rerf.cl{trial,i} = rpclassificationforest(ntrees,X,Ystr,'sparsemethod','sparse','nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
        rerf.Lhat(:,i,trial) = oobpredict(rerf.cl{trial,i},X,Ystr,'every');

        rerfdn.cl{trial,i} = rpclassificationforest(ntrees,X,Ystr,'sparsemethod','sparse','mdiff','node','nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
        rerfdn.Lhat(:,i,trial) = oobpredict(rerfdn.cl{trial,i},X,Ystr,'every');

        rf_rot.cl{trial,i} = rpclassificationforest(ntrees,X,Ystr,'RandomForest',true,'rotate',true,'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
        rf_rot.Lhat(:,i,trial) = oobpredict(rf_rot.cl{trial,i},X,Ystr,'every');

        i = i + 1;
    end

    %% Select best hyperparameters for each algorithm
    [~,rf.minIdx(trial)] = min(rf.Lhat(end,:,trial));
    [~,rerf.minIdx(trial)] = min(rerf.Lhat(end,:,trial));
    [~,rerfdn.minIdx(trial)] = min(rerfdn.Lhat(end,:,trial));
    [~,rf_rot.minIdx(trial)] = min(rf_rot.Lhat(end,:,trial));

    rf.posteriors(:,:,trial) = rerf_classprob(rf.cl{trial,rf.minIdx(trial)},[Xpost Ypost Zpost]);
    rerf.posteriors(:,:,trial) = rerf_classprob(rerf.cl{trial,rerf.minIdx(trial)},[Xpost Ypost Zpost]);
    rerfdn.posteriors(:,:,trial) = rerf_classprob(rerfdn.cl{trial,rerfdn.minIdx(trial)},[Xpost Ypost Zpost]);
    rf_rot.posteriors(:,:,trial) = rerf_classprob(rf_rot.cl{trial,rf_rot.minIdx(trial)},[Xpost Ypost Zpost]);
    
    rf = rmfield(rf,'cl');
    rerf = rmfield(rerf,'cl');
    rerfdn = rmfield(rerfdn,'cl');
    rf_rot = rmfield(rf_rot,'cl');

end

save([rerfPath 'RandomerForest/Results/Sparse_parity_posteriors_npoints_50.mat'],...
    'n','d','dgood','ntrees','mtrys','rf','rerf','rerfdn','rf_rot',...
    'Xpost','Ypost','Zpost')