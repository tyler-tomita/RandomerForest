%% Plot Sparse Parity Posterior Heat Maps 
% Trains RF, RerF, RerFdn, and Rotation RF on sparse parity and plots
% posterior heat maps for each classifier

%% Initialize parameters
close all
clear
clc

n = 1000;
d = 20;
dgood = 3;
ntrees = 1000;
NWorkers = 2;
if d <= 5
    mtrys = 1:d;
else
    mtrys = ceil(d.^[0 1/4 1/2 3/4 1]);
end
rf.Lhat = zeros(ntrees,length(mtrys));
rerf.Lhat = zeros(ntrees,length(mtrys));
rerfdn.Lhat = zeros(ntrees,length(mtrys));
rf_rot.Lhat = zeros(ntrees,length(mtrys));
Class = [0;1];

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

xmin = min(X(:,1));
xmax = max(X(:,1));
ymin = min(X(:,2));
ymax = max(X(:,2));

%% Train classifiers
i = 1;

for mtry = mtrys

    fprintf('mtry = %d\n',mtry)

    rf.cl{i} = rpclassificationforest(ntrees,X,Ystr,'RandomForest',true,'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
    rf.Lhat(:,i) = oobpredict(rf.cl{i},X,Ystr,'every');

    rerf.cl{i} = rpclassificationforest(ntrees,X,Ystr,'sparsemethod','sparse','nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
    rerf.Lhat(:,i) = oobpredict(rerf.cl{i},X,Ystr,'every');

    rerfdn.cl{i} = rpclassificationforest(ntrees,X,Ystr,'sparsemethod','sparse','mdiff','node','nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
    rerfdn.Lhat(:,i) = oobpredict(rerfdn.cl{i},X,Ystr,'every');

    rf_rot.cl{i} = rpclassificationforest(ntrees,X,Ystr,'RandomForest',true,'rotate',true,'nvartosample',mtry,'NWorkers',NWorkers,'Stratified',true);
    rf_rot.Lhat(:,i) = oobpredict(rf_rot.cl{i},X,Ystr,'every');

    i = i + 1;
end

%% Select best hyperparameters for each algorithm
[~,rf.minIdx] = min(rf.Lhat(end,:));
[~,rerf.minIdx] = min(rerf.Lhat(end,:));
[~,rerfdn.minIdx] = min(rerfdn.Lhat(end,:));
[~,rf_rot.minIdx] = min(rf_rot.Lhat(end,:));

%% Plot posterior heat maps
npoints = 200;
[xgv,ygv] = meshgrid(linspace(xmin,xmax,npoints),linspace(ymin,ymax,npoints));
X = xgv(:);
Y = ygv(:);
Z = zeros(npoints^2,d-2);

rf.posteriors = rerf_classprob(rf.cl{rf.minIdx},[X Y Z]);
p1 = posterior_map(X,Y,rf.posteriors);
xlabel('x1')
ylabel('x2')
title('RF')

figure(2)
rerf.posteriors = rerf_classprob(rerf.cl{rerf.minIdx},[X Y Z]);
p2 = posterior_map(X,Y,rerf.posteriors);
xlabel('x1')
ylabel('x2')
title('RerF')

figure(3)
rerfdn.posteriors = rerf_classprob(rerfdn.cl{rerfdn.minIdx},[X Y Z]);
p3 = posterior_map(X,Y,rerfdn.posteriors);
xlabel('x1')
ylabel('x2')
title('RerFdn')

figure(4)
rf_rot.posteriors = rerf_classprob(rf_rot.cl{rf_rot.minIdx},[X Y Z]);
p4 = posterior_map(X,Y,rf_rot.posteriors);
xlabel('x1')
ylabel('x2')
title('Rotation RF')

cmin = min([p1.CData(:);p2.CData(:);p3.CData(:);p4.CData(:)]);
cmax = max([p1.CData(:);p2.CData(:);p3.CData(:);p4.CData(:)]);

for i = 1:4
    figure(i)
    caxis([cmin cmax])
end