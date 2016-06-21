close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

rng(1);

X = loadMNISTImages('train-images-idx3-ubyte');
X = sparse(X');
Y = loadMNISTLabels('train-labels-idx1-ubyte');
X = X((Y==3|Y==7|Y==8),:);
Y = cellstr(num2str(Y((Y==3|Y==7|Y==8))));
[n,p] = size(X);
ih = sqrt(p);
iw = ih;

load Random_matrix_adjustment_factor

ntrees = 500;
Stratified = true;
NWorkers = 24;


%% Structured RerF %%
Fs = [2 4 8 12 20];
Lhat.srerf = NaN(length(Fs),7);

fprintf('Structured RerF\n')

for i = 1:length(Fs)
    F = Fs(i);
    
    fprintf('F = %d\n',F)
    
    if F^2 <= 4
        ds = [1:F^2 ceil((F^2).^[1.5 2])];
    elseif F^2 > 4 && F^2 <= 100
        ds = ceil((F^2).^[0 1/4 1/2 3/4 1 1.5 2]);
    else
        ds = [ceil((F^2).^[0 1/4 1/2 3/4 1]) 5*F^2 10*F^2];
    end
    
    for j = 1:length(ds)
        d = ds(j);
        fprintf('d = %d\n',d)
        
        dprime = ceil(d^(1/interp1(ps,slope,p)));

        srerf = rpclassificationforest(ntrees,X,Y,'Image',true,'ih',ih,'iw',iw,...
            'F',F,'sparsemethod','sparse-adjusted','dprime',dprime,...
            'nvartosample',d,'NWorkers',NWorkers,'Stratified',Stratified);
        Predictions = oobpredict(srerf,X,Y);
        Lhat.srerf(i,j) = oob_error(Predictions,Y,'last');
    end
end

%% RerF %%
fprintf('\nRerF\n')

Lhat.rerf = NaN(1,7);

ds = [ceil(p.^[0 1/4 1/2 3/4 1]) 10*p 15*p];

for j = 1:length(ds)
    d = ds(j);
    
    fprintf('d = %d\n',d)
    
    dprime = ceil(d^(1/interp1(ps,slope,p)));
    
    rerf = rpclassificationforest(ntrees,X,Y,'sparsemethod',...
        'sparse-adjusted','nvartosample',d,'dprime',dprime,'NWorkers',NWorkers,...
        'Stratified',Stratified);
    Predictions = oobpredict(rerf,X,Y);
    Lhat.rerf(j) = oob_error(Predictions,Y,'last');
end

%% RF %%
fprintf('\nRandom Forest\n')

Lhat.rf = NaN(1,5);

ds = ceil(p.^[0 1/4 1/2 3/4 1]);

for j = 1:length(ds)
    d = ds(j);
    
    fprintf('d = %d\n',d)
    
    rf = rpclassificationforest(ntrees,X,Y,'RandomForest',true,...
        'nvartosample',d,'NWorkers',NWorkers,'Stratified',Stratified);
    Predictions = oobpredict(rf,X,Y);
    Lhat.rf(j) = oob_error(Predictions,Y,'last');
end

save([rerfPath 'RandomerForest/Results/mnist_378.mat'],'ntrees','Fs','Lhat')
