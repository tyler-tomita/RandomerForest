close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

rng(1);

ns = [100 300 1000 3000 10000];
X = loadMNISTImages('train-images-idx3-ubyte');
X = sparse(X');
Y = loadMNISTLabels('train-labels-idx1-ubyte');
Ystr = cellstr(num2str(Y));
Labels = [3,7,8];
TrainIdx = [];
for l = 1:length(Labels)
    TrainIdx = [TrainIdx randsample(find(Y==Labels(l)),round(nTrain/length(Labels)))'];
end
Y = cellstr(num2str(Y));

[n,p] = size(X);
ih = sqrt(p);
iw = ih;

load Random_matrix_adjustment_factor

ntrees = 500;
Stratified = true;
NWorkers = 24;


for k = 1:length(ns)
    nTrain = ns(k);
    fprintf('n = %d\n',nTrain)

    TrainIdx = [];
    for l = 1:length(Labels)
        TrainIdx = [TrainIdx randsample(find(Y==Labels(l)),round(nTrain/length(Labels)))'];
    end

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
            ds = [ceil((F^2).^[0 1/4 1/2 3/4 1]) 5*F^2 15*F^2];
        end

        for j = 1:length(ds)
            d = ds(j);
            fprintf('d = %d\n',d)

            dprime = ceil(d^(1/interp1(ps,slope,p)));

            srerf = rpclassificationforest(ntrees,X(TrainIdx,:),Ystr(TrainIdx),...
                'Image',true,'ih',ih,'iw',iw,'F',F,'sparsemethod',...
                'sparse-adjusted','nvartosample',d,'dprime',dprime,...
                'NWorkers',NWorkers,'Stratified',Stratified);
            Predictions = oobpredict(srerf,X(TrainIdx,:),Ystr(TrainIdx));
            Lhat.srerf{k}(i,j) = oob_error(Predictions,Ystr(TrainIdx),'last');
        end
    end

    %% RerF %%
    fprintf('\nRerF\n')

    Lhat.rerf{k} = NaN(1,7);

    ds = [ceil(p.^[0 1/4 1/2 3/4 1]) 10*p 15*p];

    for j = 1:length(ds)
        d = ds(j);

        fprintf('d = %d\n',d)

        dprime = ceil(d^(1/interp1(ps,slope,p)));

        rerf = rpclassificationforest(ntrees,X(TrainIdx,:),Ystr(TrainIdx),'sparsemethod',...
            'sparse-adjusted','nvartosample',d,'dprime',dprime,'NWorkers',NWorkers,...
            'Stratified',Stratified);
        Predictions = oobpredict(rerf,X(TrainIdx,:),Ystr(TrainIdx));
        Lhat.rerf{k}(j) = oob_error(Predictions,Ystr(TrainIdx),'last');
    end

    %% RF %%
    fprintf('\nRandom Forest\n')

    Lhat.rf{k} = NaN(1,5);

    ds = ceil(p.^[0 1/4 1/2 3/4 1]);

    for j = 1:length(ds)
        d = ds(j);

        fprintf('d = %d\n',d)

        rf = rpclassificationforest(ntrees,X(TrainIdx,:),Ystr(TrainIdx),'RandomForest',true,...
            'nvartosample',d,'NWorkers',NWorkers,'Stratified',Stratified);
        Predictions = oobpredict(rf,X(TrainIdx,:),Ystr(TrainIdx));
        Lhat.rf{k}(j) = oob_error(Predictions,Ystr(TrainIdx),'last');
    end
end

save([rerfPath 'RandomerForest/Results/mnist_378_vary_n.mat'],'ns',...
    'ntrees','Fs','Lhat')
