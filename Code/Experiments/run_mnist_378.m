close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

rng(1);

X = loadMNISTImages('train-images-idx3-ubyte');
X = sparse(X');
Y = loadMNISTLabels('train-labels-idx1-ubyte');
Ystr = cellstr(num2str(Y));
Labels = [3,7,8];
[n,p] = size(X);
ih = sqrt(p);
iw = ih;

load mnist_378_data

load Random_matrix_adjustment_factor

ntrees = 500;
Stratified = true;
NWorkers = 24;


for k = 1:length(ns)
        nTrain = ns(k);
        fprintf('\nn = %d\n',nTrain)
        
        Lhat.srerf{k} = NaN(ntrials,7);
%         Lhat.srerfd{k} = NaN(ntrials,7);
        Lhat.rerf{k} = NaN(ntrials,7);
        Lhat.rf{k} = NaN(ntrials,5);
    
    for trial = 1:ntrials
        fprintf('\ntrial = %d\n',trial)

        %% Structured RerF %%

        fprintf('\nStructured RerF\n')

        ds = [ceil(p.^[0 1/4 1/2 3/4 1]) 10*p 20*p];

        for j = 1:length(ds)
            d = ds(j);
            fprintf('d = %d\n',d)

            srerf = rpclassificationforest(ntrees,X(TrainIdx{k}(trial,:),:),Ystr(TrainIdx{k}(trial,:)),...
                'Image',true,'ih',ih,'iw',iw,'nvartosample',d,'NWorkers',...
                NWorkers,'Stratified',Stratified);
            Predictions = oobpredict(srerf,X(TrainIdx{k}(trial,:),:),Ystr(TrainIdx{k}(trial,:)));
            Lhat.srerf{k}(trial,j) = oob_error(Predictions,Ystr(TrainIdx{k}(trial,:)),'last');
        end
        
%         %% Structured RerF w/ mean difference %%
% 
%         fprintf('\nStructured RerF(d)\n')
% 
%         ds = [ceil(p.^[0 1/4 1/2 3/4 1]) 10*p 20*p];
% 
%         for j = 1:length(ds)
%             d = ds(j);
%             fprintf('d = %d\n',d)
% 
%             srerfd = rpclassificationforest(ntrees,X(TrainIdx{k}(trial,:),:),Ystr(TrainIdx{k}(trial,:)),...
%                 'Image',true,'ih',ih,'iw',iw,'mdiff','node',...
%                 'nvartosample',d,'NWorkers',NWorkers,'Stratified',Stratified);
%             Predictions = oobpredict(srerfd,X(TrainIdx{k}(trial,:),:),Ystr(TrainIdx{k}(trial,:)));
%             Lhat.srerfd{k}(trial,j) = oob_error(Predictions,Ystr(TrainIdx{k}(trial,:)),'last');
%         end

        %% RerF %%
        fprintf('\nRerF\n')

        ds = [ceil(p.^[0 1/4 1/2 3/4 1]) 10*p 20*p];

        for j = 1:length(ds)
            d = ds(j);

            fprintf('d = %d\n',d)

            dprime = ceil(d^(1/interp1(ps,slope,p)));

            rerf = rpclassificationforest(ntrees,X(TrainIdx{k}(trial,:),:),Ystr(TrainIdx{k}(trial,:)),'sparsemethod',...
                'sparse-adjusted','nvartosample',d,'dprime',dprime,'NWorkers',NWorkers,...
                'Stratified',Stratified);
            Predictions = oobpredict(rerf,X(TrainIdx{k}(trial,:),:),Ystr(TrainIdx{k}(trial,:)));
            Lhat.rerf{k}(trial,j) = oob_error(Predictions,Ystr(TrainIdx{k}(trial,:)),'last');
        end

        %% RF %%
        fprintf('\nRandom Forest\n')

        ds = ceil(p.^[0 1/4 1/2 3/4 1]);

        for j = 1:length(ds)
            d = ds(j);

            fprintf('d = %d\n',d)

            rf = rpclassificationforest(ntrees,X(TrainIdx{k}(trial,:),:),Ystr(TrainIdx{k}(trial,:)),'RandomForest',true,...
                'nvartosample',d,'NWorkers',NWorkers,'Stratified',Stratified);
            Predictions = oobpredict(rf,X(TrainIdx{k}(trial,:),:),Ystr(TrainIdx{k}(trial,:)));
            Lhat.rf{k}(trial,j) = oob_error(Predictions,Ystr(TrainIdx{k}(trial,:)),'last');
        end
    end
end

save([rerfPath 'RandomerForest/Results/mnist_378_vary_n.mat'],'ns',...
    'ntrees','Lhat')