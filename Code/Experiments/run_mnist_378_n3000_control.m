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

ntrees = 500;
Stratified = true;
NWorkers = 24;

FileID = fopen('~/mnist.out','w');

for k = 4:4
        nTrain = ns(k);
        fprintf(FileID,'\nn = %d\n',nTrain);
        
        Lhat.srerf{k} = NaN(ntrials,7);
        Lhat.control{k} = NaN(ntrials,7);
    
    for trial = 1:ntrials
        fprintf(FileID,'\ntrial = %d\n',trial);
        
        %% Structured RerF %%

        fprintf(FileID,'\nStructured RerF\n');

        ds = [ceil(p.^[0 1/4 1/2 3/4 1]) 10*p 20*p];

        for j = 1:length(ds)
            d = ds(j);
            fprintf(FileID,'d = %d\n',d);

            srerf = rpclassificationforest(ntrees,X(TrainIdx{k}(trial,:),:),Ystr(TrainIdx{k}(trial,:)),...
                'Image','on','ih',ih,'iw',iw,'nvartosample',d,'NWorkers',...
                NWorkers,'Stratified',Stratified);
            Predictions = oobpredict(srerf,X(TrainIdx{k}(trial,:),:),Ystr(TrainIdx{k}(trial,:)));
            Lhat.srerf{k}(trial,j) = oob_error(Predictions,Ystr(TrainIdx{k}(trial,:)),'last');
        end        
        
        %% RerF controlled for density%%
        fprintf(FileID,'\nRerF control\n');

        ds = [ceil(p.^[0 1/4 1/2 3/4 1]) 10*p 20*p];

        for j = 1:length(ds)
            d = ds(j);

            fprintf(FileID,'d = %d\n',d);

            control = rpclassificationforest(ntrees,X(TrainIdx{k}(trial,:),:),Ystr(TrainIdx{k}(trial,:)),...
                'Image','control','ih',ih,'iw',iw,'nvartosample',d,'NWorkers',...
                NWorkers,'Stratified',Stratified);
            Predictions = oobpredict(control,X(TrainIdx{k}(trial,:),:),Ystr(TrainIdx{k}(trial,:)));
            Lhat.control{k}(trial,j) = oob_error(Predictions,Ystr(TrainIdx{k}(trial,:)),'last');
        end
    end
end

fclose(FileID);

save([rerfPath 'RandomerForest/Results/mnist_378_n3000_control.mat'],'ns',...
    'ntrees','Lhat')