close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

rng(1);

load image_shapes_data

[ih,iw,n] = size(X_image);
p = ih*iw;
X = reshape(X_image,p,n)';
Ystr = cellstr(num2str(Y));

clear X_image

load Random_matrix_adjustment_factor

ntrees = 500;
Stratified = true;
NWorkers = 2;

FileID = fopen('~/shapes.out','w');

for k = 1:length(ns)
        nTrain = ns(k);
        fprintf(FileID,'\nn = %d\n',nTrain);
        
        Lhat.srerf{k} = NaN(ntrials,7);
        Lhat.control{k} = NaN(ntrials,7);
        Lhat.rerf{k} = NaN(ntrials,7);
        Lhat.rf{k} = NaN(ntrials,5);
    
    for trial = 1:ntrials
        fprintf(FileID,'\ntrial = %d\n',trial);

        %% Structured RerF %%

        fprintf(FileID,'\nStructured RerF\n');

        ds = [ceil(p.^[0 1/4 1/2 3/4 1]) 10*p 20*p];

        for j = 1:length(ds)
            d = ds(j);
            fprintf(FileID,'d = %d\n',d);

            srerf{j} = rpclassificationforest(ntrees,X(TrainIdx{k}(trial,:),:),Ystr(TrainIdx{k}(trial,:)),...
                'Image','on','ih',ih,'iw',iw,'nvartosample',d,'NWorkers',...
                NWorkers,'Stratified',Stratified);
            Predictions = oobpredict(srerf{j},X(TrainIdx{k}(trial,:),:),Ystr(TrainIdx{k}(trial,:)));
            Lhat.srerf{k}(trial,j) = oob_error(Predictions,Ystr(TrainIdx{k}(trial,:)),'last');
        end
        
        [~,BestIdx] = min(Lhat.srerf{k}(trial,:),[],2);
        Predictions = predict(srerf{BestIdx},X(Test,:));
        TestError.srerf{k}(trial) = misclassification_rate(Predictions,Ystr(Test),false);
        
        clear srerf
        
        %% RerF controlled for density%%
        fprintf(FileID,'\nRerF control\n');

        ds = [ceil(p.^[0 1/4 1/2 3/4 1]) 10*p 20*p];

        for j = 1:length(ds)
            d = ds(j);

            fprintf(FileID,'d = %d\n',d);

            control{j} = rpclassificationforest(ntrees,X(TrainIdx{k}(trial,:),:),Ystr(TrainIdx{k}(trial,:)),...
                'Image','control','ih',ih,'iw',iw,'nvartosample',d,'NWorkers',...
                NWorkers,'Stratified',Stratified);
            Predictions = oobpredict(control{j},X(TrainIdx{k}(trial,:),:),Ystr(TrainIdx{k}(trial,:)));
            Lhat.control{k}(trial,j) = oob_error(Predictions,Ystr(TrainIdx{k}(trial,:)),'last');
        end
        
        [~,BestIdx] = min(Lhat.control{k}(trial,:),[],2);
        Predictions = predict(control{BestIdx},X(Test,:));
        TestError.control{k}(trial) = misclassification_rate(Predictions,Ystr(Test),false);
        
        clear control

        %% RerF %%
        fprintf(FileID,'\nRerF\n');

        ds = [ceil(p.^[0 1/4 1/2 3/4 1]) 10*p 20*p];

        for j = 1:length(ds)
            d = ds(j);

            fprintf(FileID,'d = %d\n',d);

            dprime = ceil(d^(1/interp1(ps,slope,p,'linear','extrap')));

            rerf{j} = rpclassificationforest(ntrees,X(TrainIdx{k}(trial,:),:),Ystr(TrainIdx{k}(trial,:)),'sparsemethod',...
                'sparse-adjusted','nvartosample',d,'dprime',dprime,'NWorkers',NWorkers,...
                'Stratified',Stratified);
            Predictions = oobpredict(rerf{j},X(TrainIdx{k}(trial,:),:),Ystr(TrainIdx{k}(trial,:)));
            Lhat.rerf{k}(trial,j) = oob_error(Predictions,Ystr(TrainIdx{k}(trial,:)),'last');
        end
        
        [~,BestIdx] = min(Lhat.rerf{k}(trial,:),[],2);
        Predictions = predict(rerf{BestIdx},X(Test,:));
        TestError.rerf{k}(trial) = misclassification_rate(Predictions,Ystr(Test),false);
        
        clear rerf

        %% RF %%
        fprintf(FileID,'\nRandom Forest\n');

        ds = ceil(p.^[0 1/4 1/2 3/4 1]);

        for j = 1:length(ds)
            d = ds(j);

            fprintf(FileID,'d = %d\n',d);

            rf{j} = rpclassificationforest(ntrees,X(TrainIdx{k}(trial,:),:),Ystr(TrainIdx{k}(trial,:)),'RandomForest',true,...
                'nvartosample',d,'NWorkers',NWorkers,'Stratified',Stratified);
            Predictions = oobpredict(rf{j},X(TrainIdx{k}(trial,:),:),Ystr(TrainIdx{k}(trial,:)));
            Lhat.rf{k}(trial,j) = oob_error(Predictions,Ystr(TrainIdx{k}(trial,:)),'last');
        end
        
        [~,BestIdx] = min(Lhat.rf{k}(trial,:),[],2);
        Predictions = predict(rf{BestIdx},X(Test,:));
        TestError.rf{k}(trial) = misclassification_rate(Predictions,Ystr(Test),false);
        
        clear rf
        
        save([rerfPath 'RandomerForest/Results/image_shapes.mat'],...
            'ntrees','Lhat','TestError')
    end
end