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

Xtrain = X(TrainIdx{1}(1,:),:);
Ytrain = Y(TrainIdx{1}(1,:));
Ytrain_str = Ystr(TrainIdx{1}(1,:));

clear X_image

ntrees = 2000;
Stratified = true;
NWorkers = 24;

nTrain = ns(1);

ds = ceil(p.^[0 1/4 1/2 3/4 1]);

BaggedError = NaN(ntrees,length(ds));
for i = 1:length(ds)
    d = ds(i);
    
    fprintf('d = %d\n',d)
    
    srerf = rpclassificationforest(ntrees,Xtrain,Ytrain_str,...
        'Image','on','ih',ih,'iw',iw,'nvartosample',d,'NWorkers',...
        NWorkers,'Stratified',Stratified);
    Predictions = oobpredict(srerf,Xtrain,Ytrain_str);
    BaggedError(:,i) = oob_error(Predictions,Ytrain_str,'every');
end

save([rerfPath 'RandomerForest/Results/image_shapes_number_of_trees.mat'],...
    'BaggedError')