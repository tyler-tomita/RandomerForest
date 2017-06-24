close all
clear
clc

rng(1);

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

ntrain = 100;
ntest = 10000;
ns = [10 20 50];
ih = 20;
iw = ih;

Xtrain_image = zeros(ih,iw,ntrain);
Ytrain = zeros(ntrain,1);
Xtest_image = zeros(ih,iw,ntest);
Ytest = zeros(ntest,1);

%Class 0
for i = 1:ntrain/2
    rowidx = randperm(ih,5);
    Xtrain_image(rowidx,:,i) = 1;
end

for i = 1:ntest/2
    rowidx = randperm(ih,5);
    Xtest_image(rowidx,:,i) = 1;
end

%Class 1
for i = ntrain/2+1:ntrain
    colidx = randperm(iw,5);
    Xtrain_image(:,colidx,i) = 1;
end
Ytrain(ntrain/2+1:end) = 1;

for i = ntest/2+1:ntest
    colidx = randperm(iw,5);
    Xtest_image(:,colidx,i) = 1;
end
Ytest(ntest/2+1:end) = 1;


NewOrdering = randperm(ntrain);
Xtrain_image = Xtrain_image(:,:,NewOrdering);
Ytrain = Ytrain(NewOrdering);
Labels = unique(Ytrain);

NewOrdering = randperm(ntest);
Xtest_image = Xtest_image(:,:,NewOrdering);
Ytest = Ytest(NewOrdering);

ntrials = 10;

for k = 1:length(ns)
        nTrain = ns(k);
    
    for trial = 1:ntrials

        Idx = [];
        for l = 1:length(Labels)
            Idx = [Idx randsample(find(Ytrain==Labels(l)),round(nTrain/length(Labels)))'];
        end
        TrainIdx{k}(trial,:) = Idx(randperm(length(Idx)));
    end
end

save([rerfPath 'RandomerForest/Data/image_stripes_data.mat'],'ns','ntrials',...
    'TrainIdx','Xtrain_image','Ytrain','Xtest_image','Ytest')