close all
clear
clc

rng(1);

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

n = 10000;
ns = [50 100 500 1000];
ih = 32;
iw = ih;
p = ih*iw;

X_image = zeros(ih,iw,n);
[rowidx,colidx] = ind2sub([ih,iw],1:p);
radius = ceil(rand(1,n)*3) + 9;
Y = zeros(n,1);

%Class 0
for i = 1:n/2
    centroid(i,:) = randi(ih-radius(i)*2,1,2) + radius(i);
    x = zeros(ih,iw);
    x(sqrt(sum(([rowidx',colidx'] - repmat(centroid(i,:),p,1)).^2,2)) <= radius(i)) = 1;
    X_image(:,:,i) = x;
end

%Class 1
for i = n/2+1:n
    centroid(i,:) = randi(ih-radius(i)*2,1,2) + radius(i);
    x = zeros(ih,iw);
    x(abs(rowidx-centroid(i,1))<=radius(i)&abs(colidx-centroid(i,2))<=radius(i)) = 1;
    X_image(:,:,i) = x;
end
Y(n/2+1:end) = 1;

NewOrdering = randperm(n);
X_image = X_image(:,:,NewOrdering);
Y = Y(NewOrdering);
Labels = unique(Y);

ntrials = 10;

for k = 1:length(ns)
        nTrain = ns(k);
    
    for trial = 1:ntrials

        Idx = [];
        for l = 1:length(Labels)
            Idx = [Idx randsample(find(Y==Labels(l)),round(nTrain/length(Labels)))'];
        end
        TrainIdx{k}(trial,:) = Idx(randperm(length(Idx)));
    end
end

save([rerfPath 'RandomerForest/Data/image_shapes_data.mat'],'ns','ntrials',...
    'TrainIdx','X_image','Y')
