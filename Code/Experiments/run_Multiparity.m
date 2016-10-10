close all
clear
clc

load Random_matrix_adjustment_factor
 
rng(123);

ntrain = 1000;
ntest = 10000;
ntrials = 1;
p = 50;
p_prime1 = 1:3;
p_prime2 = 4:p;
Sigma1 = 1/64*ones(1,length(p_prime1));
Sigma2 = ones(1,length(p_prime2));
Lambda = [0.33 0.67];
Mu2 = [-0.4*ones(1,length(p_prime2));0.4*ones(1,length(p_prime2))];
obj = gmdistribution(Mu2,Sigma2,Lambda);

Classifiers = {'rerf' 'frc'};

Xtrain = zeros(ntrain,p,ntrials);
Ytrain = cell(ntrain,ntrials);
for trial = 1:ntrials
    Mu1 = sparse(ntrain,length(p_prime1));
    xparity = zeros(ntrain,length(p_prime1));
    for i = 1:ntrain
        Mu1(i,:) = binornd(1,0.5,1,length(p_prime1));
        xparity(i,:) = mvnrnd(Mu1(i,:),Sigma1);
    end
    nOnes1 = sum(Mu1,2);
    [xgaussian,idx] = random(obj,ntrain);
    Xtrain(:,:,trial) = [xparity,xgaussian];
%     nOnes2 = sum(Mu(:,p_prime2),2);
    y = zeros(ntrain,1);
    y(mod(nOnes1,2)==0 & idx==1) = 1;
    y(mod(nOnes1,2)==1 & idx==1) = 2;
    y(mod(nOnes1,2)==0 & idx==2) = 3;
    y(mod(nOnes1,2)==1 & idx==2) = 3;
    Ytrain(:,trial) = cellstr(num2str(y));
end

Mu1 = sparse(ntest,length(p_prime1));
xparity = zeros(ntest,length(p_prime1));
for i = 1:ntest
    Mu1(i,:) = binornd(1,0.5,1,length(p_prime1));
    xparity(i,:) = mvnrnd(Mu1(i,:),Sigma1);
end
nOnes1 = sum(Mu1,2);
[xgaussian,idx] = random(obj,ntest);
Xtest = [xparity,xgaussian];
%     nOnes2 = sum(Mu(:,p_prime2),2);
y = zeros(ntest,1);
y(mod(nOnes1,2)==0 & idx==1) = 1;
y(mod(nOnes1,2)==1 & idx==1) = 2;
y(mod(nOnes1,2)==0 & idx==2) = 3;
y(mod(nOnes1,2)==1 & idx==2) = 3;
Ytest = cellstr(num2str(y));

if p <= 5
    mtrys = [1:p p^2];
elseif p > 5 && p <= 50
    mtrys = ceil(p.^[1/4 1/2 3/4 1 2]);
else
    mtrys = [ceil(p.^[1/4 1/2 3/4 1]) 10*p 20*p];
end
mtrys_rf = mtrys(mtrys<=p);

for c = 1:length(Classifiers)
    fprintf('%s start\n',Classifiers{c})

    Params.(Classifiers{c}).nTrees = 500;
    Params.(Classifiers{c}).Stratified = true;
    Params.(Classifiers{c}).NWorkers = 2;
    if strcmp(Classifiers{c},'rfr') || strcmp(Classifiers{c},...
            'rerfr') || strcmp(Classifiers{c},'frcr') || ...
            strcmp(Classifiers{c},'rr_rfr')
        Params.(Classifiers{c}).Rescale = 'rank';
    elseif strcmp(Classifiers{c},'rfn') || strcmp(Classifiers{c},...
            'rerfn') || strcmp(Classifiers{c},'frcn') || ...
            strcmp(Classifiers{c},'rr_rfn')
        Params.(Classifiers{c}).Rescale = 'normalize';
    elseif strcmp(Classifiers{c},'rfz') || strcmp(Classifiers{c},...
            'rerfz') || strcmp(Classifiers{c},'frcz') || ...
            strcmp(Classifiers{c},'rr_rfz')
        Params.(Classifiers{c}).Rescale = 'zscore';
    else
        Params.(Classifiers{c}).Rescale = 'off';
    end
    if strcmp(Classifiers{c},'rerfd')
        Params.(Classifiers{c}).mdiff = 'node';
    else
        Params.(Classifiers{c}).mdiff = 'off';
    end
    if strcmp(Classifiers{c},'rf') || strcmp(Classifiers{c},'rfr')...
            || strcmp(Classifiers{c},'rfn') || strcmp(Classifiers{c},'rfz') || ...
            strcmp(Classifiers{c},'rr_rf') || strcmp(Classifiers{c},'rr_rfr') || ...
            strcmp(Classifiers{c},'rr_rfn') || strcmp(Classifiers{c},'rr_rfz')
        Params.(Classifiers{c}).ForestMethod = 'rf';
        Params.(Classifiers{c}).d = mtrys_rf;
        Params.(Classifiers{c}).nmix = 1;
    elseif strcmp(Classifiers{c},'rerf') || strcmp(Classifiers{c},'rerfr')...
            || strcmp(Classifiers{c},'rerfn') || strcmp(Classifiers{c},'rerfz') || ...
            strcmp(Classifiers{c},'rerfd')
        Params.(Classifiers{c}).ForestMethod = 'uniform-nnzs';
        Params.(Classifiers{c}).d = mtrys;
        for j = 1:length(Params.(Classifiers{c}).d)
            Params.(Classifiers{c}).dprime(j) = ...
                ceil(Params.(Classifiers{c}).d(j)^(1/interp1(ps,...
                slope,p)));
        end
        Params.(Classifiers{c}).nmix = 1:p;
    elseif strcmp(Classifiers{c},'frc') || strcmp(Classifiers{c},'frcr') || ...
            strcmp(Classifiers{c},'frcn') || strcmp(Classifiers{c},'frcz')
        Params.(Classifiers{c}).ForestMethod = 'frc';
        Params.(Classifiers{c}).d = mtrys;
        Params.(Classifiers{c}).nmix = [3,5,10,25,50];
    end
    if strcmp(Classifiers{c},'rr_rf') || strcmp(Classifiers{c},'rr_rfr') || ...
            strcmp(Classifiers{c},'rr_rfn') || strcmp(Classifiers{c},'rr_rfz')
        Params.(Classifiers{c}).Rotate = true;
    end
    
    if strcmp(Params.(Classifiers{c}).ForestMethod,'frc')
        OOBError.(Classifiers{c}) = ...
            NaN(ntrials,length(Params.(Classifiers{c}).d)*length(Params.(Classifiers{c}).nmix));
        OOBAUC.(Classifiers{c}) = ...
            NaN(ntrials,length(Params.(Classifiers{c}).d)*length(Params.(Classifiers{c}).nmix));
        TrainTime.(Classifiers{c}) = ...
            NaN(ntrials,length(Params.(Classifiers{c}).d)*length(Params.(Classifiers{c}).nmix));
    else
        OOBError.(Classifiers{c}) = ...
            NaN(ntrials,length(Params.(Classifiers{c}).d));
        OOBAUC.(Classifiers{c}) = ...
            NaN(ntrials,length(Params.(Classifiers{c}).d));
        TrainTime.(Classifiers{c}) = ...
            NaN(ntrials,length(Params.(Classifiers{c}).d));
    end
    
    for trial = 1:ntrials
        fprintf('trial %d\n',trial)

        [Forest,~,TrainTime.(Classifiers{c})(trial,:)] = ...
            RerF_train(Xtrain(:,:,trial),...
            Ytrain(:,trial),Params.(Classifiers{c}));

        % select best hyperparameter

        for j = 1:length(Forest)
            if ~isempty(Forest{j})
                Scores = rerf_oob_classprob(Forest{j},...
                    Xtrain(:,:,trial),'last');
                Predictions = predict_class(Scores,Forest{j}.classname);
                OOBError.(Classifiers{c})(trial,j) = ...
                    misclassification_rate(Predictions,Ytrain(:,trial),...
                    false);
                if size(Scores,2) > 2
                    Yb = binarize_labels(Ytrain(:,trial),Forest{j}.classname);
                    [~,~,~,OOBAUC.(Classifiers{c})(trial,j)] = ...
                        perfcurve(Yb(:),Scores(:),'1');
                else
                    [~,~,~,OOBAUC.(Classifiers{c})(trial,j)] = ...
                        perfcurve(Ytrain(:,trial),Scores(:,2),'1');
                end
            end
        end
        NotEmptyIdx = find(~isnan(OOBError.(Classifiers{c})(trial,:)));
        BestIdx = hp_optimize(OOBError.(Classifiers{c})(trial,NotEmptyIdx),...
            OOBAUC.(Classifiers{c})(trial,NotEmptyIdx));
        if length(BestIdx)>1
            BestIdx = BestIdx(end);
        end
        BestIdx = NotEmptyIdx(BestIdx);

        if strcmp(Forest{BestIdx}.Rescale,'off')
            Scores = rerf_classprob(Forest{BestIdx},Xtest,'last');
        else
            Scores = rerf_classprob(Forest{BestIdx},Xtest,...
                'last',Xtrain(:,:,trial));
        end
        Predictions = predict_class(Scores,Forest{BestIdx}.classname);
        TestError.(Classifiers{c})(trial) = misclassification_rate(Predictions,...
            Ytest,false);

        clear Forest
        
        save('~/Multiparity_test.mat','Params','OOBError','OOBAUC',...
            'TestError','TrainTime')
    end
end