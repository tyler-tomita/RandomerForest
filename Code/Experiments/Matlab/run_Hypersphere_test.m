close all
clear
clc

load Random_matrix_adjustment_factor
 
rng(123);

ntrain = 1000;
ntest = 10000;
Class = [1;2;3];
ntrials = 3;
p = 20;
p_prime1 = 1:18;
p_prime2 = 19:20;
R = 2;
Vcube = (2*R)^length(p_prime1);
Vsphere = pi^(length(p_prime1)/2)*R^length(p_prime1)/gamma(length(p_prime1)/2+1);
ntrain_in = round(1/3*ntrain);
ntrain_out = ntrain - ntrain_in;
ntrain_out_adjusted = ceil(1.5*ntrain_out*Vcube/(Vcube-Vsphere));
ntest_in = round(1/3*ntest);
ntest_out = ntest - ntest_in;
ntest_out_adjusted = ceil(1.5*ntest_out*Vcube/(Vcube-Vsphere));

Classifiers = {'rf' 'rerf' 'frc' 'rr_rf'};

Xtrain = zeros(ntrain,p,ntrials);
Ytrain = cell(ntrain,ntrials);
for trial = 1:ntrials
    x_in = zeros(ntrain_in,p);
    x_in(:,p_prime1) = rand_hypersphere(ntrain_in,length(p_prime1),R);
    x_in(:,[p_prime2,p_prime2(end)+1:p]) = rand(ntrain_in,length([p_prime2,p_prime2(end)+1:p]))*2*R - R;
    y_in = ones(ntrain_in,1);
    x_out = rand(ntrain_out_adjusted,p)*2*R - R;
    Radii = sqrt(sum(x_out(:,p_prime1).^2,2));
    x_out(Radii<=R,:) = [];
    x_out = x_out(1:ntrain_out,:);
    y_out = zeros(ntrain_out,1);
    y_out(sum(x_out(:,p_prime2),2)>0) = 2;
    y_out(sum(x_out(:,p_prime2),2)<=0) = 3;
    y_out = y_out(1:ntrain_out);
    Xtrain(:,:,trial) = [x_in;x_out];
    Ytrain(:,trial) = cellstr(num2str([y_in;y_out]));
end

x_in = zeros(ntest_in,p);
x_in(:,p_prime1) = rand_hypersphere(ntest_in,length(p_prime1),R);
x_in(:,[p_prime2,p_prime2(end)+1:p]) = rand(ntest_in,length([p_prime2,p_prime2(end)+1:p]))*2*R - R;
y_in = ones(ntest_in,1);
x_out = rand(ntest_out_adjusted,p)*2*R - R;
Radii = sqrt(sum(x_out(:,p_prime1).^2,2));
x_out(Radii<=R,:) = [];
x_out = x_out(1:ntest_out,:);
y_out = zeros(ntest_out,1);
y_out(sum(x_out(:,p_prime2),2)>0) = 2;
y_out(sum(x_out(:,p_prime2),2)<=0) = 3;
y_out = y_out(1:ntest_out);
Xtest = [x_in;x_out];
Ytest = cellstr(num2str([y_in;y_out]));

if p <= 5
    mtrys = [1:p p^2 p^3];
elseif p > 5 && p <= 20
    mtrys = ceil(p.^[1/4 1/2 3/4 1 2]);
else
    mtrys = [ceil(p.^[1/4 1/2 3/4 1]) 10*p 20*p];
end
mtrys_rf = mtrys(mtrys<=p);

for c = 1:length(Classifiers)
    fprintf('%s start\n',Classifiers{c})

    Params.(Classifiers{c}).nTrees = 1000;
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
        Params.(Classifiers{c}).nmix = 1:20;
    elseif strcmp(Classifiers{c},'frc') || strcmp(Classifiers{c},'frcr') || ...
            strcmp(Classifiers{c},'frcn') || strcmp(Classifiers{c},'frcz')
        Params.(Classifiers{c}).ForestMethod = 'frc';
        Params.(Classifiers{c}).d = mtrys;
        Params.(Classifiers{c}).nmix = [2,10,18];
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
    end
end