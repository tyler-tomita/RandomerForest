%% Compute posterior heat map for simple "oblique" example 

%% Initialize parameters
close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

rng(1);

Classifiers = {'rf' 'rerf' 'frc' 'rr_rf'};

load Oblique_example_data

xmin = -1;
xmax = 1;
ymin = xmin;
ymax = xmax;
npoints = 50;
[xgv,ygv] = meshgrid(linspace(xmin,xmax,npoints),linspace(ymin,ymax,npoints));
Xpost = xgv(:);
Ypost = ygv(:);

mtrys = [1:p p^2];
mtrys_rf = 1:p;

for i = 1:length(ns)
    [ntrain,p] = size(Xtrain{i});
    fprintf('n = %d\n',ntrain)

    for c = 1:length(Classifiers)
        fprintf('%s start\n',Classifiers{c})

        if strcmp(Classifiers{c},'rf')
            Params.(Classifiers{c}).ForestMethod = 'rf';
            Params.(Classifiers{c}).d = mtrys_rf;
        elseif strcmp(Classifiers{c},'rerf')
            Params.(Classifiers{c}).ForestMethod = 'rerf';
            Params.(Classifiers{c}).RandomMatrix = 'binary';
            Params.(Classifiers{c}).d = mtrys;
        elseif strcmp(Classifiers{c},'frc')
            Params.(Classifiers{c}).ForestMethod = 'rerf';
            Params.(Classifiers{c}).RandomMatrix = 'frc';
            Params.(Classifiers{c}).L = 2;
            Params.(Classifiers{c}).d = mtrys;
        elseif strcmp(Classifiers{c},'rr_rf')
            Params.(Classifiers{c}).ForestMethod = 'rf';
            Params.(Classifiers{c}).Rotate = true;
            Params.(Classifiers{c}).d = mtrys_rf;
        end

        OOBError.(Classifiers{c}) = NaN(1,length(Params.(Classifiers{c}).d));
        OOBAUC.(Classifiers{c}) = NaN(1,length(Params.(Classifiers{c}).d));

        % train classifier

        tic;
        [Forest,~] = ...
            RerF_train(Xtrain{i},...
            Ytrain{i},Params.(Classifiers{c}));

        % select best hyperparameter

        for j = 1:length(Params.(Classifiers{c}).d)
            Scores = rerf_oob_classprob(Forest{j},...
                Xtrain{i},'last');
            Predictions = predict_class(Scores,Forest{j}.classname);
            OOBError.(Classifiers{c})(j) = ...
                misclassification_rate(Predictions,Ytrain{i},...
                false);
            [~,~,~,OOBAUC.(Classifiers{c})(j)] = ...
                perfcurve(Ytrain{i},Scores(:,2),'1');
        end
        BestIdx = hp_optimize(OOBError.(Classifiers{c}),...
            OOBAUC.(Classifiers{c}));
        if length(BestIdx)>1
            BestIdx = BestIdx(end);
        end

        Phats{i}.(Classifiers{c}) = rerf_classprob(Forest{BestIdx},...
            [Xpost,Ypost],'last');

        save([rerfPath 'RandomerForest/Results/2017.02.15/Oblique_example_posteriors.mat'],...
            'Phats','Xpost','Ypost')
        fprintf('%s complete\n',Classifiers{c})
    end
end