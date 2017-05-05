%% Plot benchmark classifier rank distributions

clear
close all
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

Classifiers = {'rf','rerf','rr_rf'};

inPath1 = [rerfPath 'RandomerForest/Results/2017.04.01/Benchmarks/Raw/'];
inPath2 = '~/Benchmarks/Results/R/dat/Raw/';
contents = dir([inPath1 '*.mat']);

AbsoluteError = NaN(length(contents),length(Classifiers));
NormRelativeError = NaN(length(contents),length(Classifiers)-1);
ChanceProb = NaN(length(contents),1);
p = NaN(length(contents),1);
p_lowrank = NaN(length(contents),1);
ntrain = NaN(length(contents),1);

k = 1;

Datasets = cell(length(contents),1);

for i = 1:length(contents)
    Dataset = strsplit(contents(i).name,'.');
    Dataset = Dataset{1};

    load([inPath1 contents(i).name])

    isComplete = true;

    for c = 1:length(Classifiers)
        cl = Classifiers{c};
        if ~strcmp(cl,'xgb')
            if ~isfield(TestError,cl)
                isComplete = false;
            end
        else
            if ~exist([inPath2 Dataset '_testError.dat'])
                isComplete = false;
            end
        end
    end

    if isComplete
        TrainSet = dlmread(['~/Benchmarks/Data/dat/Raw/' Dataset '_train.dat']);
        TestSet = dlmread(['~/Benchmarks/Data/dat/Raw/' Dataset '_test.dat']);
        nClasses = length(unique(TestSet(:,end)));
        ClassCounts = histcounts(TestSet(:,end),nClasses);
        ChanceProb(k) = 1 - max(ClassCounts)/sum(ClassCounts);

        for c = 1:length(Classifiers)
            cl = Classifiers{c};
            
            % select best model for test predictions
            BI = hp_optimize(OOBError.(cl)(1:length(Params.(cl).d)),...
                OOBAUC.(cl)(1:length(Params.(cl).d)));
            BI = BI(end);
            
            if ~strcmp(cl,'xgb')
                AbsoluteError(k,c) = TestError.(cl)(BI);
            else
                AbsoluteError(k,c) = dlmread([inPath2 Dataset '_testError.dat']);
            end
        end
        for c = 1:length(Classifiers)
            if c > 1
                NormRelativeError(k,c-1) = (AbsoluteError(k,c) - ...
                    AbsoluteError(k,1))/ChanceProb(k);
            end
        end
        k = k + 1;
        Datasets{i} = Dataset;
    end
end
Datasets = Datasets(~cellfun(@isempty,Datasets));
NormRelativeError(all(isnan(NormRelativeError),2),:) = [];

[~,srtidx] = sort(NormRelativeError(:,1));

cutoff = 15;

TopDatasets = Datasets(srtidx(1:cutoff));

for i = 1:cutoff
    TrainSet = dlmread(['~/Benchmarks/Data/dat/Raw/' TopDatasets{i} '_train.dat']);
    dlmwrite(['~/RandomerForest/Data/meda_run/inputs/' TopDatasets{i} '.train.features.csv'],TrainSet(:,1:end-1))
    dlmwrite(['~/RandomerForest/Data/meda_run/inputs/' TopDatasets{i} '.train.labels.csv'],TrainSet(:,end))
    TestSet = dlmread(['~/Benchmarks/Data/dat/Raw/' TopDatasets{i} '_test.dat']);
    dlmwrite(['~/RandomerForest/Data/meda_run/inputs/' TopDatasets{i} '.test.features.csv'],TestSet(:,1:end-1))
    dlmwrite(['~/RandomerForest/Data/meda_run/inputs/' TopDatasets{i} '.test.labels.csv'],TestSet(:,end))
end