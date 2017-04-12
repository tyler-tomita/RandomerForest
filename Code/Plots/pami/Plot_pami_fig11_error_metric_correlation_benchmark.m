%% Plot benchmark classifier rank distributions

clear
close all
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

load('purple2green')
Colors.rf = ColorMap(2,:);
Colors.rerf= ColorMap(10,:);
Colors.rr_rf = ColorMap(4,:);
Colors.xgb= ColorMap(8,:);
LineWidth = 2;
MarkerSize = 8;
FontSize = 0.175;
axWidth = 1.5;
axHeight = 1.5;
axLeft = [FontSize*4,FontSize*7+axWidth];
axBottom = FontSize*4*ones(1,2);
legWidth = axWidth*0.75;
legHeight = axHeight;
legLeft = axLeft(end) + axWidth;
legBottom = axBottom(end);
figWidth = legLeft + legWidth;
figHeight = axBottom(1) + axHeight + FontSize*2.5;

fig = figure;
fig.Units = 'inches';
fig.PaperUnits = 'inches';
fig.Position = [0 0 figWidth figHeight];
fig.PaperPosition = [0 0 figWidth figHeight];
fig.PaperSize = [figWidth figHeight];

Classifiers = {'rf','rerf','rr_rf'};

inPath1 = [rerfPath 'RandomerForest/Results/2017.04.01/Benchmarks/Raw/'];
inPath2 = '~/Benchmarks/Results/R/dat/Raw/';
contents = dir([inPath1 '*.mat']);

AbsoluteError = NaN(length(contents),length(Classifiers));
RelativeError = NaN(length(contents),length(Classifiers));
ChanceProb = NaN(length(contents),1);
p = NaN(length(contents),1);
p_lowrank = NaN(length(contents),1);
ntrain = NaN(length(contents),1);

k = 1;

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
        [ntrain(k),p(k)] = size(TrainSet(:,1:end-1));
        [coeff,score,lambda] = pca([TrainSet(:,1:end-1);TestSet(:,1:end-1)]);
        VarCutoff = 0.9;
        AboveCutoff = find(cumsum(lambda)/sum(lambda) >= 0.9);
        p_lowrank(k) = AboveCutoff(1);
        nClasses = length(unique(TestSet(:,end)));
        ClassCounts = histcounts(TestSet(:,end),nClasses);
        ChanceProb(k) = 1 - max(ClassCounts)/sum(ClassCounts);

        for c = 1:length(Classifiers)
            cl = Classifiers{c};
            if ~strcmp(cl,'xgb')
                AbsoluteError(k,c) = TestError.(cl)(BestIdx.(cl));
            else
                AbsoluteError(k,c) = dlmread([inPath2 Dataset '_testError.dat']);
            end
%             if c > 1
                RelativeError(k,c) = AbsoluteError(k,c)/ChanceProb(k);
%             end
        end
        k = k + 1;
    end
end

AbsoluteError(all(isnan(RelativeError),2),:) = [];
RelativeError(all(isnan(RelativeError),2),:) = [];
p(isnan(p)) = [];
p_lowrank(isnan(p_lowrank)) = [];
ntrain(isnan(ntrain)) = [];

k = 1;
ax(k) = axes;
hold on
for c = 1:length(Classifiers)
    cl = Classifiers{c};
    plot(p_lowrank./p,RelativeError(:,c),'.','MarkerSize',MarkerSize,...
        'Color',Colors.(cl))
end