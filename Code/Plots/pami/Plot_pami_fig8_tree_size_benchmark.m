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
FontSize = .2;
axWidth = 2;
axHeight = 2;
axLeft = FontSize*4;
axBottom = FontSize*4;
legWidth = axWidth;
legHeight = axHeight;
legLeft = axLeft + axWidth + FontSize;
legBottom = axBottom;
figWidth = legLeft + legWidth;
figHeight = axBottom + axHeight + FontSize*2.5;

fig = figure;
fig.Units = 'inches';
fig.PaperUnits = 'inches';
fig.Position = [0 0 figWidth figHeight];
fig.PaperPosition = [0 0 figWidth figHeight];
fig.PaperSize = [figWidth figHeight];

Classifiers = {'rf','rerf','rr_rf','xgb'};

inPath1 = [rerfPath 'RandomerForest/Results/2017.04.01/Benchmarks/Raw/'];
inPath2 = '~/Benchmarks/Results/R/dat/Raw/';
contents = dir([inPath1 '*.mat']);

AbsoluteError = NaN(length(contents),length(Classifiers));
NormError = NaN(length(contents),length(Classifiers));
ChanceProb = NaN(length(contents),1);
MeanDepth = NaN(length(contents),length(Classifiers));

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
        [ntrain,p] = size(TrainSet(:,1:end-1));
        nClasses = length(unique(TrainSet(:,end)));
        ClassCounts = histcounts(TrainSet(:,end),nClasses);
        ChanceProb(k) = 1 - max(ClassCounts)/sum(ClassCounts);

        for c = 1:length(Classifiers)
            cl = Classifiers{c};
            if ~strcmp(cl,'xgb')
                AbsoluteError(k,c) = TestError.(cl);
                MeanDepth(k,c) = mean(Depth.(cl)(:,BestIdx.(cl)));
            else
                AbsoluteError(k,c) = dlmread([inPath2 Dataset '_testError.dat']);
                MeanDepth(k,c) = dlmread([inPath2 Dataset '_depth.dat']);
            end
            NormError(k,c) = AbsoluteError(k,c)/ChanceProb(k);
        end
        k = k + 1;
    end
end

NormError(all(isnan(NormError),2),:) = [];
MeanDepth(all(isnan(MeanDepth),2),:) = [];

ax = axes;
hold on
for c = 1:length(Classifiers)
    cl = Classifiers{c};
    plot(MeanDepth(:,c),NormError(:,c),'.','MarkerSize',MarkerSize,...
        'Color',Colors.(cl))
end

xlabel('Mean Depth')
ylabel('Error Rate')
title('Raw Benchmarks')

ax.LineWidth = LineWidth;
ax.FontUnits = 'inches';
ax.FontSize = FontSize;
ax.Units = 'inches';
ax.Position = [axLeft axBottom axWidth axHeight];
ax.XScale = 'log';
ax.YScale = 'log';

lh = legend('RF','RerF','RR-RF','XGBoost');
lh.Box = 'off';
lh.Units = 'inches';
lh.Position = [legLeft legBottom legWidth legHeight];
        
save_fig(gcf,[rerfPath 'RandomerForest/Figures/pami/PAMI_fig8_tree_size_benchmark'],{'fig','pdf','png'})