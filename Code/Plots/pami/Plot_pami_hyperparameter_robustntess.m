%% Plot benchmark classifier rank distributions

clear
close all
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);
rerfPath = '~/';

load('purple2green')
ColorMap = interpolate_colormap(ColorMap(round(size(ColorMap,1)/2):end,:),64,false);
Colors.rf = ColorMap(3,:);
Colors.rerf =  ColorMap(9,:);
Marker.rf = 'o';
Marker.rerf = 'x';
LineWidth = 2;
MarkerSize = 8;
FontSize = .2;
axWidth = 2;
axHeight = 2;
cbWidth = axWidth;
cbHeight = 0.15;
axBottom = FontSize*5*ones(1,5);
axLeft = fliplr([FontSize*9+axHeight*4,FontSize*8+axHeight*3,...
    FontSize*7+axHeight*2,FontSize*6+axHeight,...
    FontSize*5]);
legWidth = axWidth*0.75;
legHeight = axHeight;
legLeft = axLeft(end) + axWidth;
legBottom = axBottom(end);
figWidth = legLeft + legWidth;
figHeight = axBottom(1) + axHeight + FontSize*2.5;

% fig = figure;
% fig.Units = 'inches';
% fig.PaperUnits = 'inches';
% fig.Position = [0 0 figWidth figHeight];
% fig.PaperPosition = [0 0 figWidth figHeight];
% fig.PaperSize = [figWidth figHeight];

Classifiers = {'rerfp','frc'};

inPath1 = [rerfPath 'RandomerForest/Results/2017.05.27/Benchmarks/'];
% inPath2 = '~/Benchmarks/Results/R/dat/Raw/';
contents = dir([inPath1 '*.mat']);

AbsoluteError = NaN(length(contents),4,2);
ChanceProb = NaN(length(contents),1);
nClasses = NaN(length(contents),1);
% ModelVariance = NaN(length(contents),length(Classifiers));
dError = NaN(length(contents),4);

DatasetNames = importdata('~/Benchmarks/Data/uci/Names.txt');

% S = load('~/Benchmarks/Results/Benchmark_untransformed.mat');

k = 1;

for i = 1:length(contents)
    fprintf('Dataset %d\n',i)
    Dataset = strsplit(contents(i).name,'_2017_05_27');
    Dataset = Dataset{1};
    DatasetIdx = find(strcmp(Dataset,DatasetNames));

    load([inPath1 contents(i).name])

    isComplete = true;

    for c = 1:length(Classifiers)
        cl = Classifiers{c};
        if ~strcmp(cl,'xgb')
            if ~isfield(TestError,cl)
                isComplete = false;
            end
        end
    end

    if isComplete
        TrainSet = dlmread(['~/Benchmarks/Data/uci/processed/' Dataset '.train.csv']);
        TestSet = dlmread(['~/Benchmarks/Data/uci/processed/' Dataset '.test.csv']);
        nClasses(k) = length(unique(TestSet(:,end)));
        ClassCounts = histcounts(TestSet(:,end),nClasses(k));
        ChanceProb(k) = 1 - max(ClassCounts)/sum(ClassCounts);

        for c = 1:length(Classifiers)
            cl = Classifiers{c};
            if strcmp(cl,'frc')
                SparseParam = 'L';
            elseif strcmp(cl,'rerfp')
                SparseParam = 'lambda';
            end
            ModelError.(cl) = zeros(1,length(Params.(cl).(SparseParam)));
            for j = 1:length(Params.(cl).(SparseParam))
                BI = hp_optimize(OOBError.(cl)(end,(1:length(Params.(cl).d)) + (j-1)*length(Params.(cl).(SparseParam))),...
                    OOBAUC.(cl)(end,(1:length(Params.(cl).d)) + (j-1)*length(Params.(cl).(SparseParam))));
                BI = BI(randperm(length(BI),1));
                BI = BI + (j-1)*length(Params.(cl).(SparseParam));
                ModelError.(cl)(j) = TestError.(cl)(BI);
                AbsoluteError(k,j,c) = TestError.(cl)(BI);
            end
%             ModelVariance(k,c) = max(min(reshape(TestError.(cl),length(Params.(cl).d),length(Params.(cl).(SparseParam))))) - ...
%                 min(min(reshape(TestError.(cl),length(Params.(cl).d),length(Params.(cl).(SparseParam)))));
%             ModelVariance(k,c) = std(min(reshape(TestError.(cl),length(Params.(cl).d),length(Params.(cl).(SparseParam)))));
        end
        dError(k,1:length(ModelError.frc)) = (ModelError.frc - ModelError.rerfp)/ChanceProb(k);
        k = k + 1;
    end
end

AbsoluteError(k,c) = TestError.(cl)(BI);
% ModelVariance(all(isnan(ModelVariance),2),:) = [];
dError(all(isnan(dError),2),:) = [];

figure;
k = 1;

BinEdges = [-1,-0.2,-0.1,-0.05:0.01:-0.01,-0.005,0,0,0.005,0.01:0.01:0.05,0.1,0.2,1];

Counts = zeros(length(BinEdges)-1,length(Classifiers)-1);

for c = 1:4
    Counts(:,c) = histcounts(dError(:,c),BinEdges)';
end
Counts(length(BinEdges)/2,:) = sum(dError==0);
Counts(length(BinEdges)/2+1,:) = Counts(length(BinEdges)/2+1,:) - Counts(length(BinEdges)/2,:);
Fractions = Counts./repmat(sum(Counts),size(Counts,1),1);

ax = axes;
YTLabel = {'L = 2','L = 3','L = 4','L = 5'};

h = heatmap(Fractions',cellstr(num2str(BinEdges')),YTLabel,ColorMap,...
    true,'horizontal');
xlabel({'(Error(F-RC) - Error(RerFp))/Chance'})
title('New Benchmark Datasets')

hold on
for c = 2:4
    plot(h.XLim,[c-0.5,c-0.5],'-k','LineWidth',LineWidth)
end
ax.XTick = [0.5,10,19.5];
ax.XTickLabel = {'-1';'0';'1'};
ax.XTickLabelRotation = 0;
ax.TickLength = [0 0];
ax.LineWidth = LineWidth;


% figure;
% k = 1;
% ax(k) = axes;
% hold on
% plot(sqrt(ModelVariance(:,1)),sqrt(ModelVariance(:,2)),'.','MarkerSize',MarkerSize,...
%     'MarkerSize',14)
% 
% ax(k).Box = 'off';
% xlabel('Model Variance of RerFp')
% ylabel('Model Variance of F-RC')
% title('Updated Benchmark Datasets')
% minx = 0;
% maxx = max(sqrt(ModelVariance(:,1)));
% miny = 0;
% maxy = max(sqrt(ModelVariance(:,2)));
% maxx = min(maxx,maxy);
% maxy = maxx;
% plot([minx maxx],[miny maxy],'k--','LineWidth',2)
save_fig(gcf,[rerfPath 'RandomerForest/Figures/pami/pami_benchmark_hyperparameter_robustness_2'],{'fig','pdf','png'})