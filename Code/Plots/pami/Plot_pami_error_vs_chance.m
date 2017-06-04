%% Plot benchmark classifier rank distributions

clear
close all
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

load('purple2green')
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

Classifiers = {'rf','rerf','frc','rr_rf'};

inPath1 = [rerfPath 'RandomerForest/Results/2017.05.27/Benchmarks/'];
% inPath2 = '~/Benchmarks/Results/R/dat/Raw/';
contents = dir([inPath1 '*.mat']);

AbsoluteError = NaN(length(contents),length(Classifiers));
NormalizedRelativeError = NaN(length(contents),length(Classifiers)-1);
ChanceProb = NaN(length(contents),1);
nClasses = NaN(length(contents),1);
d = NaN(length(contents),length(Classifiers));
p = NaN(length(contents),1);
p_lowrank = NaN(length(contents),1);
lambda = cell(length(contents),1);   % eigenvalues of data matrix
ntrain = NaN(length(contents),1);
TraceNorm = NaN(length(contents),1);

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
%         else
%             if ~exist([inPath2 Dataset '_testError.dat'])
%                 isComplete = false;
%             end
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
            if ~strcmp(cl,'xgb')
                BI = hp_optimize(OOBError.(cl)(end,1:length(Params.(cl).d)),...
                    OOBAUC.(cl)(end,1:length(Params.(cl).d)));
                BI = BI(randperm(length(BI),1));
                AbsoluteError(k,c) = TestError.(cl)(BI);
%             else
%                 AbsoluteError(k,c) = dlmread([inPath2 Dataset '_testError.dat']);
            end
        end
        k = k + 1;
    end
end

AbsoluteError(all(isnan(AbsoluteError),2),:) = [];
ChanceProb(all(isnan(ChanceProb),2),:) = [];

figure;
k = 1;
ax(k) = axes;
hold on
for c = 1:2
    cl = Classifiers{c};
    plot(ChanceProb,AbsoluteError(:,c),Marker.(cl),'MarkerSize',MarkerSize,...
        'MarkerEdgeColor',Colors.(cl),'MarkerSize',10)
end
% logx = log10(p);
% linmod = fitlm(logx,NormalizedRelativeError(:,1),'linear');
% beta = linmod.Coefficients.Estimate;
% yfit = beta(1) + beta(2)*logx;
% hold on
% plot(logx,yfit,'k','LineWidth',2)
% plot(logx,yfit,'LineStyle','none','Marker','none')

% ax(k).XScale = 'log';
ax(k).Box = 'off';
xlabel('Chance Error')
ylabel('Classifier Error')
title('Updated Benchmark Datasets')
plot([0 1],[0 1],'k--','LineWidth',2)
legend('RF','RerF')
save_fig(gcf,[rerfPath 'RandomerForest/Figures/pami/pami_benchmark_chance_vs_error'],{'fig','pdf','png'})