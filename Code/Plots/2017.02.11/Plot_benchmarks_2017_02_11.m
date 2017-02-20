%% Plot distributions of errors for benchmark data

clear
close all
clc

load('purple2green')
ColorMap = interpolate_colormap(ColorMap(round(size(ColorMap,1)/2):end,:),64,false);
% ColorMap = flipud(ColorMap);
% Colormap = parula;
FontSize = 12;
LineWidth = 2;

% LineWidth = 2;
% LineWidth_box = 4;
% LineWidth_whisker = 1.5;
% MarkerSize = 6;
% FontSize = .18;
% axWidth = 1.5;
% axHeight = 1.4;
% axLeft = FontSize*5*ones(1,5);
% axBottom = [FontSize*12+axHeight*4,FontSize*7+axHeight*3,...
%     FontSize*5+axHeight*2,FontSize*3+axHeight,...
%     FontSize];
% figWidth = axLeft(end) + axWidth + FontSize;
% figHeight = axBottom(1) + axHeight + FontSize*2;

% fig = figure;
% fig.Units = 'inches';
% fig.PaperUnits = 'inches';
% fig.Position = [0 0 figWidth figHeight];
% fig.PaperPosition = [0 0 figWidth figHeight];
% fig.PaperSize = [figWidth figHeight];

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);


BinEdges = [-1,-0.2,-0.1,-0.05:0.005:0,0:0.005:0.05,0.1,0.2,1];
% BinEdges = [-1:0.1:0,0:0.1:1];

load aggregated_results_2017_02_11
nDatasets = length(Results);

Classifiers = fieldnames(Results(1).TestError);
Classifiers(strcmp(Classifiers,'rr_rf')) = [];

ChanceProb = NaN(nDatasets,1);
AbsoluteError = NaN(nDatasets,length(Classifiers));
NormRelativeError = NaN(nDatasets,length(Classifiers)-1);

for j = 1:nDatasets
    DatasetName = Results(j).Name(1:regexp(Results(j).Name,'\.mat')-1);
    Xtrain = dlmread(['~/Benchmarks/Data/dat/' DatasetName '_train.dat']);
    Ytrain = Xtrain(:,end);
    Xtrain(:,end) = [];
    ClassCounts = histcounts(Ytrain);
    ChanceProb(j) = 1 - max(ClassCounts)/sum(ClassCounts);
    for k = 1:length(Classifiers)
        if strcmp(Classifiers{k},'rerf_big')
            AbsoluteError(j,k) = Results(j).TestError.(Classifiers{k})(Results(j).BestIdx.(Classifiers{k}));
        else
            AbsoluteError(j,k) = Results(j).TestError.(Classifiers{k});
        end
    end
    for k = 1:length(Classifiers)
        if k > 1
            NormRelativeError(j,k-1) = (AbsoluteError(j,k) ...
                - AbsoluteError(j,1))/ChanceProb(j);
        end
    end
end
    
Counts = zeros(length(BinEdges)-1,length(Classifiers)-1);
    
for k = 1:length(Classifiers)-1
%         h = histogram(RelativeError(:,k),BinEdges);
%         Counts(:,k) = h.Values';
    Counts(:,k) = histcounts(NormRelativeError(:,k),BinEdges)';
end
Counts(length(BinEdges)/2,:) = sum(NormRelativeError==0);
Counts(length(BinEdges)/2+1,:) = Counts(length(BinEdges)/2+1,:) - Counts(length(BinEdges)/2,:);

p.rerf = signrank(AbsoluteError(:,1),AbsoluteError(:,2),'tail','right');
p.rerf_big = signrank(AbsoluteError(:,1),AbsoluteError(:,3),'tail','right');

figure;
h = heatmap(flipud(Counts),{['RerF\newlinep = ',num2str(p.rerf,'%0.4f')],['RerF+\newlinep = ',num2str(p.rerf_big,'%0.4f')]},cellstr(num2str(flipud(BinEdges'))),ColorMap,true);
ylabel('($\hat{L}_X - \hat{L}_{RF})/Chance$','Interpreter','latex')
title(sprintf('%d Benchmark Datasets',nDatasets))
colorbar;
h.FontSize = FontSize;
hold on
for k = 2:length(Classifiers)-1
    plot([k-0.5,k-0.5],h.YLim,'-k','LineWidth',LineWidth+1)
end
save_fig(gcf,'~/RandomerForest/Figures/2017.02.11/Benchmark_heatmap_histogram_normalized_by_chance_2017_02_17')