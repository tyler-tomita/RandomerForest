%% Plot distributions of errors for benchmark data

clear
close all
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

LineWidth = 2;
FontSize = .12;
axWidth = 3.75;
axHeight = 0.75;
axLeft = FontSize*3;
axBottom = FontSize*4;
legWidth = 1;
legHeight = 1;
legLeft = axLeft + axWidth - 4*FontSize;
legBottom = axBottom - (legHeight - axHeight)/2;
figWidth = legLeft + legWidth;
figHeight = axBottom + axHeight + FontSize*2;

fig = figure;
fig.Units = 'inches';
fig.PaperUnits = 'inches';
fig.Position = [0 0 figWidth figHeight];
fig.PaperPosition = [0 0 figWidth figHeight];
fig.PaperSize = [figWidth figHeight];

load aggregated_results_2017_01_16
nDatasets = length(Results);

Classifiers = fieldnames(Results(1).TestError);
Classifiers{end+1} = 'xgb';

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
        if strcmp(Classifiers{k},'xgb')
            fh = fopen(['~/Benchmarks/Results/dat/' DatasetName, '_testError.dat']);
            if fh==-1
                AbsoluteError(j,k) = NaN;
            else
                e = textscan(fh,'%f');
                AbsoluteError(j,k) = e{1};
                fclose(fh);
            end
        else
            AbsoluteError(j,k) = Results(j).TestError.(Classifiers{k});
        end
    end
end

% h = heatmap(flipud(Counts),{'RerF','RR-RF','XGBoost'},cellstr(num2str(flipud(BinEdges'))),ColorMap,true);
% ylabel('($\hat{L}_X - \hat{L}_{RF})/Chance$','Interpreter','latex')
% title('119 Benchmark Datasets')
% colorbar;
% h.FontSize = FontSize;
% hold on
% for k = 2:length(Classifiers)-1
%     plot([k-0.5,k-0.5],h.YLim,'-k','LineWidth',LineWidth+1)
% end


ClRanks = tiedrank(AbsoluteError')';
IntRanks = floor(ClRanks);


RankCounts = NaN(length(Classifiers));
for j = 1:length(Classifiers)
    RankCounts(j,:) = sum(IntRanks==j); 
%     MeanRank.(Classifiers{j}).(Transformations{i}) = mean(IntRanks(:,j));
end

ax = axes;
bar(RankCounts')
Bars = allchild(ax);
for j = 1:length(Bars)
    Bars(j).EdgeColor = 'k';
    Bars(j).BarWidth = 1;
end

ylabel('Frequency')
text(0.5,1.05,'119 Benchmark Datasets','FontSize',12,'FontWeight','bold','Units',...
    'normalized','HorizontalAlignment','center','VerticalAlignment'...
    ,'bottom')
    [lh,objh] = legend('1st','2nd','3rd','4th');
    lh.Box = 'off';
    lh.FontSize = 9;
    lh.Units = 'inches';
    lh.Position = [legLeft legBottom legWidth legHeight];
    BarWidth = (objh(5).Children.YData(2) - objh(5).Children.YData(1))/2;
    for j = 5:8
        objh(j).Children.XData(1) = objh(j).Children.XData(3) - BarWidth;
        objh(j).Children.XData(2) = objh(j).Children.XData(3) - BarWidth;
        objh(j).Children.XData(5) = objh(j).Children.XData(3) - BarWidth;
    end

ax.LineWidth = LineWidth;
ax.FontUnits = 'inches';
ax.FontSize = FontSize;
ax.Units = 'inches';
ax.Position = [axLeft axBottom axWidth axHeight];
xlabel('Rank')
ax.XTickLabel = {'RF' 'RerF' 'RR-RF' 'XGBoost'};
ax.XLim = [0.5 4.5];
ax.YLim = [0 60];

ColoredIdx = [1,3];
for j = ColoredIdx
    p = patch([j-0.5 j+0.5 j+0.5 j-0.5],[0 0 ax.YLim(2) ax.YLim(2)],...
        [0.9 0.9 0.9]);
    p.EdgeColor = 'none';
end

ColoredIdx = [2,4];
for j = ColoredIdx
    p = patch([j-0.5 j+0.5 j+0.5 j-0.5],[0 0 ax.YLim(2) ax.YLim(2)],...
        [0.8 0.8 0.8]);
    p.EdgeColor = 'none';
end

ch = ax.Children;
ch(1:4) = [];
ch(end+1:end+4) = ax.Children(1:4);
ax.Children = ch;

save_fig(gcf,'~/RandomerForest/Figures/Benchmark_ranks_2017_01_17',{'fig','pdf','png'})