%% Plot Performance Profiles for Corrupted Benchmarks
close all
clear
clc

LineWidth = 2;
FontSize = .12;
axWidth = 2.75;
axHeight = 1.25;
axLeft = FontSize*3;
axBottom = FontSize*3;
legWidth = 0.75*axHeight;
legHeight = 0.75*axHeight;
legLeft = axLeft(end) + axWidth + FontSize;
legBottom = axBottom - (legHeight - axHeight)/2;
figWidth = legLeft + legWidth + FontSize;
figHeight = axBottom(1) + axHeight + FontSize;

fig = figure;
fig.Units = 'inches';
fig.PaperUnits = 'inches';
fig.Position = [0 0 figWidth figHeight];
fig.PaperPosition = [0 0 figWidth figHeight];
fig.PaperSize = [figWidth figHeight];

runSims = false;

load('~/Benchmarks/Results/Benchmark_outlier.mat')

TestError = TestError(~cellfun(@isempty,TestError));
Classifiers = fieldnames(TestError{1});
Classifiers(~ismember(Classifiers,{'frc','frcr','frcn','frcz'})) = [];

ErrorMatrix = [];
for i = 1:length(TestError)
    for j = 1:length(Classifiers)
        ErrorMatrix(i,j) = TestError{i}.(Classifiers{j});
    end
end

ClRanks = tiedrank(ErrorMatrix')';
IntRanks = floor(ClRanks);

RankCounts = NaN(length(Classifiers));
for i = 1:length(Classifiers)
    RankCounts(i,:) = sum(IntRanks==i);
end

bar(RankCounts')
ax = gca;
Bars = allchild(ax);
for i = 1:length(Bars)
    Bars(i).EdgeColor = 'k';
    Bars(i).BarWidth = 1;
end

xlabel('Rank')
ylabel('Frequency')
l = legend('1st place','2nd place','3rd place','4th place');
l.Box = 'off';
l.FontSize = 9;
l.Units = 'inches';
l.Position = [legLeft legBottom legWidth legHeight];
ax.LineWidth = LineWidth;
ax.FontUnits = 'inches';
ax.FontSize = FontSize;
ax.Units = 'inches';
ax.Position = [axLeft axBottom axWidth axHeight];
ax.XTickLabel = {'F-RC' 'F-RC(r)' 'F-RC(n)' 'F-RC(z)'};
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
t = text(2,ax.YLim(2),'\bf{+}','HorizontalAlignment','center',...
    'VerticalAlignment','top','FontSize',14,'Color','r');

save_fig(gcf,'~/RandomerForest/Figures/FigS2_benchmark_scaling_methods')