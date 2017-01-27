%% Plot Performance Profiles for Corrupted Benchmarks
close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

load('purple2green')
ColorMap = interpolate_colormap(ColorMap(round(size(ColorMap,1)/2):end,:),64,false);

LineWidth = 2;
FontSize = .11;
axWidth = 2;
axHeight = 2;
cbWidth = .15;
cbHeight = axHeight;
axLeft = FontSize*5.5;
axBottom = FontSize*2;
cbLeft = axLeft + axWidth + FontSize/2;
cbBottom = axBottom;
figWidth = cbLeft(1) + cbWidth + FontSize*1.5;
figHeight = axBottom(1) + axHeight + FontSize*2;

fig = figure;
fig.Units = 'inches';
fig.PaperUnits = 'inches';
fig.Position = [0 0 figWidth figHeight];
fig.PaperPosition = [0 0 figWidth figHeight];
fig.PaperSize = [figWidth figHeight];

load ~/Benchmarks/Results/Benchmark_outlier
load Benchmark_plus_20_percent_outliers_data

BinEdges = [-1,-0.2,-0.1,-0.05:0.01:-0.01,-0.005,0,0,0.005,0.01:0.01:0.05,0.1,0.2,1];

Classifiers = fieldnames(TestError{1});
Classifiers(~ismember(Classifiers,{'frc','frcr','frcn','frcz'})) = [];

NotEmpty = find(~cellfun(@isempty,TestError));

ChanceProb = NaN(length(NotEmpty),1);
AbsoluteError = NaN(length(NotEmpty),length(Classifiers));
NormRelativeError = NaN(length(NotEmpty),length(Classifiers)-1);

for j = 1:length(NotEmpty)
    ClassCounts = histcounts(grp2idx(Datasets(NotEmpty(j)).Ytrain));
    ChanceProb(j) = 1 - max(ClassCounts)/sum(ClassCounts);
    for k = 1:length(Classifiers)
        AbsoluteError(j,k) = TestError{NotEmpty(j)}.(Classifiers{k});
    end
    for k = 1:length(Classifiers)
        if k > 1
            NormRelativeError(j,k-1) = (TestError{NotEmpty(j)}.(Classifiers{k})...
                - TestError{NotEmpty(j)}.(Classifiers{1}))/ChanceProb(j);
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

ax = axes;
h = heatmap(flipud(Counts),{'Frank','F-RC(n)','F-RC(z)'},cellstr(num2str(flipud(BinEdges'))),ColorMap,true);
ylabel('Normalized Error Relative to F-RC')
title('Corrupted')
cb = colorbar;
cb.Units = 'inches';
cb.Position = [cbLeft cbBottom cbWidth cbHeight];
cb.Box = 'off';
h.FontSize = FontSize;
hold on
for k = 2:length(Classifiers)-1
    plot([k-0.5,k-0.5],h.YLim,'-k','LineWidth',LineWidth+1)
end
ax.LineWidth = LineWidth;
ax.FontUnits = 'inches';
ax.FontSize = FontSize;
ax.Units = 'inches';
ax.Position = [axLeft axBottom axWidth axHeight];
 
% runSims = false;
% 
% load('~/Benchmarks/Results/Benchmark_outlier.mat')
% 
% TestError = TestError(~cellfun(@isempty,TestError));
% Classifiers = fieldnames(TestError{1});
% Classifiers(~ismember(Classifiers,{'frc','frcr','frcn','frcz'})) = [];
% 
% ErrorMatrix = [];
% for i = 1:length(TestError)
%     for j = 1:length(Classifiers)
%         ErrorMatrix(i,j) = TestError{i}.(Classifiers{j});
%     end
% end
% 
% ClRanks = tiedrank(ErrorMatrix')';
% IntRanks = floor(ClRanks);
% 
% RankCounts = NaN(length(Classifiers));
% for i = 1:length(Classifiers)
%     RankCounts(i,:) = sum(IntRanks==i);
% end
% 
% bar(RankCounts')
% ax = gca;
% Bars = allchild(ax);
% for i = 1:length(Bars)
%     Bars(i).EdgeColor = 'k';
%     Bars(i).BarWidth = 1;
% end
% 
% xlabel('Rank')
% ylabel('Frequency')
% l = legend('1st place','2nd place','3rd place','4th place');
% l.Box = 'off';
% l.FontSize = 9;
% l.Units = 'inches';
% l.Position = [legLeft legBottom legWidth legHeight];
% ax.LineWidth = LineWidth;
% ax.FontUnits = 'inches';
% ax.FontSize = FontSize;
% ax.Units = 'inches';
% ax.Position = [axLeft axBottom axWidth axHeight];
% ax.XTickLabel = {'F-RC' 'Frank' 'F-RC(n)' 'F-RC(z)'};
% ax.XLim = [0.5 4.5];
% ax.YLim = [0 60];
% 
% ColoredIdx = [1,3];
% for j = ColoredIdx
%     p = patch([j-0.5 j+0.5 j+0.5 j-0.5],[0 0 ax.YLim(2) ax.YLim(2)],...
%         [0.9 0.9 0.9]);
%     p.EdgeColor = 'none';
% end
% 
% ColoredIdx = [2,4];
% for j = ColoredIdx
%     p = patch([j-0.5 j+0.5 j+0.5 j-0.5],[0 0 ax.YLim(2) ax.YLim(2)],...
%         [0.8 0.8 0.8]);
%     p.EdgeColor = 'none';
% end
% 
% ch = ax.Children;
% ch(1:4) = [];
% ch(end+1:end+4) = ax.Children(1:4);
% ax.Children = ch;
% t = text(2,ax.YLim(2),'\bf{+}','HorizontalAlignment','center',...
%     'VerticalAlignment','top','FontSize',14,'Color','r');
% 
save_fig(gcf,[rerfPath 'RandomerForest/Figures/ROFLMAO_figS2_benchmark_scaling_methods_2017_01_23'],{'fig','pdf','png'})