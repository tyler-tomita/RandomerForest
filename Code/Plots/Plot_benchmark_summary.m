%% Plot misclassification rate, strength, bias, and variance for each benchmark

% Iterate through each benchmark summary file, plot each metric against
% mtry for all classifiers

close all
clear
clc

C = [0 1 1;0 1 0;1 0 1;1 0 0;0 0 0;1 .5 0];
Colors.rf = C(1,:);
Colors.rerf = C(2,:);
Colors.rf_rot = C(3,:);
Colors.rerfr = C(4,:);
Colors.frc = C(5,:);
Colors.rerfdn = C(6,:);

LineWidth.rf = 2;
LineWidth.rerf = 2;
LineWidth.rf_rot = 2;
LineWidth.rerfr = 2;
LineWidth.frc = 2;

FontSize = .16;
axWidth = 1.5;
axHeight = 1.5;
axLeft = [FontSize*5,FontSize*10+axWidth,FontSize*15+axWidth*2,FontSize*20+axWidth*3];
axBottom = [FontSize*4,FontSize*4,FontSize*4,FontSize*4];
figWidth = axLeft(end) + axWidth + FontSize*4;
figHeight = axBottom(1) + axHeight + FontSize*4;

fig = figure;
fig.Units = 'inches';
fig.PaperUnits = 'inches';
fig.Position = [0 0 figWidth figHeight];
fig.PaperPosition = [0 0 figWidth figHeight];
fig.PaperSize = [figWidth figHeight];

InPath = '../Results/Summary/Untransformed/';
OutPath = '../Figures/Untransformed/';

load('../Data/Matlab/Benchmark_data.mat')

Metrics = {'MR','S','V','B'};
MetricNames = {'Forest Error','Tree Error','Tree Variance',...
    'Tree Bias'};

ClassifierNames = containers.Map({'rf','rerf','rf_rot','rerfr','frc'},...
    {'RF','RerF','RotRF','RerF(rank)','F-RC'});

HasSummary = [];
ii = 1;
for i = 1:length(Datasets)
    InFile = [InPath,Datasets(i).Name,'_untransformed_summary.mat'];
    if exist(InFile)
        fprintf('Dataset %d: %s\n',i,Datasets(i).Name)
        load(InFile)
        HasSummary = [HasSummary i];
        fig = figure;
        fig.Units = 'inches';
        fig.PaperUnits = 'inches';
        fig.Position = [0 0 figWidth figHeight];
        fig.PaperPosition = [0 0 figWidth figHeight];
        fig.PaperSize = [figWidth figHeight];    
        for m = 1:length(Metrics)
            Metric = Metrics{m};
            Classifiers = fieldnames(Summary.(Metric));
            ax = subplot(1,length(Metrics),m);
            hold on
            LineNames = [];
            min_x = [];
            min_y = [];
            max_x = [];
            max_y = [];
            for c = 1:length(Classifiers)
                cl = Classifiers{c};
                if length(Summary.NMIX.(cl)) <= 0
                    plot(Summary.MTRY.(cl),Summary.(Metric).(cl),...
                        'Color',Colors.(cl),'LineWidth',LineWidth.(cl))
                    min_x = [min_x min(Summary.MTRY.(cl))];
                    max_x = [max_x max(Summary.MTRY.(cl))];
                    min_y = [min_y min(Summary.(Metric).(cl))];
                    max_y = [max_y max(Summary.(Metric).(cl))];
                else
                    plot(Summary.MTRY.(cl),Summary.(Metric).(cl)(:,1),...
                        'Color',Colors.(cl),'LineWidth',LineWidth.(cl))
                    min_x = [min_x min(Summary.MTRY.(cl))];
                    max_x = [max_x max(Summary.MTRY.(cl))];
                    min_y = [min_y min(Summary.(Metric).(cl)(:,1))];
                    max_y = [max_y max(Summary.(Metric).(cl)(:,1))];
                end
                LineNames = [LineNames,{ClassifierNames(cl)}];
                min_x = min(min_x);
                max_x = max(max_x);
                min_y = min(min_y);
                max_y = max(max_y);
            end
            if m == 1
                title(sprintf('%s (n = %d, p = %d)',Datasets(i).Name,...
                    Datasets(i).n,Datasets(i).p))
            end
            xlabel('mtry')
            ylabel(MetricNames(m))
            ax.LineWidth = 2;
            ax.FontUnits = 'inches';
            ax.FontSize = FontSize;
            ax.Units = 'inches';
            ax.Position = [axLeft(m),axBottom(m),axWidth,axHeight];
            ax.Box = 'off';
            ax.XScale = 'log';
        ax.XLim = [min_x,max_x];
            if min_y==0 && max_y==0
            ax.YLim = [0 1];
            else
            ax.YLim = [min_y,max_y];
            end
        end
        save_fig(gcf,[OutPath,Datasets(i).Name,'_summary'])
        close

        Classifiers = fieldnames(Summary.MR);
        for c = 1:length(Classifiers)
            cl = Classifiers{c};
            if isempty(Summary.NMIX.(cl))
                [MinMR.(cl)(ii),MinIdx] = min(Summary.MR.(cl));
                BestMTRY.(cl)(ii) = Summary.MTRY.(cl)(MinIdx);
            else
                [MinMR.(cl)(ii),MinIdx] = min(Summary.MR.(cl)(:,1));
                BestMTRY.(cl)(ii) = Summary.MTRY.(cl)(MinIdx);
            end
        end
        ii = ii + 1;
    end
end

%Plot distributions of differences in Lhats
MRDiff.rerf = MinMR.rerf - MinMR.rf;
MRDiff.frc = MinMR.frc - MinMR.rf;

p = plotSpread({MRDiff.rerf,MRDiff.frc},'distributionMarkers',{'.','.'},...
    'distributionColors',{'g','k'},'xNames',{'RerF','F-RC'});
ax = p{end};
ch = allchild(ax);
ch(1).MarkerSize = 10;
ch(2).MarkerSize = 10;
hold on
plot(ax.XLim,[0,0],'r--')
ylabel('Error (relative to RF)')
save_fig(gcf,[OutPath,'Error_difference_histogram'])

%Fraction of times mtry > p resulted in best performance binned by p
Win.rerf = MRDiff.rerf < 0;
Win.frc = MRDiff.frc < 0;
Subset = [Datasets(HasSummary).p];

binWidth = 20;
bins.rerf = floor(min(Subset(Win.rerf))/binWidth)*binWidth:binWidth:...
    ceil(max(Subset(Win.rerf))/binWidth)*binWidth;
bins.frc = floor(min(Subset(Win.frc))/binWidth)*binWidth:binWidth:...
    ceil(max(Subset(Win.frc))/binWidth)*binWidth;
bins = min(bins.rerf(1),bins.frc(1)):binWidth:max(bins.rerf(end),bins.frc(end));
for i = 1:length(bins)-1
    lb = bins(i);
    ub = bins(i+1);
    inBin = Subset > lb & Subset <= ub;
    Fraction.rerf(i) = sum(Win.rerf & BestMTRY.rerf > Subset & inBin)/...
        sum(Win.rerf & inBin);
    Fraction.frc(i) = sum(Win.frc & BestMTRY.frc > Subset & inBin)/...
        sum(Win.frc & inBin);
    xName{i} = sprintf('%d - %d',lb+1,ub);
end

figure;
b = bar([Fraction.rerf',Fraction.frc']);
b(1).FaceColor = 'g';
b(1).EdgeColor = 'g';
b(2).FaceColor = 'k';
b(2).EdgeColor = 'k';
ax = gca;
ax.Box = 'off';
ax.XTickLabel = xName;
xlabel('p')
ylabel('Fraction')
title('Fraction of cases where mtry > p is optimal')
l = legend('RerF','F-RC','Location','northwest');
l.Box = 'off';
save_fig(gcf,[OutPath,'Mtry_binned'])

% Unbinned fraction of times mtry > p resulted in best performance
Fraction.rerf = sum(Win.rerf & BestMTRY.rerf > Subset)/sum(Win.rerf);
Fraction.frc = sum(Win.frc & BestMTRY.frc > Subset)/sum(Win.frc);
figure;
ax = gca;
ax.Box = 'off';
hold on
YData = [Fraction.rerf,Fraction.frc];
bar(1,YData(1),'FaceColor','g','EdgeColor','g');
bar(2,YData(2),'FaceColor','k','EdgeColor','k');
ax.XTick = 1:length(YData);
ax.XTickLabel = {'RerF','F-RC'};
ylabel('Fraction')
title('Fraction of cases where mtry > p is optimal')
save_fig(gcf,[OutPath,'Mtry'])

pval.rerf = signrank(MinMR.rf,MinMR.rerf,'tail','right');
pval.frc = signrank(MinMR.rf,MinMR.frc,'tail','right');

save('../Results/pvalues_untransformed','pval')