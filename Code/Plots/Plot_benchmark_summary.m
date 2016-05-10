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
                        'Color',Colors.(cl),...
                'LineWidth',LineWidth.(cl))
                    min_x = [min_x min(Summary.MTRY.(cl)(:,1))];
                    max_x = [max_x max(Summary.MTRY.(cl)(:,1))];
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
                [MinMR.(cl)(i),MinIdx] = min(Summary.MR.(cl));
                BestMTRY.(cl)(i) = Summary.MTRY.(cl)(MinIdx);
            else
                [MinMR.(cl)(i),MinIdx] = min(Summary.MR.(cl)(:,1));
                BestMTRY.(cl)(i) = Summary.MTRY.(cl)(MinIdx,1);
            end
        end
    end
end

Win.rerf = MinMR.rerf < MinMR.rf;
Win.frc = MinMR.frc < MinMR.rf;
WinFraction.rerf = mean(Win.rerf);
WinFraction.frc = mean(Win.frc);
Subset = [Datasets(HasSummary).p];
BigMTRY.rerf = Win.rerf && BestMTRY.rerf > Subset;
BigMTRY.frc = Win.frc && BestMTRY.frc > Subset;

%Bar plot of fraction of times each of Rerf and FRC did better than RF
bar([WinFraction.rerf,WinFraction,frc],0.5,'FaceColor',[0 .5 .5],...
    'EdgeColor',[0 .5 .5])
ylabel('Fraction')
ax = gca;
ax.XTickLabel = {'RerF','F-RC'};
title('Fraction of times > RF')
save_fig(gcf,[OutPath,'Win_fraction'])

%Histogram of p when mtry > p resulted in best performance
histogram(Subset(BigMTRY.rerf),'FaceColor',[0 .5 .5],...
    'EdgeColor',[0 .5 .5])
xlabel('p')
ylabel('frequency')
title('Cases where mtry > p is optimal (RerF)')
save_fig(gcf,[OutPath,'Mtry_RerF'])

histogram(Subset(BigMTRY.frc),'FaceColor',[0 .5 .5],...
    'EdgeColor',[0 .5 .5])
xlabel('p')
ylabel('frequency')
title('Cases where mtry > p is optimal (F-RC)')
save_fig(gcf,[OutPath,'Mtry_FRC'])