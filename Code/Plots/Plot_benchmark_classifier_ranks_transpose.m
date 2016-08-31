% %% Plot Performance Profiles for Benchmark Transformations
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

FontSize = 16;

figWidth = 11;
figHeight = 8.5;

fig = figure;
fig.Units = 'inches';
fig.PaperUnits = 'inches';
fig.Position = [0 0 figWidth figHeight];
fig.PaperPosition = [0 0 figWidth figHeight];
fig.PaperSize = [figWidth figHeight];

runSims = false;

load('~/Benchmarks/Results/Benchmark_untransformed_2016_08_05.mat')

Classifiers = fieldnames(TestError{1});
Classifiers(ismember(Classifiers,{'rerf' 'rerfr' 'rerfd'})) = [];

TestError = TestError(~cellfun(@isempty,TestError));

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
    Bars(i).EdgeColor = 'w';
    Bars(i).BarWidth = 1;
end

xlabel('Rank')
ylabel('Frequency')
title('Untransformed')

ax.FontSize = FontSize;
ax.XTickLabel = {'RF' 'RF(r)' 'F-RC' 'F-RC(r)' 'RR-RF' 'RR-RF(r)'};
ax.YLim = [0 60];
ax.LineWidth = 2;

l = legend('1st place','2nd place','3rd place','4th place','5th place',...
    '6th place');
l.Location = 'northwest';
l.Box = 'off';

line([1 1],[30 40],'Color','k','LineWidth',3)
line([3 3],[29 40],'Color','k','LineWidth',3)
line([1 3],[40 40],'Color','k','LineWidth',3)
t = text(2,40,'*','HorizontalAlignment','center',...
    'VerticalAlignment','bottom','FontSize',FontSize+2);

line([1 1],[41 45],'Color','k','LineWidth',3)
line([4 4],[29 45],'Color','k','LineWidth',3)
line([1 4],[45 45],'Color','k','LineWidth',3)
t = text(2.5,45,'*','HorizontalAlignment','center',...
    'VerticalAlignment','bottom','FontSize',FontSize+2);

% line([1 1],[21.5 23.5],'Color','k','LineWidth',3)
% line([7 7],[22.5 23.5],'Color','k','LineWidth',3)
% line([1 7],[23.5 23.5],'Color','k','LineWidth',3)
% t = text(4,23.5,'**','HorizontalAlignment','center',...
%     'VerticalAlignment','bottom','FontSize',FontSize+2);

save_fig(gcf,'~/Benchmarks/Figures/Benchmark_untransformed_ranks')

%% Plot Performance Profiles for Rotated Benchmarks
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

FontSize = 16;

figWidth = 11;
figHeight = 8.5;

fig = figure;
fig.Units = 'inches';
fig.PaperUnits = 'inches';
fig.Position = [0 0 figWidth figHeight];
fig.PaperPosition = [0 0 figWidth figHeight];
fig.PaperSize = [figWidth figHeight];

runSims = false;

load('~/Benchmarks/Results/Benchmark_rotated_2016_08_05.mat')

Classifiers = fieldnames(TestError{1});
Classifiers(ismember(Classifiers,{'rerf' 'rerfr' 'rerfd'})) = [];

TestError = TestError(~cellfun(@isempty,TestError));

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
    Bars(i).EdgeColor = 'w';
    Bars(i).BarWidth = 1;
end

xlabel('Rank')
ylabel('Frequency')
title('Rotated')

ax.FontSize = FontSize;
ax.XTickLabel = {'RF' 'RF(r)' 'F-RC' 'F-RC(r)' 'RR-RF' 'RR-RF(r)'};
ax.YLim = [0 60];
ax.LineWidth = 2;

l = legend('1st place','2nd place','3rd place','4th place','5th place',...
    '6th place');
l.Location = 'northwest';
l.Box = 'off';

line([1 1],[29 36],'Color','k','LineWidth',3)
line([3 3],[35 36],'Color','k','LineWidth',3)
line([1 3],[36 36],'Color','k','LineWidth',3)
t = text(2,36,'*','HorizontalAlignment','center',...
    'VerticalAlignment','bottom','FontSize',FontSize+2);

line([1 1],[37 42],'Color','k','LineWidth',3)
line([4 4],[41 42],'Color','k','LineWidth',3)
line([1 4],[42 42],'Color','k','LineWidth',3)
t = text(2.5,42,'*','HorizontalAlignment','center',...
    'VerticalAlignment','bottom','FontSize',FontSize+2);

line([1 1],[43 45],'Color','k','LineWidth',3)
line([5 5],[43 45],'Color','k','LineWidth',3)
line([1 5],[45 45],'Color','k','LineWidth',3)
t = text(3,45,'*','HorizontalAlignment','center',...
    'VerticalAlignment','bottom','FontSize',FontSize+2);

line([1 1],[46 48],'Color','k','LineWidth',3)
line([6 6],[27 48],'Color','k','LineWidth',3)
line([1 6],[48 48],'Color','k','LineWidth',3)
t = text(3.5,48,'*','HorizontalAlignment','center',...
    'VerticalAlignment','bottom','FontSize',FontSize+2);

save_fig(gcf,'~/Benchmarks/Figures/Benchmark_rotated_ranks')

%% Plot Performance Profiles for Scaled Benchmarks
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

FontSize = 16;

figWidth = 11;
figHeight = 8.5;

fig = figure;
fig.Units = 'inches';
fig.PaperUnits = 'inches';
fig.Position = [0 0 figWidth figHeight];
fig.PaperPosition = [0 0 figWidth figHeight];
fig.PaperSize = [figWidth figHeight];

Transformations = {'Untransformed' 'Rotated' 'Scaled' 'Affine'};

runSims = false;

load('~/Benchmarks/Results/Benchmark_scaled_2016_08_05.mat')

Classifiers = fieldnames(TestError{1});
Classifiers(ismember(Classifiers,{'rerf' 'rerfr' 'rerfd'})) = [];

TestError = TestError(~cellfun(@isempty,TestError));

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
    Bars(i).EdgeColor = 'w';
    Bars(i).BarWidth = 1;
end

xlabel('Rank')
ylabel('Frequency')
title('Scaled')

ax.FontSize = FontSize;
ax.XTickLabel = {'RF' 'RF(r)' 'F-RC' 'F-RC(r)' 'RR-RF' 'RR-RF(r)'};
ax.YLim = [0 90];
ax.LineWidth = 2;

l = legend('1st place','2nd place','3rd place','4th place','5th place',...
    '6th place');
l.Location = 'northwest';
l.Box = 'off';

save_fig(gcf,'~/Benchmarks/Figures/Benchmark_scaled_ranks')

%% Plot Performance Profiles for Affine Benchmarks
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

FontSize = 16;

figWidth = 11;
figHeight = 8.5;

fig = figure;
fig.Units = 'inches';
fig.PaperUnits = 'inches';
fig.Position = [0 0 figWidth figHeight];
fig.PaperPosition = [0 0 figWidth figHeight];
fig.PaperSize = [figWidth figHeight];

Transformations = {'Untransformed' 'Rotated' 'Scaled' 'Affine'};

runSims = false;

load('~/Benchmarks/Results/Benchmark_affine_2016_08_05.mat')

Classifiers = fieldnames(TestError{1});
Classifiers(ismember(Classifiers,{'rerf' 'rerfr' 'rerfd'})) = [];

TestError = TestError(~cellfun(@isempty,TestError));

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
    Bars(i).EdgeColor = 'w';
    Bars(i).BarWidth = 1;
end

xlabel('Rank')
ylabel('Frequency')
title('Affine')

ax.FontSize = FontSize;
ax.XTickLabel = {'RF' 'RF(r)' 'F-RC' 'F-RC(r)' 'RR-RF' 'RR-RF(r)'};
ax.YLim = [0 95];
ax.LineWidth = 2;

l = legend('1st place','2nd place','3rd place','4th place','5th place',...
    '6th place');
l.Location = 'northwest';
l.Box = 'off';

line([1 1],[33 39],'Color','k','LineWidth',3)
line([4 4],[38 39],'Color','k','LineWidth',3)
line([1 4],[39 39],'Color','k','LineWidth',3)
t = text(2.5,39,'*','HorizontalAlignment','center',...
    'VerticalAlignment','bottom','FontSize',FontSize+2);

line([1 1],[40 43],'Color','k','LineWidth',3)
line([6 6],[38 43],'Color','k','LineWidth',3)
line([1 6],[43 43],'Color','k','LineWidth',3)
t = text(3.5,43,'*','HorizontalAlignment','center',...
    'VerticalAlignment','bottom','FontSize',FontSize+2);

save_fig(gcf,'~/Benchmarks/Figures/Benchmark_affine_ranks')

%% Plot Performance Profiles for Outlier Benchmarks
% close all
% clear
% clc
% 
% C = [0 1 1;0 1 0;1 0 1;1 0 0;0 0 0;1 .5 0];
% Colors.rf = C(1,:);
% Colors.rerf = C(2,:);
% Colors.rf_rot = C(3,:);
% Colors.rerfr = C(4,:);
% Colors.frc = C(5,:);
% Colors.rerfdn = C(6,:);
% 
% FontSize = 16;
% 
% % figWidth = 16;
% % figHeight = 7;
% figWidth = 11;
% figHeight = 8.5;
% 
% fig = figure;
% fig.Units = 'inches';
% fig.PaperUnits = 'inches';
% fig.Position = [0 0 figWidth figHeight];
% fig.PaperPosition = [0 0 figWidth figHeight];
% fig.PaperSize = [figWidth figHeight];
% 
% Transformations = {'Untransformed' 'Rotated' 'Scaled' 'Affine'};
% 
% runSims = false;
% 
% load('~/Benchmarks/Results/Benchmark_outlier_p_lte_7_2016_08_25.mat')
% 
% TestError = TestError(~cellfun(@isempty,TestError));
% Classifiers = fieldnames(TestError{1});
% Classifiers(~ismember(Classifiers,{'frc','frcr','frcn','frcz'})) = [];
% % Classifiers(~ismember(Classifiers,{'rerf','rerfr','rerfn','rerfz'})) = [];
% % Classifiers(~ismember(Classifiers,{'rf','rfr','rfn','rfz'})) = [];
% % Classifiers(~ismember(Classifiers,{'rr_rf','rr_rfr','rr_rfn','rr_rfz'})) = [];
% % Classifiers(~ismember(Classifiers,{'rerf','rerfr','rerfn','rerfz',...
% %     'frc','frcr','frcn','frcz',...
% %     'rr_rf','rr_rfr','rr_rfn','rr_rfz'})) = [];
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
% % subplot(1,2,1)
% bar(RankCounts')
% ax = gca;
% Bars = allchild(ax);
% for i = 1:length(Bars)
%     Bars(i).EdgeColor = 'w';
%     Bars(i).BarWidth = 1;
% end
% 
% % xlabel('Rank')
% % xlabel('RerF')
% xlabel('F-RC')
% % xlabel('RF')
% % xlabel('RR-RF')
% ylabel('Frequency')
% title('25 Benchmark Datasets w/ Outliers added')
% 
% ax.FontSize = FontSize;
% % ax.XTickLabel = {'F-RC' 'F-RC(r)' 'F-RC(n)' 'F-RC(z)'};
% % ax.XTickLabel = {'RerF' 'RerF(rank)' 'RerF(normalize)' 'RerF(z-score)'};
% ax.XTickLabel = {'raw' 'rank' 'normalize' 'z-score'};
% % ax.XTickLabel = {'RerF' 'RerF(r)' 'RerF(n)' 'RerF(z)','F-RC' 'F-RC(r)' 'F-RC(n)' 'F-RC(z)','RR-RF' 'RR-RF(r)' 'RR-RF(n)' 'RR-RF(z)'};
% ax.YLim = [0 20];
% ax.LineWidth = 2;
% 
% % l = legend('1st place','2nd place','3rd place','4th place');
% % l.Location = 'northwest';
% % l.Box = 'off';
% 
% % line([1 1],[14.5 19.5],'Color','k','LineWidth',3)
% % line([3 3],[18.5 19.5],'Color','k','LineWidth',3)
% % line([1 3],[19.5 19.5],'Color','k','LineWidth',3)
% % t = text(2,19.5,'*','HorizontalAlignment','center',...
% %     'VerticalAlignment','bottom','FontSize',FontSize+2);
% % 
% % line([1 1],[20 21],'Color','k','LineWidth',3)
% % line([5 5],[16.5 21],'Color','k','LineWidth',3)
% % line([1 5],[21 21],'Color','k','LineWidth',3)
% % t = text(3,21,'*','HorizontalAlignment','center',...
% %     'VerticalAlignment','bottom','FontSize',FontSize+2);
% % 
% % line([1 1],[21.5 23.5],'Color','k','LineWidth',3)
% % line([7 7],[22.5 23.5],'Color','k','LineWidth',3)
% % line([1 7],[23.5 23.5],'Color','k','LineWidth',3)
% % t = text(4,23.5,'**','HorizontalAlignment','center',...
% %     'VerticalAlignment','bottom','FontSize',FontSize+2);
% 
% % save_fig(gcf,'~/Benchmarks/Figures/Classifier_ranks_outlier_transpose')

%% Plot Performance Profiles for Affine Rescale Methods Benchmarks
% % close all
% % clear
% % clc
% % 
% % C = [0 1 1;0 1 0;1 0 1;1 0 0;0 0 0;1 .5 0];
% % Colors.rf = C(1,:);
% % Colors.rerf = C(2,:);
% % Colors.rf_rot = C(3,:);
% % Colors.rerfr = C(4,:);
% % Colors.frc = C(5,:);
% % Colors.rerfdn = C(6,:);
% % 
% % FontSize = 16;
% % 
% % figWidth = 11;
% % figHeight = 8.5;
% % 
% % fig = figure;
% % fig.Units = 'inches';
% % fig.PaperUnits = 'inches';
% % fig.Position = [0 0 figWidth figHeight];
% % fig.PaperPosition = [0 0 figWidth figHeight];
% % fig.PaperSize = [figWidth figHeight];
% % 
% % Transformations = {'Untransformed' 'Rotated' 'Scaled' 'Affine'};
% % 
% % runSims = false;
% % 
% load('~/Benchmarks/Results/Benchmark_affine_rescale_method_comparison_p_lte_7_2016_08_23.mat')
% 
% TestError = TestError(~cellfun(@isempty,TestError));
% Classifiers = fieldnames(TestError{1});
% % Classifiers(~ismember(Classifiers,{'frc','frcr','frcn','frcz'})) = [];
% Classifiers(~ismember(Classifiers,{'rerf','rerfr','rerfn','rerfz'})) = [];
% % Classifiers(~ismember(Classifiers,{'rerf','rerfr','rerfn','rerfz',...
% %     'frc','frcr','frcn','frcz',...
% %     'rr_rf','rr_rfr','rr_rfn','rr_rfz'})) = [];
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
% subplot(1,2,2)
% bar(RankCounts')
% ax = gca;
% Bars = allchild(ax);
% for i = 1:length(Bars)
%     Bars(i).EdgeColor = 'w';
%     Bars(i).BarWidth = 1;
% end
% 
% % xlabel('Rank')
% xlabel('RerF')
% ylabel('Frequency')
% title('25 Benchmark Datasets Affine-Transformed')
% 
% ax.FontSize = FontSize;
% % ax.XTickLabel = {'F-RC' 'F-RC(r)' 'F-RC(n)' 'F-RC(z)'};
% % ax.XTickLabel = {'RerF' 'RerF(rank)' 'RerF(normalize)' 'RerF(z-score)'};
% ax.XTickLabel = {'raw' 'rank' 'normalize' 'z-score'};
% % ax.XTickLabel = {'RerF' 'RerF(r)' 'RerF(n)' 'RerF(z)','F-RC' 'F-RC(r)' 'F-RC(n)' 'F-RC(z)','RR-RF' 'RR-RF(r)' 'RR-RF(n)' 'RR-RF(z)'};
% ax.YLim = [0 12];
% ax.LineWidth = 2;
% 
% % l = legend('1st place','2nd place','3rd place','4th place');
% % l.Location = 'northwest';
% % l.Box = 'off';
% % 
% % % line([1 1],[14.5 19.5],'Color','k','LineWidth',3)
% % % line([3 3],[18.5 19.5],'Color','k','LineWidth',3)
% % % line([1 3],[19.5 19.5],'Color','k','LineWidth',3)
% % % t = text(2,19.5,'*','HorizontalAlignment','center',...
% % %     'VerticalAlignment','bottom','FontSize',FontSize+2);
% % % 
% % % line([1 1],[20 21],'Color','k','LineWidth',3)
% % % line([5 5],[16.5 21],'Color','k','LineWidth',3)
% % % line([1 5],[21 21],'Color','k','LineWidth',3)
% % % t = text(3,21,'*','HorizontalAlignment','center',...
% % %     'VerticalAlignment','bottom','FontSize',FontSize+2);
% % % 
% % % line([1 1],[21.5 23.5],'Color','k','LineWidth',3)
% % % line([7 7],[22.5 23.5],'Color','k','LineWidth',3)
% % % line([1 7],[23.5 23.5],'Color','k','LineWidth',3)
% % % t = text(4,23.5,'**','HorizontalAlignment','center',...
% % %     'VerticalAlignment','bottom','FontSize',FontSize+2);
% % 
% % save_fig(gcf,'~/Benchmarks/Figures/Classifier_ranks_affine_rescale_methods_transpose')
% save_fig(gcf,'~/Benchmarks/Figures/Classifier_ranks_outliers_and_affine_rescale_methods')