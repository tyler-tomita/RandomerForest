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
FontSize = 0.175;
axWidth = 1.5;
axHeight = 1.5;
axLeft = [FontSize*4,FontSize*7+axWidth];
axBottom = FontSize*4*ones(1,2);
legWidth = axWidth*0.75;
legHeight = axHeight;
legLeft = axLeft(end) + axWidth;
legBottom = axBottom(end);
figWidth = legLeft + legWidth;
figHeight = axBottom(1) + axHeight + FontSize*2.5;

fig = figure;
fig.Units = 'inches';
fig.PaperUnits = 'inches';
fig.Position = [0 0 figWidth figHeight];
fig.PaperPosition = [0 0 figWidth figHeight];
fig.PaperSize = [figWidth figHeight];

Classifiers = {'rf','rerf','rr_rf'};

inPath1 = [rerfPath 'RandomerForest/Results/2017.04.01/Benchmarks/Raw/'];
inPath2 = '~/Benchmarks/Results/R/dat/Raw/';
contents = dir([inPath1 '*.mat']);

AbsoluteError = NaN(length(contents),length(Classifiers));
RelativeError = NaN(length(contents),length(Classifiers)-1);
ChanceProb = NaN(length(contents),1);
MeanDepth = NaN(length(contents),length(Classifiers));
Strength = NaN(length(contents),length(Classifiers));
Diversity = NaN(length(contents),length(Classifiers));
RelativeStrength = NaN(length(contents),length(Classifiers)-1);
RelativeDiversity = NaN(length(contents),length(Classifiers)-1);

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
                AbsoluteError(k,c) = TestError.(cl)(BestIdx.(cl));
                MeanDepth(k,c) = mean(Depth.(cl)(:,BestIdx.(cl)));
                Strength(k,c) = TreeStrength.(cl)(BestIdx.(cl));
                Diversity(k,c) = TreeDiversity.(cl)(BestIdx.(cl));
            else
                AbsoluteError(k,c) = dlmread([inPath2 Dataset '_testError.dat']);
                MeanDepth(k,c) = dlmread([inPath2 Dataset '_depth.dat']);
            end
            if c > 1
                RelativeError(k,c-1) = (AbsoluteError(k,c)-AbsoluteError(k,1))/ChanceProb(k);
                RelativeStrength(k,c-1) = Strength(k,c)-Strength(k,1);
                RelativeDiversity(k,c-1) = Diversity(k,c)-Diversity(k,1);
            end
        end
        k = k + 1;
    end
end

RelativeError(all(isnan(RelativeError),2),:) = [];
MeanDepth(all(isnan(MeanDepth),2),:) = [];
Strength(all(isnan(Strength),2),:) = [];
Diversity(all(isnan(Diversity),2),:) = [];
RelativeStrength(all(isnan(RelativeStrength),2),:) = [];
RelativeDiversity(all(isnan(RelativeDiversity),2),:) = [];

k = 1;
ax(k) = axes;
hold on
for c = 2:length(Classifiers)
    cl = Classifiers{c};
    plot(RelativeStrength(:,c-1),RelativeError(:,c-1),'.','MarkerSize',MarkerSize,...
        'Color',Colors.(cl))
end

xlabel('Relative Strength')
ylabel('Relative Error')

ax(k).LineWidth = LineWidth;
ax(k).FontUnits = 'inches';
ax(k).FontSize = FontSize;
ax(k).Units = 'inches';
ax(k).Position = [axLeft(k) axBottom(k) axWidth axHeight];
ax(k).XLim = [min(RelativeStrength(:)) max(RelativeStrength(:))+0.02];
ax(k).YLim = [min(RelativeError(:))-0.04 max(RelativeError(:))+0.02];
% ax(k).XScale = 'log';
% ax(k).YScale = 'log';

k = 2;
ax(k) = axes;
hold on
for c = 2:length(Classifiers)
    cl = Classifiers{c};
    plot(RelativeDiversity(:,c-1),RelativeError(:,c-1),'.','MarkerSize',MarkerSize,...
        'Color',Colors.(cl))
end

xlabel('Relative Diversity')

ax(k).LineWidth = LineWidth;
ax(k).FontUnits = 'inches';
ax(k).FontSize = FontSize;
ax(k).Units = 'inches';
ax(k).Position = [axLeft(k) axBottom(k) axWidth axHeight];
ax(k).XLim = [min(RelativeDiversity(:)) max(RelativeDiversity(:))];
ax(k).YLim = [min(RelativeError(:)) max(RelativeError(:))];
ax(k).XTick([2,4]) = [];
% ax(k).XScale = 'log';
% ax(k).YScale = 'log';

% k =3;
% ax(k) = axes;
% hold on
% for c = 2:length(Classifiers)
%     cl = Classifiers{c};
%     plot(RelativeDiversity(:,c-1)./RelativeStrength(:,c-1),RelativeError(:,c-1),'.','MarkerSize',MarkerSize,...
%         'Color',Colors.(cl))
% end
% 
% xlabel('Diversity/Strength')
% ylabel('Error Rate')
% 
% ax(k).LineWidth = LineWidth;
% ax(k).FontUnits = 'inches';
% ax(k).FontSize = FontSize;
% ax(k).Units = 'inches';
% ax(k).Position = [axLeft(k) axBottom(k) axWidth axHeight];
% % ax(k).XScale = 'log';
% % ax(k).YScale = 'log';
% 
lh = legend('RerF','RR-RF');
lh.Box = 'off';
lh.Units = 'inches';
lh.Position = [legLeft legBottom legWidth legHeight];

k = 3;
ax(k) = axes;
ax(k).FontUnits = 'inches';
ax(k).FontSize = FontSize;
ax(k).Units = 'inches';
ax(k).Position = [mean(axLeft) axBottom(1) axWidth axHeight];
ht = title('Raw Benchmarks');
ax(k).Visible = 'off';
ht.Visible = 'on';
        
save_fig(gcf,[rerfPath 'RandomerForest/Figures/pami/PAMI_fig10_strength_diversity_benchmark'],{'fig','pdf','png'})