%% Plot benchmark classifier rank distributions

clear
close all
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

load('purple2green')
ColorMap = interpolate_colormap(ColorMap(round(size(ColorMap,1)/2):end,:),64,false);

LineWidth = 2;
FontSize = .2;
axWidth = 2;
axHeight = 2;
cbWidth = axWidth;
cbHeight = 0.15;
axBottom = FontSize*5*ones(1,3);
axLeft = fliplr([FontSize*7+axHeight*2,FontSize*6+axHeight,...
    FontSize*5]);
cbLeft = axLeft;
cbBottom = FontSize*2*ones(1,3);
% figWidth = cbLeft(1) + cbWidth + FontSize*2;
figWidth = axLeft(end) + axWidth + FontSize;
figHeight = axBottom(1) + axHeight + FontSize*2;

fig = figure;
fig.Units = 'inches';
fig.PaperUnits = 'inches';
fig.Position = [0 0 figWidth figHeight];
fig.PaperPosition = [0 0 figWidth figHeight];
fig.PaperSize = [figWidth figHeight];

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

Transformations = {'Untransformed','Affine','Outlier'};

BinEdges = [-1,-0.2,-0.1,-0.05:0.01:-0.01,-0.005,0,0,0.005,0.01:0.01:0.05,0.1,0.2,1];

load Benchmark_plus_20_percent_outliers_data
% for i = 1
for i = 1:length(Transformations)
    load(['~/Benchmarks/Results/Benchmark_' lower(Transformations{i}) '.mat'])
    Classifiers = fieldnames(TestError{1});
    Classifiers(~ismember(Classifiers,{'rf','frc','frcr','rr_rf','rr_rfr'})) = [];

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
    Fractions = Counts./repmat(sum(Counts),size(Counts,1),1);
        
    ax(i) = axes;
    if i == 1
        YTLabel = {'F-RC','Frank','RR-RF','RR-RF(r)'};
    else
        YTLabel = {''};
    end
    h = heatmap(Fractions',cellstr(num2str(BinEdges')),YTLabel,...
        ColorMap,true,'horizontal');
    if i==1
        xlabel({'Normalized Error';'Relative to RF'})
        title('Raw')
    elseif i==3
        title('Corrupted')
    else
        title(Transformations{i})
    end
    
    if i == 2
        cb = colorbar;
        cb.Location = 'southoutside';
        xlh = xlabel(cb,'Fraction of Datasets');
        cb.Ticks = [];
        cb.Units = 'inches';
        cb.Position = [cbLeft(i) cbBottom(i) cbWidth cbHeight];
        cb.Box = 'off';
        cb.FontSize = 16;
        xlh.Position = [0.2237 -0.5 0];
    end
%     h.FontSize = FontSize;
    hold on
    for k = 2:length(Classifiers)-1
%         plot([k-0.5,k-0.5],h.YLim,'-k','LineWidth',LineWidth+1)
        plot(h.XLim,[k-0.5,k-0.5],'-k','LineWidth',LineWidth)
    end
    ax(i).XTick = [0.5,10,19.5];
    ax(i).XTickLabel = {'-1';'0';'1'};
    ax(i).XTickLabelRotation = 0;
    ax(i).TickLength = [0 0];
    ax(i).LineWidth = LineWidth;
    ax(i).FontUnits = 'inches';
    ax(i).FontSize = FontSize;
    ax(i).Units = 'inches';
    ax(i).Position = [axLeft(i) axBottom(i) axWidth axHeight];
end

save_fig(gcf,[rerfPath 'RandomerForest/Figures/ROFLMAO_fig5_benchmark_heatmap_slide_version'],{'fig','pdf','png'})