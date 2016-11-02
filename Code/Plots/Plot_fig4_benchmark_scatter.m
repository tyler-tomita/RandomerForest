%% Plot benchmark classifier rank distributions

clear
close all
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);
LineWidth = 2;
FontSize = .12;
MarkerSize = 8;
axWidth = 1.75;
axHeight = 1.75;
axLeft = FontSize*4*ones(1,4);
axBottom = [FontSize*15+axHeight*3,...
    FontSize*11+axHeight*2,FontSize*7+axHeight,...
    FontSize*3];
figWidth = axLeft(end) + axWidth + FontSize*2;
figHeight = axBottom(1) + axHeight + FontSize*2;

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

Transformations = {'Untransformed','Rotated','Scaled','Affine','Outlier'};
Names = {'RF','RerF','RerFr','F-RC','Frank'};

for i = 1:length(Transformations)
    load(['~/Benchmarks/Results/Benchmark_' lower(Transformations{i}) '.mat'])
    Classifiers = fieldnames(TestError{1});
    Classifiers(~ismember(Classifiers,{'rf','rerf','rerfr','frc','frcr'})) = [];

    TestError = TestError(~cellfun(@isempty,TestError));
    
    ErrorMatrix = [];
    RelativeError = {};
    for j = 1:length(TestError)
        for k = 1:length(Classifiers)
            ErrorMatrix(j,k) = TestError{j}.(Classifiers{k});
            if k > 2
                RelativeError{k-1}(j,1) = TestError{j}.(Classifiers{k}) - TestError{j}.(Classifiers{1});
            end
        end
    end
    
    fig = figure;
    fig.Units = 'inches';
    fig.PaperUnits = 'inches';
    fig.Position = [0 0 figWidth figHeight];
    fig.PaperPosition = [0 0 figWidth figHeight];
    fig.PaperSize = [figWidth figHeight];
   
    for j = 2:length(Classifiers)
        ax(j-1) = axes;
        plot(ErrorMatrix(:,1),ErrorMatrix(:,j),'.','MarkerSize',MarkerSize);
        hold on
        plot([0 1],[0 1],'-k','LineWidth',1.5)
        xlabel(sprintf('Error Rate (%s)',Names{1}))
        ylabel(sprintf('Error Rate (%s)',Names{j}))
        ax(j-1).LineWidth = LineWidth;
        ax(j-1).FontUnits = 'inches';
        ax(j-1).FontSize = FontSize;
        ax(j-1).Units = 'inches';
        ax(j-1).Position = [axLeft(j-1) axBottom(j-1) axWidth axHeight];
        ax(j-1).Box = 'off';
        ax(j-1).XLim = [0 1];
        ax(j-1).YLim = ax(j-1).XLim;
        ax(j-1).XTick = ax(j-1).YTick;
        if j-1 == 1
            if i == 1
                title('Raw')
            else
                title(Transformations{i})
            end
        end
    end
    save_fig(gcf,[rerfPath sprintf('RandomerForest/Figures/Fig4_benchmark_%s_scatter',Transformations{i})])
end