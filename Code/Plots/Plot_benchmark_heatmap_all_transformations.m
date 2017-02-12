%% Plot benchmark classifier rank distributions

clear
close all
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

load('purple2green')
ColorMap = interpolate_colormap(ColorMap(round(size(ColorMap,1)/2):end,:),64,false);

scaleFactor = .75;
LineWidth = 2*scaleFactor;
FontSize = .115*scaleFactor;
titleFontSize = 10*scaleFactor;
axWidth = 2*scaleFactor;
axHeight = 2*scaleFactor;
cbWidth = .15*scaleFactor;
cbHeight = axHeight;
axLeft = FontSize*5.5*ones(1,5);
axBottom = [FontSize*23+axHeight*4,FontSize*18+axHeight*3,...
    FontSize*13+axHeight*2,FontSize*8+axHeight,...
    FontSize*3];
axLeft2 = fliplr(axBottom);
axBottom2 = axLeft;
axLeft = axLeft2;
axBottom = axBottom2;
cbLeft = axLeft + axWidth + FontSize/2;
cbBottom = axBottom;
figWidth = cbLeft(end) + cbWidth + FontSize*2;
figHeight = axBottom(1) + axHeight + FontSize*2;

% % FontSize = .2;
% axWidth = 1.5;
% axHeight = 1.75;
% cbWidth = .1;
% cbHeight = axHeight;
% axLeft = repmat([FontSize*5.5,FontSize*12+axWidth],1,3);
% axLeft(end-1) = mean(axLeft(end-1:end));
% axLeft(end) = [];
% axBottom = [(FontSize*8+axHeight*2)*ones(1,2),...
%     (FontSize*5+axHeight)*ones(1,2),FontSize*2*ones(1,2)];
% cbLeft = axLeft + axWidth + FontSize/2;
% cbBottom = axBottom;
% % legWidth = axWidth;
% % legHeight = axHeight;
% % legLeft = axLeft(end) + axWidth*2/3 + FontSize;
% % legBottom = axBottom(end);
% figWidth = cbLeft(2) + cbWidth + FontSize*1.5;
% figHeight = axBottom(1) + axHeight + FontSize*2;

fig = figure;
fig.Units = 'inches';
fig.PaperUnits = 'inches';
fig.Position = [0 0 figWidth figHeight];
fig.PaperPosition = [0 0 figWidth figHeight];
fig.PaperSize = [figWidth figHeight];

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

Transformations = {'Untransformed','Rotated','Scaled','Affine','Outlier'};

BinEdges = [-1,-0.2,-0.1,-0.05:0.01:-0.01,-0.005,0,0,0.005,0.01:0.01:0.05,0.1,0.2,1];

load Benchmark_plus_20_percent_outliers_data

for i = 1:length(Transformations)
    load(['~/Benchmarks/Results/Benchmark_' lower(Transformations{i}) '.mat'])
    Classifiers = fieldnames(TestError{1});
    Classifiers(~ismember(Classifiers,{'rf','rerf','rerfr','frc','frcr','rr_rf','rr_rfr'})) = [];

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
        
    ax(i) = axes;
    h = heatmap(flipud(Counts),{'RerF','RerF(r)','F-RC','Frank','RR-RF','RR-RF(r)'},cellstr(num2str(flipud(BinEdges'))),ColorMap,true);
    if i==1
        ylabel('Normalized Error Relative to RF')
        title('Raw','FontSize',titleFontSize)
    elseif i==5
        title('Corrupted','FontSize',titleFontSize)
    else
        title(Transformations{i},'FontSize',titleFontSize)
    end
    cb = colorbar;
    cb.Units = 'inches';
    cb.Position = [cbLeft(i) cbBottom(i) cbWidth cbHeight];
    cb.Box = 'off';
    cb.FontSize = 8;
    h.FontSize = FontSize;
    hold on
    for k = 2:length(Classifiers)-1
        plot([k-0.5,k-0.5],h.YLim,'-k','LineWidth',LineWidth+1)
    end
    ax(i).LineWidth = LineWidth;
    ax(i).FontUnits = 'inches';
    ax(i).FontSize = FontSize;
    ax(i).Units = 'inches';
    ax(i).Position = [axLeft(i) axBottom(i) axWidth axHeight];
    if i>1
        ax(i).YTickLabel = {};
    end
end

save_fig(gcf,[rerfPath 'RandomerForest/Figures/Benchmark_heatmap_all_transformations_2017_02_01'],{'fig','pdf','png'})