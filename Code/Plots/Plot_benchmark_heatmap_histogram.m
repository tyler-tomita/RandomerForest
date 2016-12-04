%% Plot benchmark classifier rank distributions

clear
close all
clc

load('purple2green')
ColorMap = interpolate_colormap(ColorMap(round(size(ColorMap,1)/2):end,:),64,false);
% ColorMap = flipud(ColorMap);
% Colormap = parula;
FontSize = 12;
LineWidth = 2;

% LineWidth = 2;
% LineWidth_box = 4;
% LineWidth_whisker = 1.5;
% MarkerSize = 6;
% FontSize = .18;
% axWidth = 1.5;
% axHeight = 1.4;
% axLeft = FontSize*5*ones(1,5);
% axBottom = [FontSize*12+axHeight*4,FontSize*7+axHeight*3,...
%     FontSize*5+axHeight*2,FontSize*3+axHeight,...
%     FontSize];
% figWidth = axLeft(end) + axWidth + FontSize;
% figHeight = axBottom(1) + axHeight + FontSize*2;

% fig = figure;
% fig.Units = 'inches';
% fig.PaperUnits = 'inches';
% fig.Position = [0 0 figWidth figHeight];
% fig.PaperPosition = [0 0 figWidth figHeight];
% fig.PaperSize = [figWidth figHeight];

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

Transformations = {'Untransformed','Rotated','Scaled','Affine','Outlier'};

BinEdges = [-1,-0.2,-0.1,-0.05:0.005:0,0:0.005:0.05,0.1,0.2,1];

for i = 1
    load(['~/Benchmarks/Results/Benchmark_' lower(Transformations{i}) '.mat'])
    Classifiers = fieldnames(TestError{1});
    Classifiers(~ismember(Classifiers,{'rf','rerf','rerfr','frc','frcr','rr_rf','rr_rfr'})) = [];

    TestError = TestError(~cellfun(@isempty,TestError));
    
    AbsoluteError = NaN(length(TestError),length(Classifiers));
    RelativeError = NaN(length(TestError),length(Classifiers)-1);

    for j = 1:length(TestError)
        for k = 1:length(Classifiers)
            AbsoluteError(j,k) = TestError{j}.(Classifiers{k});
            if k > 1
                RelativeError(j,k-1) = TestError{j}.(Classifiers{k}) - TestError{j}.(Classifiers{1});
            end
        end
    end
    
    Counts = zeros(length(BinEdges)-1,length(Classifiers)-1);
    
    for k = 1:length(Classifiers)-1
        h = histogram(RelativeError(:,k),BinEdges);
        Counts(:,k) = h.Values';
    end
        Counts(length(BinEdges)/2,:) = sum(RelativeError==0);
        Counts(length(BinEdges)/2+1,:) = Counts(length(BinEdges)/2+1,:) - Counts(length(BinEdges)/2,:);
end

h = heatmap(flipud(Counts),{'RerF','RerF(r)','F-RC','Frank','RR-RF','RR-RF(r)'},cellstr(num2str(flipud(BinEdges'))),ColorMap,true);
ylabel('Relative Error')
colorbar;
h.FontSize = FontSize;
hold on
for i = 2:length(Classifiers)-1
    plot([i-0.5,i-0.5],h.YLim,'-k','LineWidth',LineWidth+1)
end

save_fig(gcf,'~/RandomerForest/Figures/Benchmark_heatmap_histogram')