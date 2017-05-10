%% Plot benchmark classifier rank distributions

clear
close all
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

load('purple2green')
ColorMap = interpolate_colormap(ColorMap(round(size(ColorMap,1)/2):end,:),64,false);
% ColorMap = flipud(interpolate_colormap(ColorMap(3:end-2,:),64,true));

LineWidth = 2;
FontSize = .2;
axWidth = 2;
axHeight = 2;
cbWidth = axWidth;
cbHeight = 0.15;
axBottom = FontSize*5*ones(1,5);
axLeft = fliplr([FontSize*9+axHeight*4,FontSize*8+axHeight*3,...
    FontSize*7+axHeight*2,FontSize*6+axHeight,...
    FontSize*5]);
cbLeft = axLeft;
cbBottom = FontSize*2*ones(1,5);
% figWidth = cbLeft(1) + cbWidth + FontSize*2;
figWidth = axLeft(end) + axWidth + FontSize;
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

Transformations = {'Raw','Rotated','Scaled','Affine','Corrupted'};
Classifiers = {'rf','rerf','rerfr','rr_rf','rr_rfr','xgb'};

BinEdges = [-1,-0.2,-0.1,-0.05:0.01:-0.01,-0.005,0,0,0.005,0.01:0.01:0.05,0.1,0.2,1];

DatasetNames = importdata('~/Benchmarks/Data/Names.txt');

for t = 1:length(Transformations)
    Classifiers = {'rf','rerf','rerfr','frc','frcr','rr_rf','rr_rfr','xgb'};
    
    if strcmp(Transformations{t},'Raw')
        S = load(['~/Benchmarks/Results/Benchmark_untransformed.mat']);
    elseif strcmp(Transformations{t},'Corrupted')
        S = load(['~/Benchmarks/Results/Benchmark_outlier.mat']);
    else
        S = load(['~/Benchmarks/Results/Benchmark_' lower(Transformations{t}) '.mat']);
    end
    
    fprintf('t = %d\n',t)
    
    inPath1 = [rerfPath 'RandomerForest/Results/2017.04.01/Benchmarks/' Transformations{t} '/'];
    inPath2 = ['~/Benchmarks/Results/R/dat/' Transformations{t} '/'];
    contents = dir([inPath1 '*.mat']);

    AbsoluteError = NaN(length(contents),length(Classifiers));
    NormRelativeError = NaN(length(contents),length(Classifiers)-1);
    ChanceProb = NaN(length(contents),1);
    
    k = 1;
    
    for i = 1:length(contents)
        Dataset = strsplit(contents(i).name,'.');
        Dataset = Dataset{1};
        DatasetIdx = find(strcmp(Dataset,DatasetNames));
        
        load([inPath1 contents(i).name])
        
        isComplete = true;
        
        for c = 1:length(Classifiers)
            cl = Classifiers{c};
            if ~strcmp(cl,'xgb') && ~strcmp(cl,'frc') && ~strcmp(cl,'frcr')
                if ~isfield(TestError,cl)
                    isComplete = false;
                end
            elseif strcmp(cl,'frc') || strcmp(cl,'frcr')
                if isempty(S.TestError{DatasetIdx})
                    isComplete = false;
                end
            else
                if ~exist([inPath2 Dataset '_testError.dat'])
                    isComplete = false;
                end
            end
        end
        
        if isComplete
            TrainSet = dlmread(['~/Benchmarks/Data/dat/' Transformations{t} '/' Dataset '_train.dat']);
            [ntrain,p] = size(TrainSet(:,1:end-1));
            nClasses = length(unique(TrainSet(:,end)));
            ClassCounts = histcounts(TrainSet(:,end),nClasses);
            ChanceProb(k) = 1 - max(ClassCounts)/sum(ClassCounts);
            
            for c = 1:length(Classifiers)
                cl = Classifiers{c};
                if ~strcmp(cl,'xgb') && ~strcmp(cl,'frc') && ~strcmp(cl,'frcr')
                    BI = hp_optimize(OOBError.(cl)(end,1:length(Params.(cl).d)),...
                        OOBAUC.(cl)(end,1:length(Params.(cl).d)));
                    BI = BI(end);
                    AbsoluteError(k,c) = TestError.(cl)(BI);
                elseif strcmp(cl,'frc') || strcmp(cl,'frcr')
                    AbsoluteError(k,c) = S.TestError{DatasetIdx}.(cl);
                else
                    AbsoluteError(k,c) = dlmread([inPath2 Dataset '_testError.dat']);
                end
            end
            for c = 1:length(Classifiers)
                if c > 1
                    NormRelativeError(k,c-1) = (AbsoluteError(k,c) - ...
                        AbsoluteError(k,1))/ChanceProb(k);
                end
            end
            k = k + 1;
        end
    end
    
    NormRelativeError(all(isnan(NormRelativeError),2),:) = [];
    
    Counts = zeros(length(BinEdges)-1,length(Classifiers)-1);

    for c = 1:length(Classifiers)-1
        Counts(:,c) = histcounts(NormRelativeError(:,c),BinEdges)';
    end
    Counts(length(BinEdges)/2,:) = sum(NormRelativeError==0);
    Counts(length(BinEdges)/2+1,:) = Counts(length(BinEdges)/2+1,:) - Counts(length(BinEdges)/2,:);
    Fractions = Counts./repmat(sum(Counts),size(Counts,1),1);

    ax(t) = axes;
    if t == 1
        YTLabel = {'RerF','RerF(r)','F-RC','Frank','RR-RF','RR-RF(r)','XGBoost'};
    else
        YTLabel = {''};
    end
    
    h = heatmap(Fractions',cellstr(num2str(BinEdges')),YTLabel,ColorMap,...
        true,'horizontal');
    if t==1
        xlabel({'Normalized Error';'Relative to RF'})
        title('Raw')
    elseif t==5
        title('Corrupted')
    else
        title(Transformations{t})
    end

    if t == 3
        cb = colorbar;
        cb.Location = 'southoutside';
        xlh = xlabel(cb,'Fraction of Datasets');
        cb.Ticks = [];
        cb.Units = 'inches';
        cb.Position = [cbLeft(t) cbBottom(t) cbWidth cbHeight];
        cb.Box = 'off';
        cb.FontSize = 16;
        xlh.Position = [0.2237 -0.5 0];
    end
    h.FontSize = FontSize;
    hold on
    for c = 2:length(Classifiers)-1
        plot(h.XLim,[c-0.5,c-0.5],'-k','LineWidth',LineWidth)
    end
    ax(t).XTick = [0.5,10,19.5];
    ax(t).XTickLabel = {'-1';'0';'1'};
    ax(t).XTickLabelRotation = 0;
    ax(t).TickLength = [0 0];
    ax(t).LineWidth = LineWidth;
    ax(t).FontUnits = 'inches';
    ax(t).FontSize = FontSize;
    ax(t).Units = 'inches';
    ax(t).Position = [axLeft(t) axBottom(t) axWidth axHeight];
        
end

save_fig(gcf,[rerfPath 'RandomerForest/Figures/pami/PAMI_fig7_error_heatmap_benchmark'],{'fig','pdf','png'})