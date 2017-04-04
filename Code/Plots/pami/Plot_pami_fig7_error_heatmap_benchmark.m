%% Plot benchmark classifier rank distributions

clear
close all
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

load('purple2green')
ColorMap = interpolate_colormap(ColorMap(round(size(ColorMap,1)/2):end,:),64,false);
% ColorMap = flipud(interpolate_colormap(ColorMap(3:end-2,:),64,true));

LineWidth = 1.5;
FontSize = 0.1;
axWidth = 2.4;
axHeight = 1.6;
cbWidth = .1;
cbHeight = axHeight;
axLeft = FontSize*5.5*ones(1,5);
axBottom = [FontSize*14+axHeight*4,FontSize*11+axHeight*3,...
    FontSize*8+axHeight*2,FontSize*5+axHeight,...
    FontSize*2];
cbLeft = axLeft + axWidth + FontSize/2;
cbBottom = axBottom;
figWidth = cbLeft(1) + cbWidth + FontSize*3;
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

for t = 1:length(Transformations)
    if t ==1
        Classifiers = {'rf','rerf','rr_rf','xgb'};
    else
        Classifiers = {'rf','rerf','rerfr','rr_rf','rr_rfr','xgb'};
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
            TrainSet = dlmread(['~/Benchmarks/Data/dat/' Transformations{t} '/' Dataset '_train.dat']);
            [ntrain,p] = size(TrainSet(:,1:end-1));
            nClasses = length(unique(TrainSet(:,end)));
            ClassCounts = histcounts(TrainSet(:,end),nClasses);
            ChanceProb(k) = 1 - max(ClassCounts)/sum(ClassCounts);
            
            for c = 1:length(Classifiers)
                cl = Classifiers{c};
                if ~strcmp(cl,'xgb')
                    if t > 1
                        TestError.(cl) = TestError.(cl)(BestIdx.(cl));
                    end
                    AbsoluteError(k,c) = TestError.(cl);
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
        h = heatmap(flipud(Fractions),{'RerF','RR-RF','XGBoost'},cellstr(num2str(flipud(BinEdges'))),ColorMap,true);
    else
        h = heatmap(flipud(Fractions),{'RerF','RerF(r)','RR-RF','RR-RF(r)','XGBoost'},cellstr(num2str(flipud(BinEdges'))),ColorMap,true);
    end
    if t==1
        ylabel('Normalized Error Relative to RF')
        title('Raw','FontSize',8)
    elseif t==5
        title('Corrupted','FontSize',8)
    else
        title(Transformations{t},'FontSize',8)
    end
    cb = colorbar;
    cb.Units = 'inches';
    cb.Position = [cbLeft(t) cbBottom(t) cbWidth cbHeight];
    cb.Box = 'off';
    cb.FontSize = 8;
    h.FontSize = FontSize;
    hold on
    for c = 2:length(Classifiers)-1
        plot([c-0.5,c-0.5],h.YLim,'-k','LineWidth',LineWidth+1)
    end
    ax(t).LineWidth = LineWidth;
    ax(t).FontUnits = 'inches';
    ax(t).FontSize = FontSize;
    ax(t).Units = 'inches';
    ax(t).Position = [axLeft(t) axBottom(t) axWidth axHeight];
        
end

save_fig(gcf,[rerfPath 'RandomerForest/Figures/pami/PAMI_fig7_error_heatmap_benchmark'],{'fig','pdf','png'})