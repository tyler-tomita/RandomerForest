%% Plot benchmark classifier rank distributions

clear
close all
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

load('purple2green')
Colors.rf = ColorMap(2,:);
Colors.rerf = 'k';
Colors.frc= ColorMap(10,:);
Colors.rr_rf = ColorMap(4,:);
Colors.xgb= ColorMap(8,:);
LineWidth = 2;
MarkerSize = 8;
FontSize = .2;
axWidth = 2;
axHeight = 2;
cbWidth = axWidth;
cbHeight = 0.15;
axBottom = FontSize*5*ones(1,5);
axLeft = fliplr([FontSize*9+axHeight*4,FontSize*8+axHeight*3,...
    FontSize*7+axHeight*2,FontSize*6+axHeight,...
    FontSize*5]);
legWidth = axWidth*0.75;
legHeight = axHeight;
legLeft = axLeft(end) + axWidth;
legBottom = axBottom(end);
figWidth = legLeft + legWidth;
figHeight = axBottom(1) + axHeight + FontSize*2.5;

% fig = figure;
% fig.Units = 'inches';
% fig.PaperUnits = 'inches';
% fig.Position = [0 0 figWidth figHeight];
% fig.PaperPosition = [0 0 figWidth figHeight];
% fig.PaperSize = [figWidth figHeight];

Classifiers = {'rf','rerf','frc','rr_rf'};

inPath1 = [rerfPath 'RandomerForest/Results/2017.04.01/Benchmarks/Raw/'];
inPath2 = '~/Benchmarks/Results/R/dat/Raw/';
contents = dir([inPath1 '*.mat']);

AbsoluteError = NaN(length(contents),length(Classifiers));
NormalizedRelativeError = NaN(length(contents),length(Classifiers)-1);
ChanceProb = NaN(length(contents),1);
nClasses = NaN(length(contents),1);
d = NaN(length(contents),length(Classifiers));
p = NaN(length(contents),1);
p_lowrank = NaN(length(contents),1);
lambda = cell(length(contents),1);   % eigenvalues of data matrix
ntrain = NaN(length(contents),1);
TraceNorm = NaN(length(contents),1);

DatasetNames = importdata('~/Benchmarks/Data/Names.txt');

S = load('~/Benchmarks/Results/Benchmark_untransformed.mat');

k = 1;

for i = 1:length(contents)
    fprintf('Dataset %d\n',i)
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
        TrainSet = dlmread(['~/Benchmarks/Data/dat/Raw/' Dataset '_train.dat']);
        TestSet = dlmread(['~/Benchmarks/Data/dat/Raw/' Dataset '_test.dat']);
        [ntrain(k),p(k)] = size(TrainSet(:,1:end-1));
        [coeff,score,lambda{k}] = pca([TrainSet(:,1:end-1);TestSet(:,1:end-1)]);
        TraceNorm(k) = sum(sqrt(lambda{k}));
        VarCutoff = 0.9;
        AboveCutoff = find(cumsum(lambda{k})/sum(lambda{k}) >= 0.9);
        p_lowrank(k) = AboveCutoff(1);
        nClasses(k) = length(unique(TestSet(:,end)));
        ClassCounts = histcounts(TestSet(:,end),nClasses(k));
        ChanceProb(k) = 1 - max(ClassCounts)/sum(ClassCounts);

        for c = 1:length(Classifiers)
            cl = Classifiers{c};
            if ~strcmp(cl,'xgb') && ~strcmp(cl,'frc') && ~strcmp(cl,'frcr')
                BI = hp_optimize(OOBError.(cl)(end,1:length(Params.(cl).d)),...
                    OOBAUC.(cl)(end,1:length(Params.(cl).d)));
                BI = BI(randperm(length(BI),1));
                AbsoluteError(k,c) = TestError.(cl)(BI);
                d(k,c) = Params.(cl).d(BI);
            elseif strcmp(cl,'frc') || strcmp(cl,'frcr')
                BI = hp_optimize(S.OOBError{DatasetIdx}.(cl),S.OOBAUC{DatasetIdx}.(cl));
                BI = BI(end);                
                AbsoluteError(k,c) = S.TestError{DatasetIdx}.(cl);
                d(k,c) = S.Params{DatasetIdx}.(cl).d(BI);
            else
                AbsoluteError(k,c) = dlmread([inPath2 Dataset '_testError.dat']);
            end
            
%             if ~strcmp(cl,'xgb')
%                 BI = hp_optimize(OOBError.(cl)(end,1:length(Params.(cl).d)),...
%                     OOBAUC.(cl)(end,1:length(Params.(cl).d)));
%                 BI = BI(end);
%                 AbsoluteError(k,c) = TestError.(cl)(BI);
%                 d(k,c) = Params.(cl).d(BI);
%             else
%                 AbsoluteError(k,c) = dlmread([inPath2 Dataset '_testError.dat']);
%             end
            if c > 1
                NormalizedRelativeError(k,c-1) = (AbsoluteError(k,c)-AbsoluteError(k,1))/ChanceProb(k);
            end
        end
        k = k + 1;
    end
end

AbsoluteError(all(isnan(NormalizedRelativeError),2),:) = [];
NormalizedRelativeError(all(isnan(NormalizedRelativeError),2),:) = [];
d(all(isnan(d),2),:) = [];
p(isnan(p)) = [];
p_lowrank(isnan(p_lowrank)) = [];
ntrain(isnan(ntrain)) = [];
TraceNorm(isnan(TraceNorm)) = [];
lambda(isnan(p_lowrank)) = [];
nClasses(isnan(nClasses)) = [];

[~,srtidx] = sort(NormalizedRelativeError(:,1));
Top15 = srtidx(1:15);
threshold = 0;
RerFIsBetter = NormalizedRelativeError(:,1)<threshold;

% scree plots for top 15 RerF-favored datasets

% figure;
% for i = 1:length(Top15)
%     axh(i) = subplot(5,3,i);
%     csum = cumsum(lambda{Top15(i)});
%     plot(csum/csum(end))
%     xlabel('no. PCs retained')
%     ylabel('explained variance')
%     axis square
%     axh(i).FontSize = 6;
% end
% save_fig(gcf,[rerfPath 'RandomerForest/Figures/pami/benchmark_scree_plots'],{'fig','pdf','png'})

% figure;
% k = 1;
% ax(k) = axes;
% hold on
% for c = 2:length(Classifiers)
%     cl = Classifiers{c};
%     plot(p(RerFIsBetter),NormalizedRelativeError(RerFIsBetter,c-1),'.','MarkerSize',MarkerSize,...
%         'Color',Colors.(cl),'MarkerSize',14)
% end
% ax(k).XScale = 'log';
% ax(k).Box = 'off';
% xlabel('p')
% ylabel('Normalized Error Relative to RF')
% title('Datasets for Which RerF Wins')
% l = legend('RerF','F-RC','RR-RF');
% % l.Box = 'off';
% save_fig(gcf,[rerfPath 'RandomerForest/Figures/pami/pami_benchmark_p_vs_error'],{'fig','pdf','png'})
% 
% figure;
% k = 2;
% ax(k) = axes;
% hold on
% for c = 2:length(Classifiers)
%     cl = Classifiers{c};
%     plot(d(RerFIsBetter,c)./p(RerFIsBetter),NormalizedRelativeError(RerFIsBetter,c-1),'.','MarkerSize',MarkerSize,...
%         'Color',Colors.(cl),'MarkerSize',14)
% end
% ax(k).XScale = 'log';
% ax(k).Box = 'off';
% xlabel('d/p')
% ylabel('Normalized Error Relative to RF')
% l = legend('RerF','F-RC','RR-RF');
% % l.Box = 'off';
% title('Datasets for Which RerF Wins')
% save_fig(gcf,[rerfPath 'RandomerForest/Figures/pami/pami_benchmark_d_over_p_vs_error'],{'fig','pdf','png'})
% 
% figure;
% k = 3;
% ax(k) = axes;
% hold on
% for c = 2:length(Classifiers)
%     cl = Classifiers{c};
%     plot(ntrain(RerFIsBetter),NormalizedRelativeError(RerFIsBetter,c-1),'.','MarkerSize',MarkerSize,...
%         'Color',Colors.(cl),'MarkerSize',14)
% end
% ax(k).XScale = 'log';
% ax(k).Box = 'off';
% xlabel('n_{train}')
% ylabel('Normalized Error Relative to RF')
% l = legend('RerF','F-RC','RR-RF');
% % l.Box = 'off';
% title('Datasets for Which RerF Wins')
% save_fig(gcf,[rerfPath 'RandomerForest/Figures/pami/pami_benchmark_ntrain_vs_error'],{'fig','pdf','png'})
% 
% figure;
% k = 4;
% ax(k) = axes;
% hold on
% for c = 2:length(Classifiers)
%     cl = Classifiers{c};
%     plot(p(RerFIsBetter)./ntrain(RerFIsBetter),NormalizedRelativeError(RerFIsBetter,c-1),'.','MarkerSize',MarkerSize,...
%         'Color',Colors.(cl),'MarkerSize',14)
% end
% ax(k).XScale = 'log';
% ax(k).Box = 'off';
% xlabel('p/n_{train}')
% ylabel('Normalized Error Relative to RF')
% l = legend('RerF','F-RC','RR-RF');
% % l.Box = 'off';
% title('Datasets for Which RerF Wins')
% save_fig(gcf,[rerfPath 'RandomerForest/Figures/pami/pami_benchmark_p_over_ntrain_vs_error'],{'fig','pdf','png'})
% 
% figure;
% k = 5;
% ax(k) = axes;
% hold on
% for c = 2:length(Classifiers)
%     cl = Classifiers{c};
%     plot(p_lowrank(RerFIsBetter),NormalizedRelativeError(RerFIsBetter,c-1),'.','MarkerSize',MarkerSize,...
%         'Color',Colors.(cl),'MarkerSize',14)
% end
% ax(k).XScale = 'log';
% ax(k).Box = 'off';
% xlabel('p^*')
% ylabel('Normalized Error Relative to RF')
% l = legend('RerF','F-RC','RR-RF');
% % l.Box = 'off';
% title('Datasets for Which RerF Wins')
% save_fig(gcf,[rerfPath 'RandomerForest/Figures/pami/pami_benchmark_pstar_vs_error'],{'fig','pdf','png'})
% 
% figure;
% k = 6;
% ax(k) = axes;
% hold on
% for c = 2:length(Classifiers)
%     cl = Classifiers{c};
%     plot(p_lowrank(RerFIsBetter)./p(RerFIsBetter),NormalizedRelativeError(RerFIsBetter,c-1),'.','MarkerSize',MarkerSize,...
%         'Color',Colors.(cl),'MarkerSize',14)
% end
% ax(k).XScale = 'log';
% ax(k).Box = 'off';
% xlabel('p^*/p')
% ylabel('Normalized Error Relative to RF')
% l = legend('RerF','F-RC','RR-RF');
% % l.Box = 'off';
% title('Datasets for Which RerF Wins')
% save_fig(gcf,[rerfPath 'RandomerForest/Figures/pami/pami_benchmark_pstar_over_p_vs_error'],{'fig','pdf','png'})
% 
% figure;
% k = 7;
% ax(k) = axes;
% hold on
% for c = 2:length(Classifiers)
%     cl = Classifiers{c};
%     plot(TraceNorm(RerFIsBetter),NormalizedRelativeError(RerFIsBetter,c-1),'.','MarkerSize',MarkerSize,...
%         'Color',Colors.(cl),'MarkerSize',14)
% end
% ax(k).XScale = 'log';
% ax(k).Box = 'off';
% xlabel('Trace-norm')
% ylabel('Normalized Error Relative to RF')
% l = legend('RerF','F-RC','RR-RF');
% % l.Box = 'off';
% title('Datasets for Which RerF Wins')
% save_fig(gcf,[rerfPath 'RandomerForest/Figures/pami/pami_benchmark_trace_norm_vs_error'],{'fig','pdf','png'})

figure;
k = 8;
ax(k) = axes;
hold on
for c = 2:length(Classifiers)
    cl = Classifiers{c};
    plot(nClasses,NormalizedRelativeError(:,c-1),'.','MarkerSize',MarkerSize,...
        'Color',Colors.(cl),'MarkerSize',14)
end
ax(k).XScale = 'log';
ax(k).Box = 'off';
xlabel('Number of Classes')
ylabel('Normalized Error Relative to RF')
l = legend('RerF','F-RC','RR-RF');
% l.Box = 'off';
title('Datasets for Which RerF Wins')
% save_fig(gcf,[rerfPath 'RandomerForest/Figures/pami/pami_benchmark_nclasses_vs_error'],{'fig','pdf','png'})