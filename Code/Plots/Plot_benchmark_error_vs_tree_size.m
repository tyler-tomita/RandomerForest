clear
close all
clc

FontSize = 14;

S = load('Benchmark_untransformed_n_lt_5000_datasets_1_13_2016_11_28.mat');

Colors = get(gca,'ColorOrder');

TreeDepth = [];
nNodes = [];
nSplits = [];
ErrorRate = [];
Clr = [];

for i = 1:length(S.TestError)
    if ~isempty(S.TestError{i})
        Classifiers = fieldnames(S.TestError{i});
        Classifiers(~ismember(Classifiers,{'rf','rerf','rerf2','frc','rr_rf'})) = [];
        for c = 1:length(Classifiers)
            BestIdx = hp_optimize(S.OOBError{i}.(Classifiers{c})(end,:),...
                S.OOBAUC{i}.(Classifiers{c})(end,:));
            if length(BestIdx)>1
                BestIdx = BestIdx(end);
            end
            TreeDepth = [TreeDepth,mean(S.Depth{i}.(Classifiers{c})(:,BestIdx))];
            nNodes = [nNodes,mean(S.Depth{i}.(Classifiers{c})(:,BestIdx))];
            nSplits = [nSplits,mean(S.Depth{i}.(Classifiers{c})(:,BestIdx))];
            ErrorRate = [ErrorRate,S.TestError{i}.(Classifiers{c})];
            Clr = [Clr;Colors(c,:)];
        end
    end
end
% plot(TreeDepth,ErrorRate,Markers{i},'MarkerEdgeColor',Colors(c,:))
gscatter(TreeDepth,ErrorRate,Clr)
save_fig(gcf,'~/RandomerForest/Figures/Benchmark_tree_size')

xlabel('Tree Depth')
ylabel('Error Rate')
title('Benchmarks')
legend('RF','RerF','RerF2','F-RC','RR-RF','Location','southeast')
ax = gca;
ax.YScale = 'log';
ax.FontSize = FontSize;