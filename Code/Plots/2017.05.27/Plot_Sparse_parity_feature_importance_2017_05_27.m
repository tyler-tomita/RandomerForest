close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

if isempty(rerfPath)
    rerfPath = '~/';
end

S = load([rerfPath 'RandomerForest/Results/2017.05.27/Feature_importance/Sparse_parity_feature_importance_permutation_2017_05_27.mat']);

plot(S.importance(1:20),'-o','MarkerSize',10,'LineWidth',2)
title({'Sparse Parity';'Top 20 Ranked Features'})
xlabel('Feature Index')
ylabel('Permutation Importance')

save_fig(gcf,[rerfPath 'RandomerForest/Figures/2017.05.27/Sparse_parity_feature_importance_permutation_2017_05_27'],{'fig' 'pdf' 'png'})

S = load([rerfPath 'RandomerForest/Results/2017.05.27/Feature_importance/Sparse_parity_feature_importance_gini_2017_05_27.mat']);

figure;
plot(S.importance(1:20),'-o','MarkerSize',10,'LineWidth',2)
title({'Sparse Parity';'Top 20 Ranked Features'})
xlabel('Feature Index')
ylabel('Gini Importance')

save_fig(gcf,[rerfPath 'RandomerForest/Figures/2017.05.27/Sparse_parity_feature_importance_gini_2017_05_27'],{'fig' 'pdf' 'png'})