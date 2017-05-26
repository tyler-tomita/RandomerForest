close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

if isempty(rerfPath)
    rerfPath = '~/';
end

S = load([rerfPath 'RandomerForest/Results/2017.05.27/Sparse_parity_feature_importance_2017_05_27.mat']);

plot(S.importance(1:20),'-o','MarkerSize',10,'LineWidth',2)
title({'Sparse Parity';'Top 25 Ranked Features'})
xlabel('Feature Index')
ylabel('Importance')

save_fig(gcf,[rerfPath 'RandomerForest/Results/2017.05.27/Sparse_parity_feature_importance_2017_05_27'],{'fig' 'pdf' 'png'})