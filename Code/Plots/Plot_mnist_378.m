%% Plot out of bag error vs n on mnist dataset

clear
close all
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

load mnist_378_vary_n

fn = fieldnames(Lhat);

keys = {'srerf' 'rerf' 'rf'};
values = {'Structured RerF' 'RerF' 'RF'};
clNames = containers.Map(keys,values);

for i = 1:length(clNames.keys)
    cl = clNames.keys{i};
    meanError = NaN(1,length(ns));
    semError = NaN(1,length(ns));
    for j = 1:length(ns)
        minError = min(Lhat.(cl){j},[],2);
        meanError(j) = mean(minError);
        semError(j) = std(minError)/sqrt(length(minError));
    end
    errorbar(ns,meanError,semError,'LineWidth',2)
    hold on
end
legend(clNames.values)
ax = gca;
ax.XScale = 'log';
xlabel('n')
ylabel('Out-of-Bag Error')
title('MNIST (Digits 3, 7, & 8)')

save_fig(gcf,[rerfPath 'RandomerForest/Figures/MNIST'])