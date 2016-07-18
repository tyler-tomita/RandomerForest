%% Plot out of bag error vs n on mnist dataset

clear
close all
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

load mnist_378_vary_n

fn = fieldnames(Lhat);

fn = {'srerf' 'rerf' 'rf' 'control'};
clNames = {'Structured RerF' 'RerF' 'RF' 'Control'};

Colors.srerf = 'g';
Colors.rerf = 'c';
Colors.rf = 'm';
Colors.control = 'k';

for i = 1:length(fn)-1
    cl = fn{i};
    meanError = NaN(1,length(ns));
    semError = NaN(1,length(ns));
    for j = 1:length(ns)
        minError = min(Lhat.(cl){j},[],2);
        meanError(j) = mean(minError);
        semError(j) = std(minError)/sqrt(length(minError));
    end
    errorbar(ns,meanError,semError,'LineWidth',2,'Color',Colors.(cl))
    hold on
end

load mnist_378_vary_n_control

i = i + 1;

cl = fn{i};
meanError = NaN(1,length(ns));
semError = NaN(1,length(ns));
for j = 1:length(ns)
    minError = min(Lhat.(cl){j},[],2);
    meanError(j) = mean(minError);
    semError(j) = std(minError)/sqrt(length(minError));
end
errorbar(ns,meanError,semError,'LineWidth',2,'Color',Colors.(cl))
hold on


l = legend(clNames);
l.Box = 'off';
ax = gca;
ax.XScale = 'log';
ax.XLim = [10^(log10(min(ns))-.1),10^(log10(max(ns))+.1)];
ax.XTick = ns;
ax.XTickLabel = cellstr(num2str(ns'))';
xlabel('n','FontSize',14)
ylabel('Out-of-Bag Error','FontSize',14)
ax.Box = 'off';
ax.FontSize = 14;

save_fig(gcf,[rerfPath 'RandomerForest/Figures/MNIST_378'])