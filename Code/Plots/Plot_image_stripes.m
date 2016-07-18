%% Plot out of bag error vs n on mnist dataset

clear
close all
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

load ('image_stripes_data','ns')
load image_stripes

fn = fieldnames(Lhat);

fn = {'srerf' 'rerf' 'rf' 'control'};
clNames = {'Structured RerF' 'RerF' 'RF' 'Control'};

Colors.srerf = 'g';
Colors.rerf = 'c';
Colors.rf = 'm';
Colors.control = 'k';

for i = 1:length(fn)
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

l = legend(clNames);
l.Box = 'off';
ax = gca;
ax.XScale = 'log';
ax.XLim = [10^(log10(min(ns))-0.1),10^(log10(max(ns))+0.1)];
ax.XTick = [10 20 50];
ax.XTickLabel = cellstr(num2str(ax.XTick'))';
xlabel('n')
ylabel('Out-of-Bag Error')
ax.Box = 'off';
ax.FontSize = 14;

save_fig(gcf,[rerfPath 'RandomerForest/Figures/Image_stripes'])