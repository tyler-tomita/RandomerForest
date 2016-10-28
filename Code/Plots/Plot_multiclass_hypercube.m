clear
close all
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

Colors.rf = 'b';
Colors.rerf = 'g';
Colors.frc = 'c';
Colors.rr_rf = 'm';
LineWidth = 2;

load Multiclass_hypercube

ntrials = length(TestError{1}.rf);
Classifiers = fieldnames(TestError{1});

ErrorMatrix = zeros(ntrials,length(ns),length(Classifiers));

for i = 1:length(ns)
    n = ns(i);
    
    for c = 1:length(Classifiers)
        cl = Classifiers{c};
        ErrorMatrix(:,i,c) = TestError{i}.(cl)';
    end
end

for c = 1:length(Classifiers)
    errorbar(ns,mean(ErrorMatrix(:,:,c)),std(ErrorMatrix(:,:,c))/sqrt(ntrials),...
        'LineWidth',LineWidth,'Color',Colors.(Classifiers{c}))
    hold on
end

ax = gca;

ax.XScale = 'log';
ax.FontSize = 16;
ax.XLim = [90 2000];
ax.XTick = ns;
ax.XTickLabel = {'100','200','400','800','1600'};

xlabel('n')
ylabel('Misclassification Rate')
title('Multiclass Hypercube')
lh = legend('RF','RerF','F-RC','RR-RF');
lh.Box = 'off';

save_fig(gcf,[rerfPath 'RandomerForest/Figures/Multiclass_hypercube'])