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

load Multiparity
load Multiparity_bayes_error

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
ax.XLim = [40 1500];
ax.XTick = [50 100 500 1000];
ax.XTickLabel = {'50','100','500','1000'};

plot(ax.XLim,BayesError*ones(1,2),'--k','LineWidth',LineWidth)

xlabel('n')
ylabel('Misclassification Rate')
lh = legend('RF','RerF','F-RC','RR-RF','Bayes');
lh.Box = 'off';

save_fig(gcf,[rerfPath 'RandomerForest/Figures/Multiparity'])