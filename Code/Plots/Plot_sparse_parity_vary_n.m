clear
close all
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

Colors.rf = 'b';
Colors.rerfb = 'g';
Colors.rerfc = 'c';
Colors.frc2 = 'y';
Colors.frc3 = 'k';
Colors.frc4 = 'r';
Colors.rerf2 = 'm';
LineWidth = 2;

load ~/Sparse_parity_vary_n

ntrials = length(TestError{1}.rf);

for j = 1:3
    p = ps(j);
    
    Classifiers = fieldnames(TestError{1,j});

    ErrorMatrix = zeros(ntrials,length(ns),length(Classifiers));

    for i = 1:length(ns{j})
        n = ns{j}(i);

        for c = 1:length(Classifiers)
            cl = Classifiers{c};
            ErrorMatrix(:,i,c) = TestError{i,j}.(cl)';
        end
    end

    figure;
    hold on
    for c = 1:length(Classifiers)
        errorbar(ns{j},mean(ErrorMatrix(:,:,c)),std(ErrorMatrix(:,:,c))/sqrt(ntrials),...
            'LineWidth',LineWidth,'Color',Colors.(Classifiers{c}))
    end

    ax = gca;

    ax.XScale = 'log';
    ax.FontSize = 16;
    ax.XLim = [10^(log10(min(ns{j}))-0.1) 10^(log10(max(ns{j}))+0.1)];
    ax.YLim = [min(ErrorMatrix(:)) max(ErrorMatrix(:))];
    ax.XTick = ns{j};
    ax.XTickLabel = cellstr(num2str(ns{j}'))';

    xlabel('n')
    ylabel('Misclassification Rate')
    title(sprintf('Sparse Parity (p = %d)',p))
    lh = legend('RF','RerF(bin)','RerF(cont)','F-RC(L=2)','F-RC(L=3)','F-RC(L=4)','RerF2');
    lh.Box = 'off';

%     save_fig(gcf,[rerfPath sprintf('RandomerForest/Figures/Sparse_parity_p_%d',p)])
end