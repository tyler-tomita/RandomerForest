clear
close all
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

load('purple2green')
Colors.rf = ColorMap(1,:);
Colors.rerf = ColorMap(8,:);
Colors.rr_rf = ColorMap(4,:);
LineWidth = 2;
FontSize = .15;
axWidth = 1.5;
axHeight = 1.5;
axLeft = [FontSize*4,FontSize*8+axWidth,FontSize*12+2*axWidth];
axBottom = FontSize*3*ones(1,3);
legWidth = axWidth;
legHeight = axHeight;
legLeft = axLeft(end) + axWidth + FontSize/2;
legBottom = axBottom(end);
figWidth = legLeft(end) + legWidth + FontSize*3;
figHeight = axBottom(1) + axHeight + FontSize*2;

fig = figure;
fig.Units = 'inches';
fig.PaperUnits = 'inches';
fig.Position = [0 0 figWidth figHeight];
fig.PaperPosition = [0 0 figWidth figHeight];
fig.PaperSize = [figWidth figHeight];

load ~/RandomerForest/Results/2017.02.19/Orthant_aggregated_results_2017_02_19

ntrials = length(TestError{1}.rf);

for j = 1:length(ps)
    p = ps(j);
    
    Classifiers = fieldnames(TestError{1,j});

    ErrorMatrix = NaN(ntrials,length(ns),length(Classifiers));

    for i = 1:length(ns{j})
        n = ns{j}(i);

        for c = 1:length(Classifiers)
            cl = Classifiers{c};
            if ~isempty(TestError{i,j}.(cl))
                ErrorMatrix(:,i,c) = TestError{i,j}.(cl)';
            end
        end
    end

    ax(j) = axes;
    hold on
    for c = 1:length(Classifiers)
        errorbar(ns{j},mean(ErrorMatrix(:,:,c)),std(ErrorMatrix(:,:,c))/sqrt(ntrials),...
            'LineWidth',LineWidth,'Color',Colors.(Classifiers{c}))
    end

    ax(j).XScale = 'log';
    ax(j).FontSize = 16;
    ax(j).XLim = [10^(log10(min(ns{j}))-0.1) 10^(log10(max(ns{j}))+0.1)];
    ax(j).YLim = [min(ErrorMatrix(:)) max(ErrorMatrix(:))];
    ax(j).XTick = ns{j};
    ax(j).XTickLabel = cellstr(num2str(ns{j}'))';
    ax(j).LineWidth = LineWidth;
    ax(j).FontUnits = 'inches';
    ax(j).FontSize = FontSize;
    ax(j).Units = 'inches';
    ax(j).Position = [axLeft(j) axBottom(j) axWidth axHeight];

    xlabel('n')
    ylabel('Error Rate')
    title(sprintf('Orthant (p = %d)',p))
    if j ==1
        [lh,objh] = legend('RF','RerF+','RR-RF');
        lh.Units = 'inches';
        lh.Position = [legLeft legBottom legWidth legHeight];
        lh.Box = 'off';
        for i = 4:length(objh)
            objh(i).Children.Children(2).XData = [(objh(i).Children.Children(2).XData(2)-objh(i).Children.Children(2).XData(1))*0.5+objh(i).Children.Children(2).XData(1),objh(i).Children.Children(2).XData(2)];
        end
    end
end

save_fig(gcf,[rerfPath 'RandomerForest/Figures/2017.02.19/Orthant_vary_n_2017_02_19'],{'fig','pdf','png'})