clear
close all
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

load('purple2green')

Colors.rf = ColorMap(2,:);
Colors.rerf = 'k';
Colors.frc= ColorMap(10,:);
Colors.rr_rf = ColorMap(4,:);
Colors.xgb= ColorMap(8,:);
LineWidth = 2;
MarkerSize = 12;
FontSize = .2;
axWidth = 2;
axHeight = 1.5;
axLeft = repmat([FontSize*4,FontSize*8+axWidth,FontSize*12+axWidth*2],1,3);
axBottom = [...
    (FontSize*14+axHeight*2)*ones(1,3),(FontSize*9+axHeight)*ones(1,3),...
    FontSize*4*ones(1,3)];
legWidth = axWidth;
legHeight = axHeight;
legLeft = axLeft(end) + axWidth + FontSize;
legBottom = axBottom(5);
figWidth = legLeft + legWidth;
figHeight = axBottom(1) + axHeight + FontSize*2.5;

fig = figure;
fig.Units = 'inches';
fig.PaperUnits = 'inches';
fig.Position = [0 0 figWidth figHeight];
fig.PaperPosition = [0 0 figWidth figHeight];
fig.PaperSize = [figWidth figHeight];

% Plot Sparse Parity

k = 1;

% load([rerfPath 'RandomerForest/Results/pami/Sparse_parity/mat/Sparse_parity_vary_n.mat'])
load([rerfPath 'RandomerForest/Results/2017.04.13/Sparse_parity/Sparse_parity_raw_vary_n_aggregated.mat'])

ntrials = length(TestError{1}.rf);
Classifiers = {'rf','rerf','frc','rr_rf','xgb'};

for j = 1:3
    p = ps(j);
    
%     Classifiers = [fieldnames(TestError{1,j});'xgb'];
    
    TimeMatrix = NaN(ntrials,length(ns),length(Classifiers));

    for i = 1:length(ns{j})
        n = ns{j}(i);

        for c = 1:length(Classifiers)
            cl = Classifiers{c};
            if ~strcmp(cl,'xgb')
                if ~isempty(TrainTime{i,j}.(cl))
                    OE = OOBError{i,j}.(cl);
                    OA = OOBAUC{i,j}.(cl);
                    ntrials = size(OE,1);
                    T = NaN(ntrials,1);
                    for trial = 1:ntrials
                        BI = hp_optimize(OE(trial,1:length(Params{i,j}.(cl).d)),OA(trial,1:length(Params{i,j}.(cl).d)));
                        T(trial) = TrainTime{i,j}.(cl)(trial,BI(end));
                    end
%                     TimeMatrix(:,i,c) = T;
                    TimeMatrix(:,i,c) = sum(TrainTime{i,j}.(cl),2);
%                     TimeMatrix(:,i,c) = TrainTime{i,j}.(cl)';
                else
                    TimeMatrix(:,i,c) = NaN;
                end
            else
                try
                    TimeMatrix(:,i,c) = ...
                        dlmread([rerfPath ...
                        sprintf('RandomerForest/Results/pami/Sparse_parity/dat/Sparse_parity_vary_n_trainTime_n%d_p%d.dat',n,p)]);
                catch ME
                    TimeMatrix(:,i,c) = NaN;
                end
            end
        end
    end

    ax((j-1)*3+k) = axes;
    hold on
    ymax = zeros(length(Classifiers),1);
    ymin = zeros(length(Classifiers),1);
    for c = 1:length(Classifiers)
        mn = mean(TimeMatrix(:,:,c));
        sem = std(TimeMatrix(:,:,c))/sqrt(ntrials);
        errorbar(ns{j},mn,sem,...
            'LineWidth',LineWidth,'Color',Colors.(Classifiers{c}))
        ymax(c) = mn(end)+2*sem(end);
        ymin(c) = mn(1)-2*sem(1);
    end

    ax((j-1)*3+k).LineWidth = LineWidth;
    ax((j-1)*3+k).FontUnits = 'inches';
    ax((j-1)*3+k).FontSize = FontSize;
    ax((j-1)*3+k).Units = 'inches';
    ax((j-1)*3+k).Position = [axLeft((j-1)*3+k) axBottom((j-1)*3+k) axWidth axHeight];
    ax((j-1)*3+k).Box = 'off';
    ax((j-1)*3+k).XLim = [10^(log10(min(ns{j}))-0.1) 10^(log10(max(ns{j}))+0.1)];
    ax((j-1)*3+k).XScale = 'log';
    ax((j-1)*3+k).XTick = ns{j};
    ax((j-1)*3+k).XTickLabel = cellstr(num2str(ns{j}'))';
    ax((j-1)*3+k).XTickLabelRotation = 45;
    ax((j-1)*3+k).YLim = [min(TimeMatrix(:)) max(TimeMatrix(:))];
    
%     ax((j-1)*3+k).XScale = 'log';
%     ax((j-1)*3+k).FontSize = 16;
%     ax((j-1)*3+k).XLim = [10^(log10(min(ns{j}))-0.1) 10^(log10(max(ns{j}))+0.1)];
%     ax((j-1)*3+k).YLim = [min(TimeMatrix(:)) max(TimeMatrix(:))];
%     ax((j-1)*3+k).XTick = ns{j};
%     ax((j-1)*3+k).XTickLabel = cellstr(num2str(ns{j}'))';

    xlabel('n')
    ylabel('Train Time (s)')
    if j == 1
        title({'Sparse Parity';sprintf('p = %d',p)})
    else
        title(sprintf('p = %d',p))
    end
%     lh = legend('RF','RerF','RR-RF');
    lh.Box = 'off';

end

% Plot Trunk

k = 2;

% load([rerfPath 'RandomerForest/Results/pami/Trunk/mat/Trunk_vary_n.mat'])
load([rerfPath 'RandomerForest/Results/2017.04.12/Trunk/Trunk_raw_vary_n_aggregated.mat'])

ntrials = length(TestError{1}.rf);

for j = 1:3
    p = ps(j);
    
%     Classifiers = [fieldnames(TestError{1,j});'xgb'];
    
    TimeMatrix = NaN(ntrials,length(ns),length(Classifiers));

    for i = 1:length(ns{j})
        n = ns{j}(i);

        for c = 1:length(Classifiers)
            cl = Classifiers{c};
            if ~strcmp(cl,'xgb')
                if ~isempty(TrainTime{i,j}.(cl))
                    OE = OOBError{i,j}.(cl);
                    OA = OOBAUC{i,j}.(cl);
                    ntrials = size(OE,1);
                    T = NaN(ntrials,1);
                    for trial = 1:ntrials
                        BI = hp_optimize(OE(trial,1:length(Params{i,j}.(cl).d)),OA(trial,1:length(Params{i,j}.(cl).d)));
                        T(trial) = TrainTime{i,j}.(cl)(trial,BI(end));
                    end
%                     TimeMatrix(:,i,c) = T;
                    TimeMatrix(:,i,c) = sum(TrainTime{i,j}.(cl),2);
%                     TimeMatrix(:,i,c) = TrainTime{i,j}.(cl)';
                else
                    TimeMatrix(:,i,c) = NaN;
                end
            else
                try
                    TimeMatrix(:,i,c) = ...
                        dlmread([rerfPath ...
                        sprintf('RandomerForest/Results/pami/Trunk/dat/Trunk_vary_n_trainTime_n%d_p%d.dat',n,p)]);
                catch ME
                    TimeMatrix(:,i,c) = NaN;
                end
            end
        end
    end

    ax((j-1)*3+k) = axes;
    hold on
    ymax = zeros(length(Classifiers),1);
    ymin = zeros(length(Classifiers),1);
    for c = 1:length(Classifiers)
        mn = mean(TimeMatrix(:,:,c));
        sem = std(TimeMatrix(:,:,c))/sqrt(ntrials);
        errorbar(ns{j},mn,sem,...
            'LineWidth',LineWidth,'Color',Colors.(Classifiers{c}))
        ymax(c) = mn(end)+2*sem(end);
        ymin(c) = mn(1)-2*sem(1);
    end

    ax((j-1)*3+k).LineWidth = LineWidth;
    ax((j-1)*3+k).FontUnits = 'inches';
    ax((j-1)*3+k).FontSize = FontSize;
    ax((j-1)*3+k).Units = 'inches';
    ax((j-1)*3+k).Position = [axLeft((j-1)*3+k) axBottom((j-1)*3+k) axWidth axHeight];
    ax((j-1)*3+k).Box = 'off';
    ax((j-1)*3+k).XLim = [10^(log10(min(ns{j}))-0.1) 10^(log10(max(ns{j}))+0.1)];
    ax((j-1)*3+k).XScale = 'log';
    ax((j-1)*3+k).XTick = ns{j};
    ax((j-1)*3+k).XTickLabel = cellstr(num2str(ns{j}'))';
    ax((j-1)*3+k).XTickLabelRotation = 45;
    ax((j-1)*3+k).YLim = [min(TimeMatrix(:)) max(TimeMatrix(:))];
    
%     ax((j-1)*3+k).XScale = 'log';
%     ax((j-1)*3+k).FontSize = 16;
%     ax((j-1)*3+k).XLim = [10^(log10(min(ns{j}))-0.1) 10^(log10(max(ns{j}))+0.1)];
%     ax((j-1)*3+k).YLim = [min(TimeMatrix(:)) max(TimeMatrix(:))];
%     ax((j-1)*3+k).XTick = ns{j};
%     ax((j-1)*3+k).XTickLabel = cellstr(num2str(ns{j}'))';

    xlabel('n')
    ylabel('Train Time (s)')
    if j == 1
        title({'Trunk';sprintf('p = %d',p)})
    else
        title(sprintf('p = %d',p))
    end
    
    if j == 2
        lh = legend('RF','RerF','F-RC','RR-RF','XGBoost');
        lh.Box = 'off';
        lh.Units = 'inches';
        lh.Position = [legLeft legBottom legWidth legHeight];
    end

end

% Plot Orthant

k = 3;

% load([rerfPath 'RandomerForest/Results/pami/Orthant/mat/Orthant_vary_n.mat'])
load([rerfPath 'RandomerForest/Results/2017.04.13/Orthant/Orthant_raw_vary_n_aggregated.mat'])

ntrials = length(TestError{1}.rf);

for j = 1:3
    p = ps(j);
    
%     Classifiers = fieldnames(TestError{1,j});
    
    TimeMatrix = NaN(ntrials,length(ns),length(Classifiers));

    for i = 1:length(ns{j})
        n = ns{j}(i);

        for c = 1:length(Classifiers)
            cl = Classifiers{c};
            if ~strcmp(cl,'xgb')
                if ~isempty(TrainTime{i,j}.(cl))
                    OE = OOBError{i,j}.(cl);
                    OA = OOBAUC{i,j}.(cl);
                    ntrials = size(OE,1);
                    T = NaN(ntrials,1);
                    for trial = 1:ntrials
                        BI = hp_optimize(OE(trial,1:length(Params{i,j}.(cl).d)),OA(trial,1:length(Params{i,j}.(cl).d)));
                        T(trial) = TrainTime{i,j}.(cl)(trial,BI(end));
                    end
%                     TimeMatrix(:,i,c) = T;
                    TimeMatrix(:,i,c) = sum(TrainTime{i,j}.(cl),2);
%                     TimeMatrix(:,i,c) = TrainTime{i,j}.(cl)';
                else
                    TimeMatrix(:,i,c) = NaN;
                end
            else
                try
                    TimeMatrix(:,i,c) = ...
                        dlmread([rerfPath ...
                        sprintf('RandomerForest/Results/pami/Orthant/dat/Orthant_vary_n_trainTime_n%d_p%d.dat',n,p)]);
                catch ME
                    TimeMatrix(:,i,c) = NaN;
                end
            end
        end
    end

    ax((j-1)*3+k) = axes;
    hold on
    ymax = zeros(length(Classifiers),1);
    ymin = zeros(length(Classifiers),1);
    for c = 1:length(Classifiers)
        mn = mean(TimeMatrix(:,:,c));
        sem = std(TimeMatrix(:,:,c))/sqrt(ntrials);
        errorbar(ns{j},mn,sem,...
            'LineWidth',LineWidth,'Color',Colors.(Classifiers{c}))
        ymax(c) = mn(end)+2*sem(end);
        ymin(c) = mn(1)-2*sem(1);
    end

    ax((j-1)*3+k).LineWidth = LineWidth;
    ax((j-1)*3+k).FontUnits = 'inches';
    ax((j-1)*3+k).FontSize = FontSize;
    ax((j-1)*3+k).Units = 'inches';
    ax((j-1)*3+k).Position = [axLeft((j-1)*3+k) axBottom((j-1)*3+k) axWidth axHeight];
    ax((j-1)*3+k).Box = 'off';
    ax((j-1)*3+k).XLim = [10^(log10(min(ns{j}))-0.1) 10^(log10(max(ns{j}))+0.1)];
    ax((j-1)*3+k).XScale = 'log';
    ax((j-1)*3+k).XTick = ns{j};
    ax((j-1)*3+k).XTickLabel = cellstr(num2str(ns{j}'))';
    ax((j-1)*3+k).XTickLabelRotation = 45;
    ax((j-1)*3+k).YLim = [min(ymin) max(ymax)];
    
%     ax((j-1)*3+k).XScale = 'log';
%     ax((j-1)*3+k).FontSize = 16;
%     ax((j-1)*3+k).XLim = [10^(log10(min(ns{j}))-0.1) 10^(log10(max(ns{j}))+0.1)];
%     ax((j-1)*3+k).YLim = [min(TimeMatrix(:)) max(TimeMatrix(:))];
%     ax((j-1)*3+k).XTick = ns{j};
%     ax((j-1)*3+k).XTickLabel = cellstr(num2str(ns{j}'))';

    xlabel('n')
    ylabel('Train Time (s)')
    if j == 1
        title({'Orthant';sprintf('p = %d',p)})
    else
        title(sprintf('p = %d',p))
    end

end

save_fig(gcf,[rerfPath 'RandomerForest/Figures/pami/PAMI_fig3_time_synthetic_data'],{'fig','pdf','png'})