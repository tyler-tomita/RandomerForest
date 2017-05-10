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
MarkerSize = 10;
FontSize = .2;
axWidth = 1.5;
axHeight = 1.5;
axLeft = [FontSize*4,FontSize*8+axWidth,FontSize*12+axWidth*2];
axBottom = FontSize*4*ones(1,3);
legWidth = axWidth;
legHeight = axHeight;
legLeft = axLeft(end) + axWidth + FontSize;
legBottom = axBottom(end);
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

Classifiers = {'rf','rerf','frc','rr_rf'};

% Markers = {'o','x','+'};

ax(k) = axes;
hold on
maxDepth = 0;
for j = 2
    p = ps(j);
    
    %Classifiers = [fieldnames(TestError{1,j})];
    
    for i = 1:length(ns{j})
        n = ns{j}(i);
        
        if n == 1000

            for c = 1:length(Classifiers)
                cl = Classifiers{c};
                if ~strcmp(cl,'xgb')
                    if ~isempty(TestError{i,j}.(cl))
                        OE = OOBError{i,j}.(cl);
                        OA = OOBAUC{i,j}.(cl);
                        ntrials = size(OE,1);
                        E = NaN(1,ntrials);
                        BI = NaN(ntrials,1);
                        for trial = 1:ntrials
                            B = hp_optimize(OE(trial,1:length(Params{i,j}.(cl).d)),OA(trial,1:length(Params{i,j}.(cl).d)));
                            BI(trial) = B(end);
                            E(trial) = TestError{i,j}.(cl)(trial,BI(trial));
                        end
%                         E = TestError{i,j}.(cl);
                        D = reshape(mean(Depth{i,j}.(cl),2),ntrials,size(Depth{i,j}.(cl),3));
                        D = D(sub2ind(size(D),(1:ntrials)',BI));
                    else
                        E = NaN;
                        D = NaN;
                    end
                else
                    try
                        E = dlmread([rerfPath ...
                            sprintf('RandomerForest/Results/pami/Sparse_parity/dat/Sparse_parity_vary_n_testError_n%d_p%d.dat',n,p)]);
                        D = dlmread([rerfPath ...
                            sprintf('RandomerForest/Results/pami/Sparse_parity/dat/Sparse_parity_vary_n_depth_n%d_p%d.dat',n,p)]);
                    catch ME
                        E = NaN;
                        D = NaN;
                    end
                end
                plot(D,E,'.','MarkerEdgeColor',Colors.(cl),...
                    'MarkerSize',MarkerSize,'LineWidth',0.5)
                maxDepth = max([maxDepth,max(D)]);
            end
            
        end
    end

%     ax(k) = axes;
%     hold on
%     ymax = zeros(length(Classifiers),1);
%     ymin = zeros(length(Classifiers),1);
%     for c = 1:length(Classifiers)
%         mn = mean(ErrorMatrix(:,:,c));
%         sem = std(ErrorMatrix(:,:,c))/sqrt(ntrials);
%         errorbar(ns{j},mn,sem,...
%             'LineWidth',LineWidth,'Color',Colors.(Classifiers{c}))
%         ymax(c) = mn(1)+2*sem(1);
%         ymin(c) = mn(end)-2*sem(end);
%     end

end

xlabel('Tree Depth')
ylabel('Error Rate')
title({'Sparse Parity';sprintf('n = %d',1000)})

ax(k).LineWidth = LineWidth;
ax(k).FontUnits = 'inches';
ax(k).FontSize = FontSize;
ax(k).Units = 'inches';
ax(k).Position = [axLeft(k) axBottom(k) axWidth axHeight];
ax(k).Box = 'off';
ax(k).XLim = [10,maxDepth];
% ax(k).XScale = 'log';
% ax(k).XTick = ns{j};
% ax(k).XTickLabel = cellstr(num2str(ns{j}'))';
% ax(k).XTickLabelRotation = 45;
% ax(k).YLim = [min(ErrorMatrix(:)) max(ErrorMatrix(:))];
% ax(k).YScale = 'log';

% Plot Trunk

k = 2;

% load([rerfPath 'RandomerForest/Results/pami/Trunk/mat/Trunk_vary_n.mat'])
load([rerfPath 'RandomerForest/Results/2017.04.12/Trunk/Trunk_raw_vary_n_aggregated.mat'])

ntrials = length(TestError{1}.rf);

Markers = {'o','x','+'};

ax(k) = axes;
hold on
maxDepth = 0;
for j = 1
    p = ps(j);
    
    %Classifiers = [fieldnames(TestError{1,j})];
    
    for i = 1:length(ns{j})
        n = ns{j}(i);
        
        if n == 100

            for c = 1:length(Classifiers)
                cl = Classifiers{c};
                if ~strcmp(cl,'xgb')
                    if ~isempty(TestError{i,j}.(cl))
                        OE = OOBError{i,j}.(cl);
                        OA = OOBAUC{i,j}.(cl);
                        ntrials = size(OE,1);
                        E = NaN(1,ntrials);
                        BI = NaN(ntrials,1);
                        for trial = 1:ntrials
                            B = hp_optimize(OE(trial,1:length(Params{i,j}.(cl).d)),OA(trial,1:length(Params{i,j}.(cl).d)));
                            BI(trial) = B(end);
                            E(trial) = TestError{i,j}.(cl)(trial,BI(trial));
                        end
%                         E = TestError{i,j}.(cl);
                        D = reshape(mean(Depth{i,j}.(cl),2),ntrials,size(Depth{i,j}.(cl),3));
                        D = D(sub2ind(size(D),(1:ntrials)',BI));
                    else
                        E = NaN;
                        D = NaN;
                    end
                else
                    try
                        E = dlmread([rerfPath ...
                            sprintf('RandomerForest/Results/pami/Trunk/dat/Trunk_vary_n_testError_n%d_p%d.dat',n,p)]);
                        D = dlmread([rerfPath ...
                            sprintf('RandomerForest/Results/pami/Trunk/dat/Trunk_vary_n_depth_n%d_p%d.dat',n,p)]);
                    catch ME
                        E = NaN;
                        D = NaN;
                    end
                end
                plot(D,E,'.','MarkerEdgeColor',Colors.(cl),...
                    'MarkerSize',MarkerSize,'LineWidth',0.5)
                maxDepth = max([maxDepth,max(D)]);
            end
            
        end
    end

%     ax(k) = axes;
%     hold on
%     ymax = zeros(length(Classifiers),1);
%     ymin = zeros(length(Classifiers),1);
%     for c = 1:length(Classifiers)
%         mn = mean(ErrorMatrix(:,:,c));
%         sem = std(ErrorMatrix(:,:,c))/sqrt(ntrials);
%         errorbar(ns{j},mn,sem,...
%             'LineWidth',LineWidth,'Color',Colors.(Classifiers{c}))
%         ymax(c) = mn(1)+2*sem(1);
%         ymin(c) = mn(end)-2*sem(end);
%     end

end

xlabel('Tree Depth')
ylabel('Error Rate')
title({'Trunk';sprintf('n = %d',100)})

ax(k).LineWidth = LineWidth;
ax(k).FontUnits = 'inches';
ax(k).FontSize = FontSize;
ax(k).Units = 'inches';
ax(k).Position = [axLeft(k) axBottom(k) axWidth axHeight];
ax(k).Box = 'off';
ax(k).XLim = [0,10];
% ax(k).XScale = 'log';
% ax(k).XTick = ns{j};
% ax(k).XTickLabel = cellstr(num2str(ns{j}'))';
% ax(k).XTickLabelRotation = 45;
% ax(k).YLim = [min(ErrorMatrix(:)) max(ErrorMatrix(:))];
% ax(k).YScale = 'log';

lh = legend('RF','RerF','F-RC','RR-RF');
lh.Box = 'off';
lh.Units = 'inches';
lh.Position = [legLeft legBottom legWidth legHeight];

% Plot Orthant

k = 3;

% load([rerfPath 'RandomerForest/Results/pami/Orthant/mat/Orthant_vary_n.mat'])
load([rerfPath 'RandomerForest/Results/2017.04.13/Orthant/Orthant_raw_vary_n_aggregated.mat'])

ntrials = length(TestError{1}.rf);

Markers = {'o','x','+'};

ax(k) = axes;
hold on
maxDepth = 0;
for j = 1
    p = ps(j);
    
%     Classifiers = fieldnames(TestError{1,j});
%     Classifiers = [fieldnames(TestError{1,j});'xgb'];
    
    for i = 1:length(ns{j})
        n = ns{j}(i);
        
        if n == 400

            for c = 1:length(Classifiers)
                cl = Classifiers{c};
                if ~strcmp(cl,'xgb')
                    if ~isempty(TestError{i,j}.(cl))
                        OE = OOBError{i,j}.(cl);
                        OA = OOBAUC{i,j}.(cl);
                        ntrials = size(OE,1);
                        E = NaN(1,ntrials);
                        BI = NaN(ntrials,1);
                        for trial = 1:ntrials
                            B = hp_optimize(OE(trial,1:length(Params{i,j}.(cl).d)),OA(trial,1:length(Params{i,j}.(cl).d)));
                            BI(trial) = B(end);
                            E(trial) = TestError{i,j}.(cl)(trial,BI(trial));
                        end
%                         E = TestError{i,j}.(cl);
                        D = reshape(mean(Depth{i,j}.(cl),2),ntrials,size(Depth{i,j}.(cl),3));
                        D = D(sub2ind(size(D),(1:ntrials)',BI));
                    else
                        E = NaN;
                        D = NaN;
                    end
                else
                    try
                        E = dlmread([rerfPath ...
                            sprintf('RandomerForest/Results/pami/Orthant/dat/Orthant_vary_n_testError_n%d_p%d.dat',n,p)]);
                        D = dlmread([rerfPath ...
                            sprintf('RandomerForest/Results/pami/Orthant/dat/Orthant_vary_n_depth_n%d_p%d.dat',n,p)]);
                    catch ME
                        E = NaN;
                        D = NaN;
                    end
                end
                plot(D,E,'.','MarkerEdgeColor',Colors.(cl),...
                    'MarkerSize',MarkerSize,'LineWidth',0.5)
                maxDepth = max([maxDepth,max(D)]);            
            end
        end
    end

%     ax(k) = axes;
%     hold on
%     ymax = zeros(length(Classifiers),1);
%     ymin = zeros(length(Classifiers),1);
%     for c = 1:length(Classifiers)
%         mn = mean(ErrorMatrix(:,:,c));
%         sem = std(ErrorMatrix(:,:,c))/sqrt(ntrials);
%         errorbar(ns{j},mn,sem,...
%             'LineWidth',LineWidth,'Color',Colors.(Classifiers{c}))
%         ymax(c) = mn(1)+2*sem(1);
%         ymin(c) = mn(end)-2*sem(end);
%     end

end

xlabel('Tree Depth')
ylabel('Error Rate')
title({'Orthant';sprintf('n = %d',400)})

ax(k).LineWidth = LineWidth;
ax(k).FontUnits = 'inches';
ax(k).FontSize = FontSize;
ax(k).Units = 'inches';
ax(k).Position = [axLeft(k) axBottom(k) axWidth axHeight];
ax(k).Box = 'off';
ax(k).XLim = [0,20];
% ax(k).XScale = 'log';
% ax(k).XTick = ns{j};
% ax(k).XTickLabel = cellstr(num2str(ns{j}'))';
% ax(k).XTickLabelRotation = 45;
% ax(k).YLim = [min(ErrorMatrix(:)) max(ErrorMatrix(:))];
ax(k).YScale = 'log';

save_fig(gcf,[rerfPath 'RandomerForest/Figures/pami/PAMI_fig5_tree_size_synthetic_data'],{'fig','pdf','png'})