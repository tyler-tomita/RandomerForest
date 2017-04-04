clear
close all
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

load('purple2green')

Colors.rf = ColorMap(2,:);
Colors.rerf= ColorMap(10,:);
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

load([rerfPath 'RandomerForest/Results/pami/Sparse_parity/mat/Sparse_parity_vary_n.mat'])

ntrials = length(TestError{1}.rf);
Classifiers = [fieldnames(TestError{1});'xgb'];

ErrorMatrix = NaN(ntrials,length(ps),length(Classifiers));
BiasMatrix = NaN(1,length(ps),length(Classifiers));
VarianceMatrix = NaN(1,length(ps),length(Classifiers));

for j = 1:3
    p = ps(j);
    

    for i = 1:length(ns{j})
        n = ns{j}(i);
        
        if n == 1000

            for c = 1:length(Classifiers)
                cl = Classifiers{c};
                if ~strcmp(cl,'xgb')
                    if ~isempty(TestError{i,j}.(cl))
                        ErrorMatrix(:,j,c) = TestError{i,j}.(cl);
                        BiasMatrix(1,j,c) = mean(Bias{i,j}.(cl)(BestIdx{i,j}.(cl)))';
                        VarianceMatrix(1,j,c) = mean(Variance{i,j}.(cl)(BestIdx{i,j}.(cl)))';
                    else
                        ErrorMatrix(:,j,c) = NaN;
                        BiasMatrix(1,j,c) = NaN;
                        VarianceMatrix(1,j,c) = NaN;
                    end
                else
                    try
                        ErrorMatrix(:,j,c) = ...
                            dlmread([rerfPath ...
                            sprintf('RandomerForest/Results/pami/Sparse_parity/dat/Sparse_parity_vary_n_testError_n%d_p%d.dat',n,p)]);
                        BiasMatrix(1,j,c) = ...
                            dlmread([rerfPath ...
                            sprintf('RandomerForest/Results/pami/Sparse_parity/dat/Sparse_parity_vary_n_bias_n%d_p%d.dat',n,p)]);
                        VarianceMatrix(1,j,c) = ...
                            dlmread([rerfPath ...
                            sprintf('RandomerForest/Results/pami/Sparse_parity/dat/Sparse_parity_vary_n_variance_n%d_p%d.dat',n,p)]);
                    catch
                        ErrorMatrix(:,j,c) = NaN;
                        BiasMatrix(1,j,c) = NaN;
                        VarianceMatrix(1,j,c) = NaN;
                    end
                end
            end
            
        end
    end
end

% bias

j = 1;
ax((j-1)*3+k) = axes;
hold on
ymax = max(BiasMatrix(:));
ymin = min(BiasMatrix(:));
for c = 1:length(Classifiers)
    plot(ps,BiasMatrix(1,:,c),...
        'LineWidth',LineWidth,'Color',Colors.(Classifiers{c}))
end

ax((j-1)*3+k).LineWidth = LineWidth;
ax((j-1)*3+k).FontUnits = 'inches';
ax((j-1)*3+k).FontSize = FontSize;
ax((j-1)*3+k).Units = 'inches';
ax((j-1)*3+k).Position = [axLeft((j-1)*3+k) axBottom((j-1)*3+k) axWidth axHeight];
ax((j-1)*3+k).Box = 'off';
ax((j-1)*3+k).XLim = [min(ps)-1 max(ps)+1];
% ax((j-1)*3+k).XScale = 'log';
ax((j-1)*3+k).XTick = ps;
ax((j-1)*3+k).XTickLabel = cellstr(num2str(ps'))';
ax((j-1)*3+k).XTickLabelRotation = 0;
ax((j-1)*3+k).YLim = [ymin ymax];

xlabel('p')
ylabel('Bias')
if j == 1
    title({'Sparse Parity';sprintf('n = %d',1000)})
end

% variance

j = 2;
ax((j-1)*3+k) = axes;
hold on
ymax = max(VarianceMatrix(:));
ymin = min(VarianceMatrix(:));
for c = 1:length(Classifiers)
    plot(ps,VarianceMatrix(1,:,c),...
        'LineWidth',LineWidth,'Color',Colors.(Classifiers{c}))
end

ax((j-1)*3+k).LineWidth = LineWidth;
ax((j-1)*3+k).FontUnits = 'inches';
ax((j-1)*3+k).FontSize = FontSize;
ax((j-1)*3+k).Units = 'inches';
ax((j-1)*3+k).Position = [axLeft((j-1)*3+k) axBottom((j-1)*3+k) axWidth axHeight];
ax((j-1)*3+k).Box = 'off';
ax((j-1)*3+k).XLim = [min(ps)-1 max(ps)+1];
% ax((j-1)*3+k).XScale = 'log';
ax((j-1)*3+k).XTick = ps;
ax((j-1)*3+k).XTickLabel = cellstr(num2str(ps'))';
ax((j-1)*3+k).XTickLabelRotation = 0;
ax((j-1)*3+k).YLim = [ymin ymax];

xlabel('p')
ylabel('Variance')
if j == 1
    title({'Sparse Parity';sprintf('n = %d',1000)})
end

% error rate

j = 3;
ax((j-1)*3+k) = axes;
hold on
ymax = zeros(length(Classifiers),1);
ymin = zeros(length(Classifiers),1);
for c = 1:length(Classifiers)
    mn = mean(ErrorMatrix(:,:,c));
    sem = std(ErrorMatrix(:,:,c))/sqrt(ntrials);
    errorbar(ps,mn,sem,...
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
ax((j-1)*3+k).XLim = [min(ps)-1 max(ps)+1];
% ax((j-1)*3+k).XScale = 'log';
ax((j-1)*3+k).XTick = ps;
ax((j-1)*3+k).XTickLabel = cellstr(num2str(ps'))';
ax((j-1)*3+k).XTickLabelRotation = 0;
ax((j-1)*3+k).YLim = [min(ymin) max(ymax)];

xlabel('p')
ylabel('Error Rate')
if j == 1
    title({'Sparse Parity';sprintf('n = %d',1000)})
end

% Plot Trunk

k = 2;

load([rerfPath 'RandomerForest/Results/pami/Trunk/mat/Trunk_vary_n.mat'])

ntrials = length(TestError{1}.rf);
Classifiers = [fieldnames(TestError{1});'xgb'];

ErrorMatrix = NaN(ntrials,length(ps),length(Classifiers));
BiasMatrix = NaN(1,length(ps),length(Classifiers));
VarianceMatrix = NaN(1,length(ps),length(Classifiers));

for j = 1:3
    p = ps(j);
    

    for i = 1:length(ns{j})
        n = ns{j}(i);
        
        if n == 100

            for c = 1:length(Classifiers)
                cl = Classifiers{c};
                if ~strcmp(cl,'xgb')
                    if ~isempty(TestError{i,j}.(cl))
                        ErrorMatrix(:,j,c) = TestError{i,j}.(cl);
                        BiasMatrix(1,j,c) = mean(Bias{i,j}.(cl)(BestIdx{i,j}.(cl)))';
                        VarianceMatrix(1,j,c) = mean(Variance{i,j}.(cl)(BestIdx{i,j}.(cl)))';
                    else
                        ErrorMatrix(:,j,c) = NaN;
                        BiasMatrix(1,j,c) = NaN;
                        VarianceMatrix(1,j,c) = NaN;
                    end
                else
                    try
                        ErrorMatrix(:,j,c) = ...
                            dlmread([rerfPath ...
                            sprintf('RandomerForest/Results/pami/Trunk/dat/Trunk_vary_n_testError_n%d_p%d.dat',n,p)]);
                        BiasMatrix(1,j,c) = ...
                            dlmread([rerfPath ...
                            sprintf('RandomerForest/Results/pami/Trunk/dat/Trunk_vary_n_bias_n%d_p%d.dat',n,p)]);
                        VarianceMatrix(1,j,c) = ...
                            dlmread([rerfPath ...
                            sprintf('RandomerForest/Results/pami/Trunk/dat/Trunk_vary_n_variance_n%d_p%d.dat',n,p)]);
                    catch
                        ErrorMatrix(:,j,c) = NaN;
                        BiasMatrix(1,j,c) = NaN;
                        VarianceMatrix(1,j,c) = NaN;
                    end
                end
            end
            
        end
    end
end

% bias

j = 1;
ax((j-1)*3+k) = axes;
hold on
ymax = max(BiasMatrix(:));
ymin = min(BiasMatrix(:));
for c = 1:length(Classifiers)
    plot(ps,BiasMatrix(1,:,c),...
        'LineWidth',LineWidth,'Color',Colors.(Classifiers{c}))
end

ax((j-1)*3+k).LineWidth = LineWidth;
ax((j-1)*3+k).FontUnits = 'inches';
ax((j-1)*3+k).FontSize = FontSize;
ax((j-1)*3+k).Units = 'inches';
ax((j-1)*3+k).Position = [axLeft((j-1)*3+k) axBottom((j-1)*3+k) axWidth axHeight];
ax((j-1)*3+k).Box = 'off';
ax((j-1)*3+k).XLim = [min(ps)-1 max(ps)+1];
ax((j-1)*3+k).XScale = 'log';
ax((j-1)*3+k).XTick = ps;
ax((j-1)*3+k).XTickLabel = cellstr(num2str(ps'))';
ax((j-1)*3+k).XTickLabelRotation = 0;
ax((j-1)*3+k).YLim = [ymin ymax];

xlabel('p')
ylabel('Bias')
if j == 1
    title({'Trunk';sprintf('n = %d',1000)})
end

% variance

j = 2;
ax((j-1)*3+k) = axes;
hold on
ymax = max(VarianceMatrix(:));
ymin = min(VarianceMatrix(:));
for c = 1:length(Classifiers)
    plot(ps,VarianceMatrix(1,:,c),...
        'LineWidth',LineWidth,'Color',Colors.(Classifiers{c}))
end

ax((j-1)*3+k).LineWidth = LineWidth;
ax((j-1)*3+k).FontUnits = 'inches';
ax((j-1)*3+k).FontSize = FontSize;
ax((j-1)*3+k).Units = 'inches';
ax((j-1)*3+k).Position = [axLeft((j-1)*3+k) axBottom((j-1)*3+k) axWidth axHeight];
ax((j-1)*3+k).Box = 'off';
ax((j-1)*3+k).XLim = [min(ps)-1 max(ps)+1];
ax((j-1)*3+k).XScale = 'log';
ax((j-1)*3+k).XTick = ps;
ax((j-1)*3+k).XTickLabel = cellstr(num2str(ps'))';
ax((j-1)*3+k).XTickLabelRotation = 0;
ax((j-1)*3+k).YLim = [ymin ymax];

xlabel('p')
ylabel('Variance')
if j == 1
    title({'Trunk';sprintf('n = %d',1000)})
end

lh = legend('RF','RerF','RR-RF','XGBoost');
lh.Box = 'off';
lh.Units = 'inches';
lh.Position = [legLeft legBottom legWidth legHeight];

% error rate

j = 3;
ax((j-1)*3+k) = axes;
hold on
ymax = zeros(length(Classifiers),1);
ymin = zeros(length(Classifiers),1);
for c = 1:length(Classifiers)
    mn = mean(ErrorMatrix(:,:,c));
    sem = std(ErrorMatrix(:,:,c))/sqrt(ntrials);
    errorbar(ps,mn,sem,...
        'LineWidth',LineWidth,'Color',Colors.(Classifiers{c}))
    ymax(c) = max(mn)+2*max(sem);
    ymin(c) = min(mn)-2*max(sem);
end

ax((j-1)*3+k).LineWidth = LineWidth;
ax((j-1)*3+k).FontUnits = 'inches';
ax((j-1)*3+k).FontSize = FontSize;
ax((j-1)*3+k).Units = 'inches';
ax((j-1)*3+k).Position = [axLeft((j-1)*3+k) axBottom((j-1)*3+k) axWidth axHeight];
ax((j-1)*3+k).Box = 'off';
ax((j-1)*3+k).XLim = [min(ps)-1 max(ps)+1];
ax((j-1)*3+k).XScale = 'log';
ax((j-1)*3+k).XTick = ps;
ax((j-1)*3+k).XTickLabel = cellstr(num2str(ps'))';
ax((j-1)*3+k).XTickLabelRotation = 0;
ax((j-1)*3+k).YLim = [min(ymin) max(ymax)];

xlabel('p')
ylabel('Error Rate')
if j == 1
    title({'Trunk';sprintf('n = %d',1000)})
end

% Plot Orthant

k = 3;

load([rerfPath 'RandomerForest/Results/pami/Orthant/mat/Orthant_vary_n.mat'])

ntrials = length(TestError{1}.rf);
Classifiers = [fieldnames(TestError{1});'xgb'];

ErrorMatrix = NaN(ntrials,length(ps),length(Classifiers));
BiasMatrix = NaN(1,length(ps),length(Classifiers));
VarianceMatrix = NaN(1,length(ps),length(Classifiers));

for j = 1:3
    p = ps(j);
    

    for i = 1:length(ns{j})
        n = ns{j}(i);
        
        if n == 400

            for c = 1:length(Classifiers)
                cl = Classifiers{c};
                if ~strcmp(cl,'xgb')
                    if ~isempty(TestError{i,j}.(cl))
                        ErrorMatrix(:,j,c) = TestError{i,j}.(cl);
                        BiasMatrix(1,j,c) = mean(Bias{i,j}.(cl)(BestIdx{i,j}.(cl)))';
                        VarianceMatrix(1,j,c) = mean(Variance{i,j}.(cl)(BestIdx{i,j}.(cl)))';
                    else
                        ErrorMatrix(:,j,c) = NaN;
                        BiasMatrix(1,j,c) = NaN;
                        VarianceMatrix(1,j,c) = NaN;
                    end
                else
                    try
                        ErrorMatrix(:,j,c) = ...
                            dlmread([rerfPath ...
                            sprintf('RandomerForest/Results/pami/Orthant/dat/Orthant_vary_n_testError_n%d_p%d.dat',n,p)]);
                        BiasMatrix(1,j,c) = ...
                            dlmread([rerfPath ...
                            sprintf('RandomerForest/Results/pami/Orthant/dat/Orthant_vary_n_bias_n%d_p%d.dat',n,p)]);
                        VarianceMatrix(1,j,c) = ...
                            dlmread([rerfPath ...
                            sprintf('RandomerForest/Results/pami/Orthant/dat/Orthant_vary_n_variance_n%d_p%d.dat',n,p)]);
                    catch
                        ErrorMatrix(:,j,c) = NaN;
                        BiasMatrix(1,j,c) = NaN;
                        VarianceMatrix(1,j,c) = NaN;
                    end
                end
            end
            
        end
    end
end

% bias

j = 1;
ax((j-1)*3+k) = axes;
hold on
ymax = max(BiasMatrix(:));
ymin = min(BiasMatrix(:));
for c = 1:length(Classifiers)
    plot(ps,BiasMatrix(1,:,c),...
        'LineWidth',LineWidth,'Color',Colors.(Classifiers{c}))
end

ax((j-1)*3+k).LineWidth = LineWidth;
ax((j-1)*3+k).FontUnits = 'inches';
ax((j-1)*3+k).FontSize = FontSize;
ax((j-1)*3+k).Units = 'inches';
ax((j-1)*3+k).Position = [axLeft((j-1)*3+k) axBottom((j-1)*3+k) axWidth axHeight];
ax((j-1)*3+k).Box = 'off';
ax((j-1)*3+k).XLim = [min(ps)-1 max(ps)+1];
% ax((j-1)*3+k).XScale = 'log';
ax((j-1)*3+k).XTick = ps;
ax((j-1)*3+k).XTickLabel = cellstr(num2str(ps'))';
ax((j-1)*3+k).XTickLabelRotation = 0;
ax((j-1)*3+k).YLim = [ymin ymax];

xlabel('p')
ylabel('Bias')
if j == 1
    title({'Orthant';sprintf('n = %d',1000)})
end

% variance

j = 2;
ax((j-1)*3+k) = axes;
hold on
ymax = max(VarianceMatrix(:));
ymin = min(VarianceMatrix(:));
for c = 1:length(Classifiers)
    plot(ps,VarianceMatrix(1,:,c),...
        'LineWidth',LineWidth,'Color',Colors.(Classifiers{c}))
end

ax((j-1)*3+k).LineWidth = LineWidth;
ax((j-1)*3+k).FontUnits = 'inches';
ax((j-1)*3+k).FontSize = FontSize;
ax((j-1)*3+k).Units = 'inches';
ax((j-1)*3+k).Position = [axLeft((j-1)*3+k) axBottom((j-1)*3+k) axWidth axHeight];
ax((j-1)*3+k).Box = 'off';
ax((j-1)*3+k).XLim = [min(ps)-1 max(ps)+1];
% ax((j-1)*3+k).XScale = 'log';
ax((j-1)*3+k).XTick = ps;
ax((j-1)*3+k).XTickLabel = cellstr(num2str(ps'))';
ax((j-1)*3+k).XTickLabelRotation = 0;
ax((j-1)*3+k).YLim = [ymin ymax];

xlabel('p')
ylabel('Variance')
if j == 1
    title({'Orthant';sprintf('n = %d',1000)})
end

% error rate

j = 3;
ax((j-1)*3+k) = axes;
hold on
ymax = zeros(length(Classifiers),1);
ymin = zeros(length(Classifiers),1);
for c = 1:length(Classifiers)
    mn = mean(ErrorMatrix(:,:,c));
    sem = std(ErrorMatrix(:,:,c))/sqrt(ntrials);
    errorbar(ps,mn,sem,...
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
ax((j-1)*3+k).XLim = [min(ps)-1 max(ps)+1];
% ax((j-1)*3+k).XScale = 'log';
ax((j-1)*3+k).XTick = ps;
ax((j-1)*3+k).XTickLabel = cellstr(num2str(ps'))';
ax((j-1)*3+k).XTickLabelRotation = 0;
ax((j-1)*3+k).YLim = [min(ymin) max(ymax)];

xlabel('p')
ylabel('Error Rate')
if j == 1
    title({'Orthant';sprintf('n = %d',1000)})
end

save_fig(gcf,[rerfPath 'RandomerForest/Figures/pami/PAMI_fig6_bias_variance_synthetic_data'],{'fig','pdf','png'})