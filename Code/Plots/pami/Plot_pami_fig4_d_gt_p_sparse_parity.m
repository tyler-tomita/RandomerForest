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
axWidth = 1.5;
axHeight = 1.5;
axLeft = [FontSize*4,FontSize*8+axWidth];
axBottom = FontSize*4*ones(1,2);
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

load([rerfPath 'RandomerForest/Results/pami/Sparse_parity/mat/Sparse_parity_vary_n.mat'])

ntrials = length(TestError{5}.rf);

p = ps(2);
n = ns{2}(2);

Classifiers = fieldnames(TestError{5});

k = 1;
ax(k) = axes;
hold on
ymax = zeros(length(Classifiers),1);
ymin = zeros(length(Classifiers),1);
for c = 1:length(Classifiers)
    cl = Classifiers{c};
    mn = mean(OOBError{5}.(cl)(:,:,end));
    sem = std(OOBError{5}.(cl)(:,:,end))/sqrt(ntrials);
    errorbar(Params{5}.(cl).d,mn,sem,...
        'LineWidth',LineWidth,'Color',Colors.(Classifiers{c}))
    xmax(c) = max(Params{5}.(cl).d);
    xmin(c) = min(Params{5}.(cl).d);
    ymax(c) = mn(1)+2*sem(1);
    ymin(c) = mn(end)-2*sem(end);
end

ax(k).LineWidth = LineWidth;
ax(k).FontUnits = 'inches';
ax(k).FontSize = FontSize;
ax(k).Units = 'inches';
ax(k).Position = [axLeft(k) axBottom(k) axWidth axHeight];
ax(k).Box = 'off';
ax(k).XLim = [min(xmin) max(xmax)];
ax(k).XScale = 'log';
ax(k).XTick = Params{5}.rerf.d;
ax(k).XTick([2,3]) = [];
ax(k).XTickLabel = cellstr(num2str(ax(k).XTick'))';
ax(k).XTickLabelRotation = 0;
ax(k).YLim = [min(ymin) max(ymax)];

xlabel('d')
ylabel('OOB Error')

k = 2;
ax(k) = axes;
hold on
ymax = zeros(length(Classifiers),1);
ymin = zeros(length(Classifiers),1);
for c = 1:length(Classifiers)
    cl = Classifiers{c};
    mn = mean(TrainTime{5}.(cl));
    sem = std(TrainTime{5}.(cl))/sqrt(ntrials);
    errorbar(Params{5}.(cl).d,mn,sem,...
        'LineWidth',LineWidth,'Color',Colors.(Classifiers{c}))
    xmax(c) = max(Params{5}.(cl).d);
    xmin(c) = min(Params{5}.(cl).d);
    ymax(c) = mn(end)+2*sem(end);
    ymin(c) = mn(1)-2*sem(1);
end

ax(k).LineWidth = LineWidth;
ax(k).FontUnits = 'inches';
ax(k).FontSize = FontSize;
ax(k).Units = 'inches';
ax(k).Position = [axLeft(k) axBottom(k) axWidth axHeight];
ax(k).Box = 'off';
ax(k).XLim = [min(xmin) max(xmax)];
ax(k).XScale = 'log';
ax(k).XTick = Params{5}.rerf.d;
ax(k).XTick([2,3]) = [];
ax(k).XTickLabel = cellstr(num2str(ax(k).XTick'))';
ax(k).XTickLabelRotation = 0;
ax(k).YLim = [min(ymin) max(ymax)];

xlabel('d')
ylabel('Train Time (s)')

lh = legend('RF','RerF','RR-RF');
lh.Box = 'off';
lh.Units = 'inches';
lh.Position = [legLeft legBottom legWidth legHeight];

save_fig(gcf,[rerfPath 'RandomerForest/Figures/pami/PAMI_fig4_d_gt_p_sparse_parity'],{'fig','pdf','png'})