clear
close all
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);
if isempty(rerfPath)
    rerfPath = '~/';
end

load('purple2green')

Colors.rf = ColorMap(2,:);
Colors.rerf = 'k';
Colors.frc= ColorMap(10,:);
Colors.rr_rf = ColorMap(4,:);
LineWidth = 2;
MarkerSize = 12;
FontSize = .2;
axWidth = 2;
axHeight = 1.5;
axLeft = repmat([FontSize*4,FontSize*8+axWidth],1,3);
axBottom = [...
    (FontSize*14+axHeight*2)*ones(1,2),(FontSize*9+axHeight)*ones(1,2),...
    FontSize*4*ones(1,2)];
legWidth = axWidth;
legHeight = axHeight;
legLeft = axLeft(end) + axWidth + FontSize;
legBottom = axBottom(4);
figWidth = legLeft + legWidth;
figHeight = axBottom(1) + axHeight + FontSize*2.5;

fig = figure;
fig.Units = 'inches';
fig.PaperUnits = 'inches';
fig.Position = [0 0 figWidth figHeight];
fig.PaperPosition = [0 0 figWidth figHeight];
fig.PaperSize = [figWidth figHeight];

% Plot Gaussian 0 degrees

k = 1;

load([rerfPath 'RandomerForest/Results/2017.06.04/Gaussian/Gaussian_0_deg_2017_06_04.mat'])

ntrials = length(TestError{1}.rf);
ns = ns(1:end-1);
Classifiers = {'rf','rerf','frc','rr_rf',};

ErrorMatrix = NaN(ntrials,length(ns),length(Classifiers));
BiasMatrix = NaN(1,length(ns),length(Classifiers));
VarianceMatrix = NaN(1,length(ns),length(Classifiers));

p = 2;

for i = 1:length(ns)
    n = ns(i);

    for c = 1:length(Classifiers)
        cl = Classifiers{c};
        if ~isempty(TestError{i}.(cl))
            OE = OOBError{i}.(cl);
            OA = OOBAUC{i}.(cl);
            ntrials = size(OE,1);
            TE = NaN(ntrials,1);
            for trial = 1:ntrials
                TE(trial) = TestError{i}.(cl)(trial,BestIdx{i}.(cl)(trial));
            end
            ErrorMatrix(:,i,c) = TE;
            BiasMatrix(1,i,c) = Bias{i}.(cl);
            VarianceMatrix(1,i,c) = Variance{i}.(cl);
        else
            ErrorMatrix(:,i,c) = NaN;
            BiasMatrix(1,i,c) = NaN;
            VarianceMatrix(1,i,c) = NaN;
        end
    end
end

% bias

j = 1;
ax((j-1)*2+k) = axes;
hold on
ymax = max(BiasMatrix(:));
ymin = min(BiasMatrix(:));
for c = 1:length(Classifiers)
   plot(ns,BiasMatrix(1,:,c),...
        'LineWidth',LineWidth,'Color',Colors.(Classifiers{c}))
end

ax((j-1)*2+k).LineWidth = LineWidth;
ax((j-1)*2+k).FontUnits = 'inches';
ax((j-1)*2+k).FontSize = FontSize;
ax((j-1)*2+k).Units = 'inches';
ax((j-1)*2+k).Position = [axLeft((j-1)*2+k) axBottom((j-1)*2+k) axWidth axHeight];
ax((j-1)*2+k).Box = 'off';
ax((j-1)*2+k).XLim = [min(ns)-1 max(ns)+1];
ax((j-1)*2+k).XScale = 'log';
ax((j-1)*2+k).XTick = ns;
ax((j-1)*2+k).XTickLabel = cellstr(num2str(ns'))';
ax((j-1)*2+k).XTickLabelRotation = 0;
ax((j-1)*2+k).YLim = [ymin ymax];

xlabel('n')
ylabel('Bias')
title({'Binary 2D Gaussians';'Axis-aligned'})


% variance

j = 2;
ax((j-1)*2+k) = axes;
hold on
ymax = max(VarianceMatrix(:));
ymin = min(VarianceMatrix(:));
for c = 1:length(Classifiers)
   plot(ns,VarianceMatrix(1,:,c),...
        'LineWidth',LineWidth,'Color',Colors.(Classifiers{c}))
end

ax((j-1)*2+k).LineWidth = LineWidth;
ax((j-1)*2+k).FontUnits = 'inches';
ax((j-1)*2+k).FontSize = FontSize;
ax((j-1)*2+k).Units = 'inches';
ax((j-1)*2+k).Position = [axLeft((j-1)*2+k) axBottom((j-1)*2+k) axWidth axHeight];
ax((j-1)*2+k).Box = 'off';
ax((j-1)*2+k).XLim = [min(ns)-1 max(ns)+1];
ax((j-1)*2+k).XScale = 'log';
ax((j-1)*2+k).XTick = ns;
ax((j-1)*2+k).XTickLabel = cellstr(num2str(ns'))';
ax((j-1)*2+k).XTickLabelRotation = 0;
ax((j-1)*2+k).YLim = [ymin ymax];

xlabel('n')
ylabel('Variance')

% error rate

j = 3;
ax((j-1)*2+k) = axes;
hold on
ymax = zeros(length(Classifiers),1);
ymin = zeros(length(Classifiers),1);
for c = 1:length(Classifiers)
    mn = mean(ErrorMatrix(:,:,c));
    sem = std(ErrorMatrix(:,:,c))/sqrt(ntrials);
    errorbar(ns,mn,sem,...
        'LineWidth',LineWidth,'Color',Colors.(Classifiers{c}))
    ymin(c) = mn(end)-2*sem(end);
    ymax(c) = mn(1)+2*sem(1);
end

ax((j-1)*2+k).LineWidth = LineWidth;
ax((j-1)*2+k).FontUnits = 'inches';
ax((j-1)*2+k).FontSize = FontSize;
ax((j-1)*2+k).Units = 'inches';
ax((j-1)*2+k).Position = [axLeft((j-1)*2+k) axBottom((j-1)*2+k) axWidth axHeight];
ax((j-1)*2+k).Box = 'off';
ax((j-1)*2+k).XLim = [min(ns)-1 max(ns)+1];
ax((j-1)*2+k).XScale = 'log';
ax((j-1)*2+k).XTick = ns;
ax((j-1)*2+k).XTickLabel = cellstr(num2str(ns'))';
ax((j-1)*2+k).XTickLabelRotation = 0;
ax((j-1)*2+k).YLim = [min(ymin) max(ymax)];

xlabel('n')
ylabel('Error Rate')

% Plot Gaussian 45 degrees

k = 2;

load([rerfPath 'RandomerForest/Results/2017.06.04/Gaussian/Gaussian_45_deg_2017_06_04.mat'])

ntrials = length(TestError{1}.rf);
ns = ns(1:end-1);
Classifiers = {'rf','rerf','frc','rr_rf',};

ErrorMatrix = NaN(ntrials,length(ns),length(Classifiers));
BiasMatrix = NaN(1,length(ns),length(Classifiers));
VarianceMatrix = NaN(1,length(ns),length(Classifiers));

p = 2;

for i = 1:length(ns)
    n = ns(i);

    for c = 1:length(Classifiers)
        cl = Classifiers{c};
        if ~isempty(TestError{i}.(cl))
            OE = OOBError{i}.(cl);
            OA = OOBAUC{i}.(cl);
            ntrials = size(OE,1);
            TE = NaN(ntrials,1);
            for trial = 1:ntrials
                TE(trial) = TestError{i}.(cl)(trial,BestIdx{i}.(cl)(trial));
            end
            ErrorMatrix(:,i,c) = TE;
            BiasMatrix(1,i,c) = Bias{i}.(cl);
            VarianceMatrix(1,i,c) = Variance{i}.(cl);
        else
            ErrorMatrix(:,i,c) = NaN;
            BiasMatrix(1,i,c) = NaN;
            VarianceMatrix(1,i,c) = NaN;
        end
    end
end

% bias

j = 1;
ax((j-1)*2+k) = axes;
hold on
ymax = max(BiasMatrix(:));
ymin = min(BiasMatrix(:));
for c = 1:length(Classifiers)
   plot(ns,BiasMatrix(1,:,c),...
        'LineWidth',LineWidth,'Color',Colors.(Classifiers{c}))
end

ax((j-1)*2+k).LineWidth = LineWidth;
ax((j-1)*2+k).FontUnits = 'inches';
ax((j-1)*2+k).FontSize = FontSize;
ax((j-1)*2+k).Units = 'inches';
ax((j-1)*2+k).Position = [axLeft((j-1)*2+k) axBottom((j-1)*2+k) axWidth axHeight];
ax((j-1)*2+k).Box = 'off';
ax((j-1)*2+k).XLim = [min(ns)-1 max(ns)+1];
ax((j-1)*2+k).XScale = 'log';
ax((j-1)*2+k).XTick = ns;
ax((j-1)*2+k).XTickLabel = cellstr(num2str(ns'))';
ax((j-1)*2+k).XTickLabelRotation = 0;
ax((j-1)*2+k).YLim = [ymin ymax];

xlabel('n')
ylabel('Bias')
title({'Binary 2D Gaussians';'45° Oblique'})

% variance

j = 2;
ax((j-1)*2+k) = axes;
hold on
ymax = max(VarianceMatrix(:));
ymin = min(VarianceMatrix(:));
for c = 1:length(Classifiers)
   plot(ns,VarianceMatrix(1,:,c),...
        'LineWidth',LineWidth,'Color',Colors.(Classifiers{c}))
end

ax((j-1)*2+k).LineWidth = LineWidth;
ax((j-1)*2+k).FontUnits = 'inches';
ax((j-1)*2+k).FontSize = FontSize;
ax((j-1)*2+k).Units = 'inches';
ax((j-1)*2+k).Position = [axLeft((j-1)*2+k) axBottom((j-1)*2+k) axWidth axHeight];
ax((j-1)*2+k).Box = 'off';
ax((j-1)*2+k).XLim = [min(ns)-1 max(ns)+1];
ax((j-1)*2+k).XScale = 'log';
ax((j-1)*2+k).XTick = ns;
ax((j-1)*2+k).XTickLabel = cellstr(num2str(ns'))';
ax((j-1)*2+k).XTickLabelRotation = 0;
ax((j-1)*2+k).YLim = [ymin ymax];

xlabel('n')
ylabel('Variance')

lh = legend('RF','RerF','F-RC','RR-RF');
lh.Box = 'off';
lh.Units = 'inches';
lh.Position = [legLeft legBottom legWidth legHeight];

% error rate

j = 3;
ax((j-1)*2+k) = axes;
hold on
ymax = zeros(length(Classifiers),1);
ymin = zeros(length(Classifiers),1);
for c = 1:length(Classifiers)
    mn = mean(ErrorMatrix(:,:,c));
    sem = std(ErrorMatrix(:,:,c))/sqrt(ntrials);
    errorbar(ns,mn,sem,...
        'LineWidth',LineWidth,'Color',Colors.(Classifiers{c}))
    ymin(c) = mn(end)-2*sem(end);
    ymax(c) = mn(1)+2*sem(1);
end

ax((j-1)*2+k).LineWidth = LineWidth;
ax((j-1)*2+k).FontUnits = 'inches';
ax((j-1)*2+k).FontSize = FontSize;
ax((j-1)*2+k).Units = 'inches';
ax((j-1)*2+k).Position = [axLeft((j-1)*2+k) axBottom((j-1)*2+k) axWidth axHeight];
ax((j-1)*2+k).Box = 'off';
ax((j-1)*2+k).XLim = [min(ns)-1 max(ns)+1];
ax((j-1)*2+k).XScale = 'log';
ax((j-1)*2+k).XTick = ns;
ax((j-1)*2+k).XTickLabel = cellstr(num2str(ns'))';
ax((j-1)*2+k).XTickLabelRotation = 0;
ax((j-1)*2+k).YLim = [min(ymin) max(ymax)];

xlabel('n')
ylabel('Error Rate')

save_fig(gcf,[rerfPath 'RandomerForest/Figures/2017.06.04/Gaussians/Gaussian_bias_variance'],{'fig','pdf','png'})