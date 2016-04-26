%% Plot bias, variance, and generalization error for Trunk simulations

close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

rng(1);

%% Plotting parameters

LineWidth.rf = 1;
LineWidth.rerf = 1;
LineWidth.rf_rot = 1;
LineWidth.frc2 = 1;
LineWidth.frc3 = 2;
LineWidth.frc4 = 3;
LineWidth.frc5 = 4;
LineWidth.frc6 = 5;
C = [0 1 1;0 1 0;1 0 1;1 0 0;0 0 0;1 .5 0];
Colors.rf = C(1,:);
Colors.rerf = C(2,:);
Colors.rf_rot = C(3,:);
Colors.frc2 = C(4,:);
Colors.frc3 = C(4,:);
Colors.frc4 = C(4,:);
Colors.frc5 = C(4,:);
Colors.frc6 = C(4,:);

axWidth = 1;
axHeight = 1;
FontSize = 0.1;
axLeft = repmat([FontSize*5,FontSize*10+axWidth,FontSize*15+axWidth*2],1,5);
axBottom = [repmat(FontSize*19+axHeight*4,1,3),repmat(FontSize*15+axHeight*3,1,3),...
    repmat(FontSize*11+axHeight*2,1,3),repmat(FontSize*7+axHeight,1,3),...
    repmat(FontSize*3,1,3)];
% axLeft = repmat([FontSize*4,FontSize*8+axWidth,FontSize*12+axWidth*2],1,3);
% axBottom = [repmat(FontSize*12+axHeight*2,1,3),...
%     repmat(FontSize*8+axHeight,1,3),repmat(FontSize*4,1,3)];
figWidth = axLeft(end) + axWidth + FontSize*3;
figHeight = axBottom(1) + axHeight + FontSize*3;

fig = figure;
fig.Units = 'inches';
fig.PaperUnits = 'inches';
fig.Position = [0 0 figWidth figHeight];
fig.PaperPosition = [0 0 figWidth figHeight];
fig.PaperSize = [figWidth figHeight];

%% Run simulation or load data

runSims = false;

if runSims
    run_Trunk_bias_variance
else
    load Trunk_bias_variance
end

load('Trunk_partitioned_data','dims')

Results.GE = GE;
Results.V = V;
Results.B = B;
MetricNames.GE = 'Generalization Error';
MetricNames.V = 'Variance';
MetricNames.B = 'Bias^2';

Classifiers = fieldnames(Params);
Metrics = fieldnames(Results);

%% Plot each metric vs mtry

for i = 1:length(dims)
    d = dims(i);
    for j = 1:length(Metrics)
        Metric = Metrics{j};
        SubplotIdx = (i-1)*length(Metrics)+j;
        ax(SubplotIdx) = subplot(length(dims),length(Metrics),SubplotIdx);
        hold on
        for k = 1:length(Classifiers)
            cl = Classifiers{k};
            if ~isempty(Results.(Metric).(cl){i})
                plot(Params.(cl).mtry{i},Results.(Metric).(cl){i},...
                    'LineWidth',LineWidth.(cl),'Color',Colors.(cl))
            end
        end
        title(sprintf('%s (d = %d)',Metric,d))
        if i == 1
            xlabel('mtry')
            ylabel(MetricNames.(Metric))
        end
    end
end

for i = 1:length(dims)*length(Metrics)
    ax(i).LineWidth = 2;
    ax(i).FontUnits = 'inches';
    ax(i).FontSize = FontSize;
    ax(i).Units = 'inches';
    ax(i).Position = [axLeft(i) axBottom(i) axWidth axHeight];
    ax(i).Box = 'off';
    if d > 5
        ax(i).XScale = 'log';
    end
end

save_fig(gcf,[rerfPath 'RandomerForest/Figures/Trunk_bias_variance'])