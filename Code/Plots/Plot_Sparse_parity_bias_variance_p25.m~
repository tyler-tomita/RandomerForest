%% Plot bias, variance, and generalization error for Trunk simulations

close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

rng(1);

%% Plotting parameters

LineWidth.rf = 2;
LineWidth.rerf = 2;
LineWidth.rf_rot = 2;
LineWidth.frc2 = 2;
C = [0 1 1;0 1 0;1 0 1;1 0 0;0 0 0;1 .5 0];
Colors.rf = C(1,:);
Colors.rerf = C(2,:);
Colors.rf_rot = C(3,:);
Colors.frc2 = C(5,:);

axWidth = 2;
axHeight = 2;
FontSize = 0.2;
axLeft = [FontSize*5,FontSize*10+axWidth,FontSize*15+axWidth*2];
axBottom = repmat(FontSize*3,1,3);
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
    run_Sparse_parity_bias_variance
else
    load Sparse_parity_bias_variance
end

load('Sparse_parity_partitioned_data','dims')

Results.GE = GE;
Results.V = V;
Results.B = B;
MetricNames.GE = 'Generalization Error';
MetricNames.V = 'Variance';
MetricNames.B = 'Bias^2';

Classifiers = fieldnames(Params);
Classifiers(end-3:end) = [];
Metrics = fieldnames(Results);

%% Plot each metric vs mtry
i = 4;
d = dims(i);
for j = 1:length(Metrics)
    Metric = Metrics{j};
    ax = subplot(1,length(Metrics),j);
    hold on
    for k = 1:length(Classifiers)
        cl = Classifiers{k};
        if ~isempty(Results.(Metric).(cl){i})
            plot(Params.(cl).mtry{i},Results.(Metric).(cl){i},...
                'LineWidth',LineWidth.(cl),'Color',Colors.(cl))
        end
    end
    if j == 2
        title(sprintf('Sparse Parity (n = 1000, p = %d)',d))
    end
    xlabel('d')
    ylabel(MetricNames.(Metric))
    if j == length(Metrics)
        l = legend(Classifiers);
        l.Location = 'southwest';
        l.Box = 'off';
    end
    
    ax.LineWidth = 2;
    ax.FontUnits = 'inches';
    ax.FontSize = FontSize;
    ax.Units = 'inches';
    ax.Position = [axLeft(j) axBottom(j) axWidth axHeight];
    ax.Box = 'off';
    if d > 5
        ax.XScale = 'log';
    end
end

save_fig(gcf,[rerfPath 'RandomerForest/Figures/Sparse_parity_bias_variance_p25'])