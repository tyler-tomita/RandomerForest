clear
close all
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

load('purple2green')

Colors(1,:) = ColorMap(3,:);
Colors(2,:)= ColorMap(9,:);
LineWidth = 2;
MarkerSize = 12;
FontSize = .2;
axWidth = 2;
axHeight = 1.5;
axLeft = repmat(FontSize*4,1,3);
axBottom = [...
    FontSize*14+axHeight*2,FontSize*9+axHeight,...
    FontSize*4];
legWidth = axWidth;
legHeight = axHeight;
legLeft = axLeft(end) + axWidth + FontSize;
legBottom = axBottom(2);
figWidth = legLeft + legWidth;
figHeight = axBottom(1) + axHeight + FontSize*2.5;

fig = figure;
fig.Units = 'inches';
fig.PaperUnits = 'inches';
fig.Position = [0 0 figWidth figHeight];
fig.PaperPosition = [0 0 figWidth figHeight];
fig.PaperSize = [figWidth figHeight];

% Plot Sparse Parity

% load([rerfPath 'RandomerForest/Results/pami/Sparse_parity/mat/Sparse_parity_vary_n.mat'])
S1 = load([rerfPath 'RandomerForest/Results/2017.04.13/Sparse_parity/Sparse_parity_raw_vary_n_aggregated.mat']);
S2 = load([rerfPath 'RandomerForest/Results/2017.05.22/Sparse_parity_raw_vary_n_rerf.2017.05.22.mat']);

ntrials = length(S1.TestError{1}.rf);
% Classifiers = [fieldnames(TestError{1});'xgb'];
Classifiers = {'rerf'};

for j = 1
    p = S1.ps(j);
    
    ErrorMatrix = NaN(ntrials,length(S1.ns{j}),2);
    BiasMatrix = NaN(1,length(S1.ns{j}),2);
    VarianceMatrix = NaN(1,length(S1.ns{j}),2);
    
    for i = 1:length(S1.ns{j})
        n = S1.ns{j}(i);        
        cl = 'rerf';
        k = 1;
        if ~isempty(S1.TestError{i,j}.(cl))
            OE = S1.OOBError{i,j}.(cl);
            OA = S1.OOBAUC{i,j}.(cl);
            ntrials = size(OE,1);
            TE = NaN(ntrials,1);
            for trial = 1:ntrials
                B = hp_optimize(OE(trial,1:length(S1.Params{i,j}.(cl).d)),OA(trial,1:length(S1.Params{i,j}.(cl).d)));
                BI(trial) = B(end);
                TE(trial) = S1.TestError{i,j}.(cl)(trial,BI(trial));
            end
            ErrorMatrix(:,i,k) = TE;
%                     ErrorMatrix(:,j,c) = TestError{i,j}.(cl)';
            BiasMatrix(1,i,k) = mean(S1.Bias{i,j}.(cl)(BI))';
            VarianceMatrix(1,i,k) = mean(S1.Variance{i,j}.(cl)(BI))';
        else
            ErrorMatrix(:,i,k) = NaN;
            BiasMatrix(1,i,k) = NaN;
            VarianceMatrix(1,i,k) = NaN;
        end
        
        k = 2;
        if ~isempty(S2.TestError{i,j}.(cl))
            OE = S2.OOBError{i,j}.(cl);
            OA = S2.OOBAUC{i,j}.(cl);
            ntrials = size(OE,1);
            TE = NaN(ntrials,1);
            for trial = 1:ntrials
                B = hp_optimize(OE(trial,1:length(S2.Params{i,j}.(cl).d)),OA(trial,1:length(S2.Params{i,j}.(cl).d)));
                BI(trial) = B(end);
                TE(trial) = S2.TestError{i,j}.(cl)(trial,BI(trial));
            end
            ErrorMatrix(:,i,k) = TE;
%                     ErrorMatrix(:,j,c) = TestError{i,j}.(cl)';
            BiasMatrix(1,i,k) = mean(S2.Bias{i,j}.(cl)(BI))';
            VarianceMatrix(1,i,k) = mean(S2.Variance{i,j}.(cl)(BI))';
        else
            ErrorMatrix(:,i,k) = NaN;
            BiasMatrix(1,i,k) = NaN;
            VarianceMatrix(1,i,k) = NaN;
        end
    end
end

% bias

j = 1;
ax(j) = axes;
hold on
ymax = max(BiasMatrix(:));
ymin = min(BiasMatrix(:));
for c = 1:2
    plot(S1.ns{1},BiasMatrix(1,:,c),...
        'LineWidth',LineWidth,'Color',Colors(c,:))
end

ax(j).LineWidth = LineWidth;
ax(j).FontUnits = 'inches';
ax(j).FontSize = FontSize;
ax(j).Units = 'inches';
ax(j).Position = [axLeft(j) axBottom(j) axWidth axHeight];
ax(j).Box = 'off';
ax(j).XLim = [min(S1.ns{1})-1 max(S1.ns{1})+1];
ax(j).XScale = 'log';
ax(j).XTick = S1.ns{1};
ax(j).XTickLabel = cellstr(num2str(S1.ns{1}'))';
ax(j).XTickLabelRotation = 0;
ax(j).YLim = [ymin ymax];

xlabel('n')
ylabel('Bias')
if j == 1
    title({'Sparse Parity';sprintf('p = %d',3)})
end

% variance

j = 2;
ax(j) = axes;
hold on
ymax = max(VarianceMatrix(:));
ymin = min(VarianceMatrix(:));
for c = 1:2
    plot(S1.ps,VarianceMatrix(1,:,c),...
        'LineWidth',LineWidth,'Color',Colors(c,:))
end

ax(j).LineWidth = LineWidth;
ax(j).FontUnits = 'inches';
ax(j).FontSize = FontSize;
ax(j).Units = 'inches';
ax(j).Position = [axLeft(j) axBottom(j) axWidth axHeight];
ax(j).Box = 'off';
ax(j).XLim = [min(S1.ns{1})-1 max(S1.ns{1})+1];
ax(j).XScale = 'log';
ax(j).XTick = S1.ns{1};
ax(j).XTickLabel = cellstr(num2str(S1.ns{1}'))';
ax(j).XTickLabelRotation = 0;
ax(j).YLim = [ymin ymax];

xlabel('n')
ylabel('Variance')

lh = legend('RerF-PET','RerF-CT');
lh.Box = 'off';
lh.Units = 'inches';
lh.Position = [legLeft legBottom legWidth legHeight];

% error rate

j = 3;
ax(j) = axes;
hold on
ymax = zeros(length(Classifiers),1);
ymin = zeros(length(Classifiers),1);
for c = 1:2
    mn = mean(ErrorMatrix(:,:,c));
    sem = std(ErrorMatrix(:,:,c))/sqrt(ntrials);
    errorbar(S1.ns{1},mn,sem,...
        'LineWidth',LineWidth,'Color',Colors(c,:))
    ymax(c) = mn(1)+2*sem(1);
    ymin(c) = mn(end)-2*sem(end);
end

ax(j).LineWidth = LineWidth;
ax(j).FontUnits = 'inches';
ax(j).FontSize = FontSize;
ax(j).Units = 'inches';
ax(j).Position = [axLeft(j) axBottom(j) axWidth axHeight];
ax(j).Box = 'off';
ax(j).XLim = [min(S1.ns{1})-1 max(S1.ns{1})+1];
ax(j).XScale = 'log';
ax(j).XTick = S1.ns{1};
ax(j).XTickLabel = cellstr(num2str(S1.ns{1}'))';
ax(j).XTickLabelRotation = 0;
ax(j).YLim = [min(ymin) max(ymax)];

xlabel('n')
ylabel('Error Rate')

save_fig(gcf,[rerfPath 'RandomerForest/Figures/2017.05.22/Plot_ct_pet_comparison_2017_05_22'],{'fig','pdf','png'})