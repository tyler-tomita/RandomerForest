clear
close all
clc

fpath = mfilename('fullpath');
frcPath = fpath(1:strfind(fpath,'RandomerForest')-1);

load('purple2green')
Colors.rf = ColorMap(1,:);
Colors.rerf = ColorMap(4,:);
Colors.rr_rf = ColorMap(8,:);
Colors.xgb = ColorMap(10,:);

LineStyles.rf = '-';
LineStyles.rerf = '-';
LineStyles.rr_rf = '-';
LineStyles.xgb = '-';
LineWidth = 2;
FontSize = .2;
axWidth = 2;
axHeight = 2;
axLeft = [FontSize*4,FontSize*8+axWidth];
axBottom = FontSize*2.5*ones(1,2);
legWidth = axWidth;
legHeight = axHeight;
legLeft = axLeft(end) + axWidth*2/3 + FontSize;
legBottom = axBottom(end);
figWidth = legLeft + legWidth + FontSize/2;
figHeight = axBottom(1) + axHeight + FontSize*1.5;
% figWidth = axLeft(end) + axWidth + FontSize;
% figHeight = axBottom(1) + axHeight + FontSize*1.5;

fig = figure;
fig.Units = 'inches';
fig.PaperUnits = 'inches';
fig.Position = [0 0 figWidth figHeight];
fig.PaperPosition = [0 0 figWidth figHeight];
fig.PaperSize = [figWidth figHeight];

%% Plot error vs d

load ~/Sparse_parity_vary_n

[nn,np] = size(TestError);

Classifiers = {'rf','rerf','rr_rf'};

ErrorMatrix = NaN(length(Params{1,3}.rerf.d),length(Classifiers));
SEM = NaN(length(Params{1,3}.rerf.d),length(Classifiers));

ax(1) = axes;

for c = 1:length(Classifiers)
    e = OOBError{1,3}.(Classifiers{c})(:,:,end);
    ErrorMatrix(1:length(Params{1,3}.(Classifiers{c}).d),c) = mean(e)';
    SEM(1:length(Params{1,3}.(Classifiers{c}).d),c) = std(e)'/size(e,1);
    errorbar(Params{1,3}.rerf.d,ErrorMatrix(:,c),SEM(:,c),...
        'LineWidth',LineWidth,'Color',Colors.(Classifiers{c}));
    hold on
end

xlabel('d')
ylabel('OOB Error')
    
ax(1).LineWidth = LineWidth;
ax(1).FontUnits = 'inches';
ax(1).FontSize = FontSize;
ax(1).Units = 'inches';
ax(1).Position = [axLeft(1) axBottom(1) axWidth axHeight];
ax(1).Box = 'off';
% ax(1).XLim = [ns{j}(1) ns{j}(end)];
ax(1).XScale = 'log';
% ax(1).XTick = ns{j};
ax(1).YLim = [0.05 0.5];
    

clear ErrorMatrix SEM
%% Plot training time vs d

TimeMatrix = NaN(length(Params{1,3}.rerf.d),length(Classifiers));
SEM = NaN(length(Params{1,3}.rerf.d),length(Classifiers));

ax(2) = axes;

for c = 1:length(Classifiers)
    e = TrainTime{1,3}.(Classifiers{c});
    TimeMatrix(1:length(Params{1,3}.(Classifiers{c}).d),c) = mean(e)';
    SEM(1:length(Params{1,3}.(Classifiers{c}).d),c) = std(e)'/size(e,1);
    errorbar(Params{1,3}.rerf.d,TimeMatrix(:,c),SEM(:,c),...
        'LineWidth',LineWidth,'Color',Colors.(Classifiers{c}));
    hold on
end

xlabel('d')
ylabel('TrainTime (s)')
    
ax(2).LineWidth = LineWidth;
ax(2).FontUnits = 'inches';
ax(2).FontSize = FontSize;
ax(2).Units = 'inches';
ax(2).Position = [axLeft(2) axBottom(2) axWidth axHeight];
ax(2).Box = 'off';
% ax(2).XLim = [ns{j}(1) ns{j}(end)];
ax(2).XScale = 'log';
% ax(2).XTick = ns{j};
ax(2).YLim = [5 430];
ax(2).YScale = 'log';
    
[lh,objh] = legend('RF','RerF','RR-RF');
lh.Box = 'off';
lh.FontSize = 14;
lh.Units = 'inches';
lh.Position = [legLeft legBottom legWidth legHeight];
for k = 4:length(objh)
    objh(k).Children.Children(2).XData = [(objh(k).Children.Children(2).XData(2)-objh(k).Children.Children(2).XData(1))*.75+objh(k).Children.Children(2).XData(1),objh(k).Children.Children(2).XData(2)];
end

save_fig(gcf,[frcPath 'RandomerForest/Figures/Sparse_parity_error_and_time_vs_d'])