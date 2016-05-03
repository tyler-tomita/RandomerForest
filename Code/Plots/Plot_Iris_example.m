close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

load fisheriris

X = meas;
Y = grp2idx(species);

LineWidth = 2;
MarkerSize = 12;
FontSize = .2;
axHeight = 2;
axWidth = 2;
axLeft = FontSize*4;
axBottom = FontSize*4;

figWidth = axLeft + axWidth + FontSize*4;
figHeight = axBottom + axHeight + FontSize*4;

fig = figure;
fig.Units = 'inches';
fig.PaperUnits = 'inches';
fig.Position = [0 0 figWidth figHeight];
fig.PaperPosition = [0 0 figWidth figHeight];
fig.PaperSize = [figWidth figHeight];

plot(X(Y==1,1),X(Y==1,2),'.c',X(Y==3,1),X(Y==3,2),'.m',...
    'MarkerSize',MarkerSize)

xlabel('Sepal Length')
ylabel('Sepal Width')
title('Iris Dataset')
l = legend('Setosa','Virginica');
l.Box = 'off';

ax = gca;
ax.LineWidth = LineWidth;
ax.Units = 'inches';
ax.FontUnits = 'inches';
ax.FontSize = FontSize;
ax.Position = [axLeft,axBottom,axWidth,axHeight];
ax.YLim = [2,4.5];
ax.YTick = 2:4;
ax.Box = 'off';

save_fig(gcf,[rerfPath,'RandomerForest/Figures/Iris_example'])