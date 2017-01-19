close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

% set up rendering and output parameters
LineWidth = 2;
MarkerSize = 10;
FontSize = .2;
axWidth = 1.5;
axHeight = 1.5;
axLeft = [FontSize*3,FontSize*7+axWidth];
axBottom = FontSize*3*ones(1,2);
figWidth = axLeft(end) + axWidth + FontSize;
figHeight = axBottom(1) + axHeight + FontSize*1.5;
titleLeft = (figWidth - axWidth)/2;

fig = figure;
fig.Units = 'inches';
fig.PaperUnits = 'inches';
fig.Position = [0 0 figWidth figHeight];
fig.PaperPosition = [0 0 figWidth figHeight];
fig.PaperSize = [figWidth figHeight];

% load sparse parity data and add space between populations to exaggerate
% separability
X = dlmread('~/Documents/R/Data/Sparse_parity/dat/Train/Sparse_parity_train_set_n1000_p20_trial1.dat');
Y = X(:,end);
X(:,end) = [];
rmIdx = any(abs(X(:,1:3)) <= 0.25,2);
Y(rmIdx) = [];
X(rmIdx,:) = [];
[n,p] = size(X);

% plot sum of first three dimensions
ax = axes;
pp = 3;
Xp = X*[ones(pp,1);zeros(p-pp,1)];
plot(Xp(Y==0),0,'b.',Xp(Y==1),1,'r.','MarkerSize',MarkerSize)
xlabel('$X_1 + X_2 + X_3$','interpreter','latex')
ylabel('Class')
ax.XLim = [-4,4];
ax.XTick = [-4,0,4];
ax.YLim = [-1,2];
ax.YTick = [0,1];
ax.LineWidth = LineWidth;
ax.FontUnits = 'inches';
ax.FontSize = FontSize;
ax.Units = 'inches';
ax.Position = [axLeft(1) axBottom(1) axWidth axHeight];
ax.Box = 'off';

% plot sum of all 20 dimensions
ax = axes;
pp = p;
Xp = X*[ones(pp,1);zeros(p-pp,1)];
plot(Xp(Y==0),0,'b.',Xp(Y==1),1,'r.','MarkerSize',MarkerSize)
xlabel('$X_1 + \cdots + X_{20}$','interpreter','latex')
ylabel('Class')
ax.YLim = [-1,2];
ax.YTick = [0,1];
ax.LineWidth = LineWidth;
ax.FontUnits = 'inches';
ax.FontSize = FontSize;
ax.Units = 'inches';
ax.Position = [axLeft(2) axBottom(2) axWidth axHeight];
ax.Box = 'off';

% make main title of figure using invisible subplot workaround

ax = axes;
title('Sparse Parity (p = 20)')
ax.FontUnits = 'inches';
ax.FontSize = FontSize;
ax.Units = 'inches';
ax.Position = [titleLeft axBottom(1) axWidth axHeight];
ax.Visible = 'off';
set(findall(ax, 'type', 'text'), 'visible', 'on')

save_fig(gcf,[rerfPath 'RandomerForest/Figures/New_fig1_Sparse_parity_projections'])