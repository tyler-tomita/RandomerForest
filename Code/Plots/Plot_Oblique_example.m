close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

LineWidth = 2;
MarkerSize = 12;
FontSize = .2;
axWidth = 1.5;
axHeight = 1.5;
axLeft = [mean([FontSize*5,FontSize*9+axWidth]),...
    repmat([FontSize*5,FontSize*9+axWidth],1,2)];
axBottom = [FontSize*9+axHeight*2,(FontSize*6+axHeight)*ones(1,2),FontSize*3*ones(1,2)];

figWidth = axLeft(end) + axWidth + FontSize;
figHeight = axBottom(1) + axHeight + FontSize*1.5;

fig = figure;
fig.Units = 'inches';
fig.PaperUnits = 'inches';
fig.Position = [0 0 figWidth figHeight];
fig.PaperPosition = [0 0 figWidth figHeight];
fig.PaperSize = [figWidth figHeight];

load Oblique_example_posteriors
load Oblique_example_data

%Compute Bayes posterior map

Posteriors = zeros(length(Xpost),2);
Posteriors(:,2) = double(Ypost - Xpost > 0);
Posteriors(:,1) = 1 - Posteriors(:,2);

isPos = strcmp(Ytrain{2},'1');
x0 = Xtrain{2}(~isPos,:);
x1 = Xtrain{2}(isPos,:);

load('purple2green')
ColorMap2 = interpolate_colormap(ColorMap(4:end-3,:),64,true);

ax = axes;

posterior_map(Xpost,Ypost,Posteriors,false,ColorMap2);

xlabel('X_1')
ylabel('X_2')
title('Posterior')
ax.LineWidth = LineWidth;
ax.FontUnits = 'inches';
ax.FontSize = FontSize;
ax.Units = 'inches';
ax.Position = [axLeft(1) axBottom(1) axWidth axHeight];
ax.Box = 'off';

ax = axes;

posterior_map(Xpost,Ypost,Phats{1}.rf,false,ColorMap2);

% hold on;
% plot(x0(:,1),x0(:,2),'.','MarkerEdgeColor',ColorMap(2,:),'MarkerSize',MarkerSize)
% plot(x1(:,1),x1(:,2),'.','MarkerEdgeColor',ColorMap(end-1,:),'MarkerSize',MarkerSize)

xlabel('X_1')
ylabel({['\bf{n = ' num2str(size(Xtrain{1},1)) '}'];'\rm{X_2}'})
title('RF')
ax.LineWidth = LineWidth;
ax.FontUnits = 'inches';
ax.FontSize = FontSize;
ax.Units = 'inches';
ax.Position = [axLeft(2) axBottom(2) axWidth axHeight];
ax.Box = 'off';

ax = axes;

posterior_map(Xpost,Ypost,Phats{1}.rerf,false,ColorMap2);

% hold on;
% plot(x0(:,1),x0(:,2),'.','MarkerEdgeColor',ColorMap(2,:),'MarkerSize',MarkerSize)
% plot(x1(:,1),x1(:,2),'.','MarkerEdgeColor',ColorMap(end-1,:),'MarkerSize',MarkerSize)

xlabel('X_1')
ylabel('X_2')
title('RerF')
ax.LineWidth = LineWidth;
ax.FontUnits = 'inches';
ax.FontSize = FontSize;
ax.Units = 'inches';
ax.Position = [axLeft(3) axBottom(3) axWidth axHeight];
ax.Box = 'off';

ax = axes;

posterior_map(Xpost,Ypost,Phats{2}.rf,false,ColorMap2);

% hold on;
% plot(x0(:,1),x0(:,2),'.','MarkerEdgeColor',ColorMap(2,:),'MarkerSize',MarkerSize)
% plot(x1(:,1),x1(:,2),'.','MarkerEdgeColor',ColorMap(end-1,:),'MarkerSize',MarkerSize)

xlabel('X_1')
ylabel({['\bf{n = ' num2str(size(Xtrain{2},1)) '}'];'\rm{X_2}'})
ax.LineWidth = LineWidth;
ax.FontUnits = 'inches';
ax.FontSize = FontSize;
ax.Units = 'inches';
ax.Position = [axLeft(4) axBottom(4) axWidth axHeight];
ax.Box = 'off';

ax = axes;

posterior_map(Xpost,Ypost,Phats{2}.rerf,false,ColorMap2);

% hold on;
% plot(x0(:,1),x0(:,2),'.','MarkerEdgeColor',ColorMap(2,:),'MarkerSize',MarkerSize)
% plot(x1(:,1),x1(:,2),'.','MarkerEdgeColor',ColorMap(end-1,:),'MarkerSize',MarkerSize)

xlabel('X_1')
ylabel('X_2')
ax.LineWidth = LineWidth;
ax.FontUnits = 'inches';
ax.FontSize = FontSize;
ax.Units = 'inches';
ax.Position = [axLeft(5) axBottom(5) axWidth axHeight];
ax.Box = 'off';

save_fig(gcf,[rerfPath 'RandomerForest/Figures/Oblique_example_posteriors'],...
    {'fig','pdf','png'})