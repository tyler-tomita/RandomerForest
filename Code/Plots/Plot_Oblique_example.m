% Plot oblique example
close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

LineWidth = 2;
MarkerSize = 12;
FontSize = .4;
axWidth = 1.5;
axHeight = 1.5;
axLeft = repmat(FontSize*4,1,3);
axBottom = [FontSize*10+axHeight*2,FontSize*6.5+axHeight,FontSize*3];
cbHeight = axHeight;
cbWidth = axWidth/8;
cbLeft = axLeft(end) + axWidth + FontSize/2;
cbBottom = axBottom(1);
figWidth = cbLeft(end) + cbWidth + FontSize*3;
figHeight = axBottom(1) + axHeight + FontSize*2;

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
ax.LineWidth = LineWidth;
ax.FontUnits = 'inches';
ax.FontSize = FontSize;
ax.Units = 'inches';
ax.Position = [axLeft(1) axBottom(1) axWidth axHeight];
ax.Box = 'off';

save_fig(gcf,[rerfPath 'RandomerForest/Figures/Oblique_example_posteriors_1'],...
    {'fig','pdf','png'})

ax = axes;

posterior_map(Xpost,Ypost,Phats{2}.rf,true,ColorMap2);

% hold on;
% plot(x0(:,1),x0(:,2),'.','MarkerEdgeColor',ColorMap(2,:),'MarkerSize',MarkerSize)
% plot(x1(:,1),x1(:,2),'.','MarkerEdgeColor',ColorMap(end-1,:),'MarkerSize',MarkerSize)

xlabel('X_1')
ylabel('X_2')
ax.LineWidth = LineWidth;
ax.FontUnits = 'inches';
ax.FontSize = FontSize;
ax.Units = 'inches';
ax.Position = [axLeft(2) axBottom(2) axWidth axHeight];
ax.Box = 'off';

save_fig(gcf,[rerfPath 'RandomerForest/Figures/Oblique_example_posteriors_2'],...
    {'fig','pdf','png'})

ax = axes;

% posterior_map(Xpost,Ypost,Phats{2}.rerf,false,ColorMap2);
posterior_map(Xpost,Ypost,Phats{2}.frc,true,ColorMap2);

% hold on;
% plot(x0(:,1),x0(:,2),'.','MarkerEdgeColor',ColorMap(2,:),'MarkerSize',MarkerSize)
% plot(x1(:,1),x1(:,2),'.','MarkerEdgeColor',ColorMap(end-1,:),'MarkerSize',MarkerSize)

xlabel('X_1')
ylabel('X_2')
ax.LineWidth = LineWidth;
ax.FontUnits = 'inches';
ax.FontSize = FontSize;
ax.Units = 'inches';
ax.Position = [axLeft(3) axBottom(3) axWidth axHeight];
ax.Box = 'off';

save_fig(gcf,[rerfPath 'RandomerForest/Figures/Oblique_example_posteriors_3'],...
    {'fig','pdf','png'})