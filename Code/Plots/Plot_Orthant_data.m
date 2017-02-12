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
axLeft = FontSize*5;
axBottom = FontSize*3;
legWidth = axWidth;
legHeight = axHeight/2;
legLeft = axLeft + axWidth*2/3 + FontSize;
legBottom = axBottom + axHeight/2 - legHeight(1)/2;
figWidth = legLeft + legWidth + FontSize;
figHeight = axBottom + axHeight + FontSize*1.5;

fig = figure;
fig.Units = 'inches';
fig.PaperUnits = 'inches';
fig.Position = [0 0 figWidth figHeight];
fig.PaperPosition = [0 0 figWidth figHeight];
fig.PaperSize = [figWidth figHeight];

load Orthant_data

x1 = Xtrain(strcmp(Ytrain,'1'),:);
x2 = Xtrain(strcmp(Ytrain,'2'),:);
x3 = Xtrain(strcmp(Ytrain,'3'),:);
x4 = Xtrain(strcmp(Ytrain,'4'),:);

ax = axes;
hold on

plot(x1(:,1),x1(:,2),'.','MarkerSize',MarkerSize)
plot(x2(:,1),x2(:,2),'.','MarkerSize',MarkerSize)
plot(x3(:,1),x3(:,2),'.','MarkerSize',MarkerSize)
plot(x4(:,1),x4(:,2),'.','MarkerSize',MarkerSize)

xlabel('X_1')
ylabel('X_2')
title('Orthant')

ax.LineWidth = LineWidth;
ax.FontUnits = 'inches';
ax.FontSize = FontSize;
ax.Units = 'inches';
ax.Position = [axLeft axBottom axWidth axHeight];
ax.Box = 'off';
[lh1,objh1] = legend('Class 1','Class 2','Class 3','Class 4');
lh1.Units = 'inches';
lh1.Position = [legLeft legBottom(1) legWidth legHeight(1)];
lh1.Box = 'off';
lh1.FontSize = 11;
objh1(6).XData = 0.75*(objh1(5).XData(2) - objh1(5).XData(1)) + objh1(5).XData(1);
objh1(8).XData = 0.75*(objh1(7).XData(2) - objh1(7).XData(1)) + objh1(7).XData(1);
objh1(10).XData = 0.75*(objh1(9).XData(2) - objh1(9).XData(1)) + objh1(9).XData(1);
objh1(12).XData = 0.75*(objh1(11).XData(2) - objh1(11).XData(1)) + objh1(11).XData(1);

save_fig(gcf,[rerfPath 'RandomerForest/Figures/Orthant_data'],...
    {'fig','pdf','png'})