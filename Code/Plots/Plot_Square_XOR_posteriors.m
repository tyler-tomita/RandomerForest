%% Plot posterior heat maps

clear
close all
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

rng(1);

n = 500;
d = 2;

x1 = rand(n/4,2);
x2 = [-rand(n/4,1),rand(n/4,1)];
x3 = -rand(n/4,2);
x4 = [rand(n/4,1),-rand(n/4,1)];
XX = [x1;x2;x3;x4];
Y = double((XX(:,1)<0 & XX(:,2)>0)...
    | (XX(:,1)>0 & XX(:,2)<0));
theta = pi/4;
Xrot = XX*[cos(theta),sin(theta);-sin(theta),cos(theta)];
X.train{1} = XX;
X.train{2} = Xrot;

LineWidth = 2;
MarkerSize = 12;
FontSize = .2;
axWidth = 1.5;
axHeight = 1.5;
cbWidth = .2;
cbHeight = axHeight;
axLeft = [FontSize*4,FontSize*10+axWidth+cbWidth,FontSize*4,...
    FontSize*10+axWidth+cbWidth,FontSize*4,FontSize*10+axWidth+cbWidth];
axBottom = [FontSize*9+axHeight*2 FontSize*9+axHeight*2,FontSize*6+axHeight,...
    FontSize*6+axHeight,FontSize*3,FontSize*3];
cbLeft = axLeft + axWidth + FontSize;
cbBottom = axBottom;
figWidth = cbLeft(end) + cbWidth + FontSize*4;
figHeight = axBottom(1) + axHeight + FontSize*3;

fig = figure;
fig.Units = 'inches';
fig.PaperUnits = 'inches';
fig.Position = [0 0 figWidth figHeight];
fig.PaperPosition = [0 0 figWidth figHeight];
fig.PaperSize = [figWidth figHeight];

runSims = false;

if runSims
    run_Square_XOR_posteriors
else
    load Square_XOR_posteriors_npoints_50.mat
end

%%

ax = subplot(3,2,1);
plot(X.train{1}(Y==0,1),X.train{1}(Y==0,2),'b.',...
    X.train{1}(Y==1,1),X.train{1}(Y==1,2),'r.',...
    'MarkerSize',MarkerSize)
xlabel('X_1')
ylabel('X_2')
title('XOR')
ax.XLim = [-1,1];
ax.YLim = [-1,1];
ax.LineWidth = LineWidth;
ax.FontUnits = 'inches';
ax.FontSize = FontSize;
ax.Units = 'inches';
ax.Position = [axLeft(1) axBottom(1) axWidth axHeight];
ax.XTick = [];
ax.YTick = ax.XTick;
ax.TickDir = 'out';
ax.TickLength = [.02 .03];
ax.Box = 'off';

%%

ax = subplot(3,2,2);
plot(X.train{2}(Y==0,1),X.train{2}(Y==0,2),'b.',...
    X.train{2}(Y==1,1),X.train{2}(Y==1,2),'r.',...
    'MarkerSize',MarkerSize)
xlabel('X_1')
ylabel('X_2')
title('Rotated XOR')
l = legend('Class 0', 'Class 1');
l.Box = 'off';
ax.XLim = [-sqrt(2),sqrt(2)];
ax.YLim = [-sqrt(2),sqrt(2)];
ax.LineWidth = LineWidth;
ax.FontUnits = 'inches';
ax.FontSize = FontSize;
ax.Units = 'inches';
ax.Position = [axLeft(2) axBottom(2) axWidth axHeight];
ax.XTick = [];
ax.YTick = ax.XTick;
ax.TickDir = 'out';
ax.TickLength = [.02 .03];
ax.Box = 'off';

%%

ax(1) = subplot(3,2,3);
p1 = posterior_map(Xpost{1},Ypost{1},Posteriors.rf{1});
xlabel('X_1')
ylabel({'RF';'X_2'})
ax(1).XLim = [-1,1];
ax(1).YLim = [-1,1];
ax(1).LineWidth = LineWidth;
ax(1).FontUnits = 'inches';
ax(1).FontSize = FontSize;
ax(1).Units = 'inches';
ax(1).Position = [axLeft(3) axBottom(3) axWidth axHeight];
ax(1).XTick = [];
ax(1).YTick = ax(1).XTick;
ax(1).TickDir = 'out';
ax(1).TickLength = [.02 .03];

cb(1) = colorbar;
cb(1).Units = 'inches';
cb(1).Position = [cbLeft(3) cbBottom(3) cbWidth cbHeight];
cb(1).Box = 'off';
colormap(ax(1),'jet')

%%

ax(2) = subplot(3,2,4);
p2 = posterior_map(Xpost{2},Ypost{2},Posteriors.rf{2});
xlabel('X_1')
ylabel('X_2')
ax(2).XLim = [-sqrt(2),sqrt(2)];
ax(2).YLim = [-sqrt(2),sqrt(2)];
ax(2).LineWidth = LineWidth;
ax(2).FontUnits = 'inches';
ax(2).FontSize = FontSize;
ax(2).Units = 'inches';
ax(2).Position = [axLeft(4) axBottom(4) axWidth axHeight];
ax(2).XTick = [];
ax(2).YTick = ax(2).XTick;
ax(2).TickDir = 'out';
ax(2).TickLength = [.02 .03];

cb(2) = colorbar;
cb(2).Units = 'inches';
cb(2).Position = [cbLeft(4) cbBottom(4) cbWidth cbHeight];
cb(2).Box = 'off';
colormap(ax(2),'jet')

cmin = min([p1.CData(:);p2.CData(:)]);
cmax = max([p1.CData(:);p2.CData(:)]);

%%

ax(3) = subplot(3,2,5);
p1 = posterior_map(Xpost{1},Ypost{1},Posteriors.rerf{1});
xlabel('X_1')
ylabel({'RerF';'X_2'})
ax(3).XLim = [-1,1];
ax(3).YLim = [-1,1];
ax(3).LineWidth = LineWidth;
ax(3).FontUnits = 'inches';
ax(3).FontSize = FontSize;
ax(3).Units = 'inches';
ax(3).Position = [axLeft(5) axBottom(5) axWidth axHeight];
ax(3).XTick = [];
ax(3).YTick = ax(3).XTick;
ax(3).TickDir = 'out';
ax(3).TickLength = [.02 .03];

cb(3) = colorbar;
cb(3).Units = 'inches';
cb(3).Position = [cbLeft(5) cbBottom(5) cbWidth cbHeight];
cb(3).Box = 'off';
colormap(ax(3),'jet')

%%

ax(4) = subplot(3,2,6);
p2 = posterior_map(Xpost{2},Ypost{2},Posteriors.rerf{2});
xlabel('X_1')
ylabel('X_2')
ax(4).XLim = [-sqrt(2),sqrt(2)];
ax(4).YLim = [-sqrt(2),sqrt(2)];
ax(4).LineWidth = LineWidth;
ax(4).FontUnits = 'inches';
ax(4).FontSize = FontSize;
ax(4).Units = 'inches';
ax(4).Position = [axLeft(6) axBottom(6) axWidth axHeight];
ax(4).XTick = [];
ax(4).YTick = ax(4).XTick;
ax(4).TickDir = 'out';
ax(4).TickLength = [.02 .03];

cb(4) = colorbar;
cb(4).Units = 'inches';
cb(4).Position = [cbLeft(6) cbBottom(6) cbWidth cbHeight];
cb(4).Box = 'off';
colormap(ax(4),'jet')

cmin = min([p1.CData(:);p2.CData(:)]);
cmax = max([p1.CData(:);p2.CData(:)]);

for i = 1:4
    axes(ax(i))
    caxis([cmin cmax])
    cb(i).Ticks = 0:0.2:1;
end

save_fig(gcf,[rerfPath 'RandomerForest/Figures/Square_XOR_posteriors'])