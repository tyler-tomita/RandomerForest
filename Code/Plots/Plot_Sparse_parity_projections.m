close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

load('purple2green')
Colors.Class0 = ColorMap(3,:);
Colors.Class1= ColorMap(9,:);

MarkerSize = 10;
LineWidth = 2;
FontSize = .18;
axWidth = 3.75;
axHeight = 1.25;
% axLeft = FontSize*6*ones(1,6);
% axBottom = [FontSize*15+axHeight*5,FontSize*14+axHeight*4,FontSize*9+axHeight*3,...
%     FontSize*8+axHeight*2,FontSize*3+axHeight,...
%     FontSize*2];
axLeft = FontSize*3*ones(1,3);
axBottom = [FontSize*8+axHeight*2,FontSize*5+axHeight,FontSize*2];
% titleAxesLeft = axLeft(1);
% titleAxesBottom = [mean([axBottom(1),axBottom(2)]),...
%     mean([axBottom(3),axBottom(4)]),...
%     mean([axBottom(5),axBottom(6)])];
legWidth = 1;
legHeight = axHeight/2;
legLeft = axLeft(2) + axWidth - 2*FontSize;
legBottom = axBottom(2) - (legHeight - axHeight)/2;
figWidth = legLeft + legWidth + FontSize*0.1;
figHeight = axBottom(1) + axHeight + FontSize*2;

fig = figure;
fig.Color = [1 1 1];
fig.Units = 'inches';
fig.PaperUnits = 'inches';
fig.Position = [0 0 figWidth figHeight];
fig.PaperPosition = [0 0 figWidth figHeight];
fig.PaperSize = [figWidth figHeight];

X = dlmread('~/Documents/R/Data/Sparse_parity/dat/Train/Sparse_parity_train_set_n1000_p20_trial7.dat');
Y = X(:,end);
X(:,end) = [];
rmIdx = any(abs(X(:,1:3)) <= 0.25,2);
Y(rmIdx) = [];
X(rmIdx,:) = [];
[n,p] = size(X);
Yunique = [0;1];

% compute impurity

ClProb(2) = sum(Y)/n;
ClProb(1) = 1 - ClProb(2);
I = sum(ClProb.*(1 - ClProb));

%% Best RF split direction %%

Xp = X;
[Xsort,SortIdx] = sort(Xp);
Ysort = Y(SortIdx);
BestVar = find_best_split_mex(Xsort,Ysort,Yunique,I);
x = Xp(:,BestVar);
x = rescale(x,[],'zscore');
xmn = round(min(x))-1;
xmx = round(max(x))+1;
xi = linspace(xmn,xmx,100);
[f0,~] = ksdensity(x(Y==0),xi);
[f1,~] = ksdensity(x(Y==1),xi);

% ax(1) = axes;
% x0 = x(Y==0);
% x1 = x(Y==1);
% plot(x0,zeros(size(x0)),'.','MarkerSize',MarkerSize,'MarkerEdgeColor',Colors.Class0)
% hold on
% plot(x1,ones(size(x1)),'.','MarkerSize',MarkerSize,'MarkerEdgeColor',Colors.Class1)
% % xlabel('Projection onto best first split direction')
% ylabel('Class')
% title('Sparse Parity (p = 20)')
% ax(1).XLim = [xmn,xmx];
% ax(1).XTick = [];
% % ax(1).XColor = ax(1).Color;
% ax(1).YLim = [-1,2];
% ax(1).YTick = [0,1];
% ax(1).Box = 'off';
% ax(1).LineWidth = LineWidth;
% ax(1).FontUnits = 'inches';
% ax(1).FontSize = FontSize;
% ax(1).Units = 'inches';
% ax(1).Position = [axLeft(1) axBottom(1) axWidth axHeight];
% text(-0.015,1.05,'(A)','FontSize',12,'Units','normalized','HorizontalAlignment','right',...
%     'VerticalAlignment','bottom')
% 
% % create invisible axes in order to make a centered side title
% 
% titleAxes(1) = axes;
% titleAxes(1).FontUnits = 'inches';
% titleAxes(1).FontSize = FontSize;
% titleAxes(1).Units = 'inches';
% titleAxes(1).Position = [titleAxesLeft titleAxesBottom(1) axWidth axHeight];
% text(-0.15,0.5,'RF','FontSize',12,'FontWeight','bold','Units',...
%     'normalized','HorizontalAlignment','center','VerticalAlignment'...
%     ,'middle','Rotation',90)
% titleAxes(1).Visible = 'off';
% set(findall(titleAxes(1), 'type', 'text'), 'visible', 'on')

ax(1) = axes;
plot(xi,f0,'LineWidth',LineWidth,'Color',Colors.Class0)
hold on
plot(xi,f1,'LineWidth',LineWidth,'Color',Colors.Class1)
xlabel('Projection onto best first split direction')
title('Sparse Parity (p = 20)')
ylabel('Density Estimate')
ax(1).XLim = [xmn,xmx];
ax(1).XTick = linspace(xmn,xmx,3);
ax(1).YTick = [];
ax(1).Box = 'off';
ax(1).LineWidth = LineWidth;
ax(1).FontUnits = 'inches';
ax(1).FontSize = FontSize;
ax(1).Units = 'inches';
ax(1).Position = [axLeft(1) axBottom(1) axWidth axHeight];
text(-0.1,0.5,'RF','FontSize',12,'FontWeight','bold','Units',...
    'normalized','HorizontalAlignment','center','VerticalAlignment'...
    ,'middle','Rotation',90)

[lh,objh] = legend('Class 0','Class 1');
lh.Box = 'off';
lh.FontSize = 12;
lh.Units = 'inches';
lh.Position = [legLeft legBottom legWidth legHeight];

for i = 3:2:length(objh)
    objh(i).XData = [(objh(i).XData(2)-objh(i).XData(1))*.75+objh(i).XData(1),objh(i).XData(2)];
end

% annotation('line',[0,1],[2/3,2/3],'LineWidth',LineWidth)

%% Best F-RC split direction %%

R = randmat(p,p^3,'frc',3/p,3);
Xp = X*R;
DeltaImpurity = zeros(1,size(Xp,2));
[Xsort,SortIdx] = sort(Xp);
Ysort = Y(SortIdx);
BestVar = find_best_split_mex(Xsort,Ysort,Yunique,I);
x = Xp(:,BestVar);
x = rescale(x,[],'zscore');
xmn = round(min(x))-1;
xmx = round(max(x))+1;
xi = linspace(xmn,xmx,100);
[f0,~] = ksdensity(x(Y==0),xi);
[f1,~] = ksdensity(x(Y==1),xi);

% ax(3) = axes;
% x0 = x(Y==0);
% x1 = x(Y==1);
% plot(x0,zeros(size(x0)),'.','MarkerSize',MarkerSize,'MarkerEdgeColor',Colors.Class0)
% hold on
% plot(x1,ones(size(x1)),'.','MarkerSize',MarkerSize,'MarkerEdgeColor',Colors.Class1)
% ax(3).XLim = [xmn,xmx];
% ax(3).XTick = [];
% % ax(3).XColor = ax(3).Color;
% ax(3).YLim = [-1,2];
% ax(3).YTick = [0,1];
% ax(3).Box = 'off';
% ax(3).LineWidth = LineWidth;
% ax(3).FontUnits = 'inches';
% ax(3).FontSize = FontSize;
% ax(3).Units = 'inches';
% ax(3).Position = [axLeft(3) axBottom(3) axWidth axHeight];
% text(-0.015,1.05,'(B)','FontSize',12,'Units','normalized','HorizontalAlignment','right',...
%     'VerticalAlignment','bottom')
% 
% % create invisible axes in order to make a centered side title
% 
% titleAxes(2) = axes;
% titleAxes(2).FontUnits = 'inches';
% titleAxes(2).FontSize = FontSize;
% titleAxes(2).Units = 'inches';
% titleAxes(2).Position = [titleAxesLeft titleAxesBottom(2) axWidth axHeight];
% text(-0.15,0.5,'F-RC','FontSize',12,'FontWeight','bold','Units',...
%     'normalized','HorizontalAlignment','center','VerticalAlignment'...
%     ,'middle','Rotation',90)
% titleAxes(2).Visible = 'off';
% set(findall(titleAxes(2), 'type', 'text'), 'visible', 'on')

ax(2) = axes;
plot(xi,f0,'LineWidth',LineWidth,'Color',Colors.Class0)
hold on
plot(xi,f1,'LineWidth',LineWidth,'Color',Colors.Class1)
ax(2).XLim = [xmn,xmx];
ax(2).XTick = linspace(xmn,xmx,3);
ax(2).YTick = [];
ax(2).Box = 'off';
ax(2).LineWidth = LineWidth;
ax(2).FontUnits = 'inches';
ax(2).FontSize = FontSize;
ax(2).Units = 'inches';
ax(2).Position = [axLeft(2) axBottom(2) axWidth axHeight];
text(-0.1,0.5,'F-RC','FontSize',12,'FontWeight','bold','Units',...
    'normalized','HorizontalAlignment','center','VerticalAlignment'...
    ,'middle','Rotation',90)

% annotation('line',[0,1],[1/3,1/3],'LineWidth',LineWidth)

%% Best RR-RF split direction

R = random_rotation(p);
Xp = X*R;
[Xsort,SortIdx] = sort(Xp);
Ysort = Y(SortIdx);
BestVar = find_best_split_mex(Xsort,Ysort,Yunique,I);
x = Xp(:,BestVar);
x = rescale(x,[],'zscore');
xmn = round(min(x))-1;
xmx = round(max(x))+1;
xi = linspace(xmn,xmx,100);
[f0,~] = ksdensity(x(Y==0),xi);
[f1,~] = ksdensity(x(Y==1),xi);

% ax(5) = axes;
% x0 = x(Y==0);
% x1 = x(Y==1);
% plot(x0,zeros(size(x0)),'.','MarkerSize',MarkerSize,'MarkerEdgeColor',Colors.Class0)
% hold on
% plot(x1,ones(size(x1)),'.','MarkerSize',MarkerSize,'MarkerEdgeColor',Colors.Class1)
% ax(5).XLim = [xmn,xmx];
% ax(5).XTick = [];
% % ax(5).XColor = ax(1).Color;
% ax(5).YLim = [-1,2];
% ax(5).YTick = [0,1];
% ax(5).Box = 'off';
% ax(5).LineWidth = LineWidth;
% ax(5).FontUnits = 'inches';
% ax(5).FontSize = FontSize;
% ax(5).Units = 'inches';
% ax(5).Position = [axLeft(5) axBottom(5) axWidth axHeight];
% text(-0.015,1.05,'(C)','FontSize',12,'Units','normalized','HorizontalAlignment','right',...
%     'VerticalAlignment','bottom')
% 
% % create invisible axes in order to make a centered side title
% 
% titleAxes(3) = axes;
% titleAxes(3).FontUnits = 'inches';
% titleAxes(3).FontSize = FontSize;
% titleAxes(3).Units = 'inches';
% titleAxes(3).Position = [titleAxesLeft titleAxesBottom(3) axWidth axHeight];
% text(-0.15,0.5,'RR-RF','FontSize',12,'FontWeight','bold','Units',...
%     'normalized','HorizontalAlignment','center','VerticalAlignment'...
%     ,'middle','Rotation',90)
% titleAxes(3).Visible = 'off';
% set(findall(titleAxes(3), 'type', 'text'), 'visible', 'on')

ax(3) = axes;
plot(xi,f0,'LineWidth',LineWidth,'Color',Colors.Class0)
hold on
plot(xi,f1,'LineWidth',LineWidth,'Color',Colors.Class1)
ax(3).XLim = [xmn,xmx];
ax(3).XTick = linspace(xmn,xmx,3);
ax(3).YTick = [];
ax(3).Box = 'off';
ax(3).LineWidth = LineWidth;
ax(3).FontUnits = 'inches';
ax(3).FontSize = FontSize;
ax(3).Units = 'inches';
ax(3).Position = [axLeft(3) axBottom(3) axWidth axHeight];
text(-0.1,0.5,'RR-RF','FontSize',12,'FontWeight','bold','Units',...
    'normalized','HorizontalAlignment','center','VerticalAlignment'...
    ,'middle','Rotation',90)







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close all
% clear
% clc
% 
% fpath = mfilename('fullpath');
% rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);
% 
% % set up rendering and output parameters
% LineWidth = 2;
% MarkerSize = 10;
% FontSize = .2;
% axWidth = 1.5;
% axHeight = 1.5;
% axLeft = [FontSize*3,FontSize*7+axWidth];
% axBottom = FontSize*3*ones(1,2);
% figWidth = axLeft(end) + axWidth + FontSize;
% figHeight = axBottom(1) + axHeight + FontSize*1.5;
% titleLeft = (figWidth - axWidth)/2;
% 
% fig = figure;
% fig.Units = 'inches';
% fig.PaperUnits = 'inches';
% fig.Position = [0 0 figWidth figHeight];
% fig.PaperPosition = [0 0 figWidth figHeight];
% fig.PaperSize = [figWidth figHeight];
% 
% % load sparse parity data and add space between populations to exaggerate
% % separability
% X = dlmread('~/Documents/R/Data/Sparse_parity/dat/Train/Sparse_parity_train_set_n1000_p20_trial1.dat');
% Y = X(:,end);
% X(:,end) = [];
% rmIdx = any(abs(X(:,1:3)) <= 0.25,2);
% Y(rmIdx) = [];
% X(rmIdx,:) = [];
% [n,p] = size(X);
% 
% % plot sum of first three dimensions
% ax = axes;
% pp = 3;
% Xp = X*[ones(pp,1);zeros(p-pp,1)];
% plot(Xp(Y==0),0,'b.',Xp(Y==1),1,'r.','MarkerSize',MarkerSize)
% xlabel('$X_1 + X_2 + X_3$','interpreter','latex')
% ylabel('Class')
% ax.XLim = [-4,4];
% ax.XTick = [-4,0,4];
% ax.YLim = [-1,2];
% ax.YTick = [0,1];
% ax.LineWidth = LineWidth;
% ax.FontUnits = 'inches';
% ax.FontSize = FontSize;
% ax.Units = 'inches';
% ax.Position = [axLeft(1) axBottom(1) axWidth axHeight];
% ax.Box = 'off';
% 
% % plot sum of all 20 dimensions
% ax = axes;
% pp = p;
% Xp = X*[ones(pp,1);zeros(p-pp,1)];
% plot(Xp(Y==0),0,'b.',Xp(Y==1),1,'r.','MarkerSize',MarkerSize)
% xlabel('$X_1 + \cdots + X_{20}$','interpreter','latex')
% ylabel('Class')
% ax.YLim = [-1,2];
% ax.YTick = [0,1];
% ax.LineWidth = LineWidth;
% ax.FontUnits = 'inches';
% ax.FontSize = FontSize;
% ax.Units = 'inches';
% ax.Position = [axLeft(2) axBottom(2) axWidth axHeight];
% ax.Box = 'off';
% 
% % make main title of figure using invisible subplot workaround
% 
% ax = axes;
% title('Sparse Parity (p = 20)')
% ax.FontUnits = 'inches';
% ax.FontSize = FontSize;
% ax.Units = 'inches';
% ax.Position = [titleLeft axBottom(1) axWidth axHeight];
% ax.Visible = 'off';
% set(findall(ax, 'type', 'text'), 'visible', 'on')
% 
save_fig(gcf,[rerfPath 'RandomerForest/Figures/ROFLMAO_fig1_Sparse_parity_projections_2017_01_23'],{'fig','pdf','png'})