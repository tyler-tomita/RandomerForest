%% Plot Sparse Parity and Trunk

clear
close all
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

load('purple2green')
Colors.rf = ColorMap(1,:);
Colors.rfr = ColorMap(1,:);
Colors.frc= ColorMap(9,:);
Colors.frcr = ColorMap(9,:);
Colors.rr_rf = ColorMap(3,:);
Colors.rr_rfr = ColorMap(3,:);
Colors.Class0 = ColorMap(3,:);
Colors.Class1= ColorMap(9,:);
% Colors.rf = 'c';
% Colors.rfr = 'c';
% Colors.rerf = 'b';
% Colors.rerfr = 'b';
% Colors.frc = 'g';
% Colors.frcr = 'g';
% Colors.rr_rf = 'm';
% Colors.rr_rfr = 'm';
LineStyles.rf = '-';
LineStyles.rfr = ':';
LineStyles.frc = '-';
LineStyles.frcr = ':';
LineStyles.rr_rf = '-';
LineStyles.rr_rfr = ':';
LineWidth = 2;
MarkerSize = 12;
FontSize = .25;
axWidth = 1.5;
axHeight = 1.5;
axLeft = [FontSize*3.5,FontSize*13+axWidth];
axBottom = FontSize*3*ones(1,2);
legWidth = axWidth;
legHeight = axHeight/2;
legLeft = axLeft + axWidth*3/4 + FontSize;
legBottom = axBottom + axHeight/2 - legHeight/2;
figWidth = legLeft(end) + legWidth + FontSize/2;
figHeight = axBottom(1) + axHeight + FontSize*1.5;

fig = figure;
fig.Units = 'inches';
fig.PaperUnits = 'inches';
fig.Position = [0 0 figWidth figHeight];
fig.PaperPosition = [0 0 figWidth figHeight];
fig.PaperSize = [figWidth figHeight];

%% Plot Sparse Parity

ax = axes;

n = 100;
p = 2;
p_prime = min(3,p);

X = rand(n,p)*2 - 1;
Y = mod(sum(X(:,1:p_prime)>0,2),2);

plot(X(Y==1,1),X(Y==1,2),'.','Color',Colors.Class1,'MarkerSize',MarkerSize)
hold on
plot(X(Y==0,1),X(Y==0,2),'.','Color',Colors.Class0,'MarkerSize',MarkerSize)

text(0.5,1.05,'Parity','FontSize',16,'FontWeight','bold','Units',...
    'normalized','HorizontalAlignment','center','VerticalAlignment'...
    ,'bottom')
xlabel('X_1')
ylabel('X_2')
[lh1,objh1] = legend('Class 1','Class 0');
lh1.Units = 'inches';
lh1.Position = [legLeft(1) legBottom(1) legWidth legHeight(1)];
lh1.Box = 'off';
lh1.FontSize = 16;
objh1(4).XData = 0.75*(objh1(3).XData(2) - objh1(3).XData(1)) + objh1(3).XData(1);
objh1(6).XData = 0.75*(objh1(5).XData(2) - objh1(5).XData(1)) + objh1(5).XData(1);
for i = [3,5]
    objh1(i).XData = [(objh1(i).XData(2)-objh1(i).XData(1))*.75+objh1(i).XData(1),objh1(i).XData(2)];
end
ax.LineWidth = LineWidth;
ax.FontUnits = 'inches';
ax.FontSize = FontSize;
ax.Units = 'inches';
ax.Position = [axLeft(1) axBottom(1) axWidth axHeight];
ax.Box = 'off';
ax.XLim = [-1 1];
ax.YLim = [-1 1];
ax.XTick = [-1 0 1];
ax.YTick = [-1 0 1];
ax.XTickLabel = {'-1';'0';'1'};
ax.YTickLabel = {'-1';'0';'1'};

%% Plot Trunk

ax = axes;

p = 100;
d_idx = 1:p;
mu1 = 1./sqrt(d_idx);
mu0 = -1*mu1;

plot(d_idx,mu1,'Color',Colors.Class1,'LineWidth',LineWidth)
hold on
plot(d_idx,mu0,'Color',Colors.Class0,'LineWidth',LineWidth)

text(0.5,1.05,'Trunk','FontSize',16,'FontWeight','bold','Units',...
    'normalized','HorizontalAlignment','center','VerticalAlignment'...
    ,'bottom')
xlabel('Dimension Index')
ylabel('Mean')
ax.LineWidth = LineWidth;
ax.FontUnits = 'inches';
ax.FontSize = FontSize;
ax.Units = 'inches';
ax.Position = [axLeft(2) axBottom(2) axWidth axHeight];
ax.Box = 'off';
% ax.XLim = [-5 5];
% ax.YLim = [-5 5];
% ax.XTick = [-5 0 5];
%ax.YTick = [-5 0 5];
% ax.XTickLabel = {'-5';'0';'5'};
%ax.YTickLabel = {'-5';'0';'5'};
[lh2,objh2] = legend('Class 1','Class 0');
lh2.Units = 'inches';
lh2.Position = [legLeft(2) legBottom(2) legWidth legHeight];
lh2.Box = 'off';
lh2.FontSize = 14;
% objh2(4).XData = 0.75*(objh2(3).XData(2) - objh2(3).XData(1)) + objh2(3).XData(1);
% objh2(6).XData = 0.75*(objh2(5).XData(2) - objh2(5).XData(1)) + objh2(5).XData(1);
for i = [3,5]
    objh2(i).XData = [(objh2(i).XData(2)-objh2(i).XData(1))*.75+objh2(i).XData(1),objh2(i).XData(2)];
end


save_fig(gcf,[rerfPath 'RandomerForest/Figures/ROFLMAO_synthetic_data'],{'fig','pdf','png'})