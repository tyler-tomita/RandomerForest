%% Plot Performance Profiles for Benchmark Transformations
close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

C = [0 1 1;0 1 0;1 0 1;1 0 0;0 0 0];
Colors.rf = C(1,:);
Colors.rerf = C(2,:);
Colors.rf_rot = C(3,:);
Colors.rerfdn = C(4,:);
LineWidth = 2;
FontSize = .16;
axWidth = 1.5;
axHeight = 1.5;
axLeft = repmat([FontSize*4 FontSize*8+axWidth],4,1);
axBottom = repmat([FontSize*12+axHeight*3;FontSize*9+axHeight*2;FontSize*6+axHeight;FontSize*3],1,2);
figWidth = axLeft(end) + axWidth + FontSize*3;
figHeight = axBottom(1) + axHeight + FontSize*3;
fig = figure;
fig.Units = 'inches';
fig.PaperUnits = 'inches';
fig.Position = [0 0 figWidth figHeight];
fig.PaperPosition = [0 0 figWidth figHeight];
fig.PaperSize = [figWidth figHeight];

Transformations = {'Untransformed' 'Rotated' 'Scaled' 'Affine'};

runSims = false;

if runSims
    run_Performance_profile_untransformed
else
    load Performance_profile_untransformed.mat
end

ax = subplot(4,2,1);
hold on

xmin = min(tau);
xmax = max(tau);
ymin = .7;
ymax = 1.04;

for i = 1:length(clnames)-2
     plot(tau,rho_ps(:,i),'LineWidth',LineWidth,'Color',Colors.(clnames{i}))
end

title('(A)','Units','normalized','Position',[0.025 1],'HorizontalAlignment','left','VerticalAlignment','top')
text(-0.5,0.5,Transformations{1},'FontSize',14,'FontWeight','bold',...
        'Units','normalized','HorizontalAlignment','center',...
        'VerticalAlignment','middle','Rotation',90)
text(0.5,1,'Error Rate','FontSize',14,'FontWeight','bold','Units','normalized','HorizontalAlignment','center','VerticalAlignment','bottom')
xlabel('Relative Performance')
ylabel('Proportion')
ax.LineWidth = LineWidth;
ax.FontUnits = 'inches';
ax.FontSize = FontSize;
ax.Units = 'inches';
ax.Position = [axLeft(1) axBottom(1) axWidth axHeight];
ax.Box = 'off';
ax.XLim = [xmin xmax];
ax.YLim = [ymin ymax];
ax.XTick = [1 10 20 30];
ax.YTick = [.7 .8 .9 1];
l = legend(sprintf('AUC = %0.2f',AUC.rf),['\bf{}' sprintf('AUC = %0.2f',AUC.rerf)],sprintf('AUC = %0.2f',AUC.rf_rot));
l.Location = 'southeast';
l.Box = 'off';


if runSims
    run_Performance_profile_rotate
else
    load Performance_profile_rotate.mat
end

ax = subplot(4,2,3);
hold on

xmin = min(tau);
xmax = max(tau);
ymin = .7;
ymax = 1.04;

for i = 1:length(clnames)-1
     plot(tau,rho_ps(:,i),'LineWidth',LineWidth,'Color',Colors.(clnames{i}))
end

title('(B)','Units','normalized','Position',[0.025 1],'HorizontalAlignment','left','VerticalAlignment','top')
text(-0.5,0.5,Transformations{2},'FontSize',14,'FontWeight','bold',...
        'Units','normalized','HorizontalAlignment','center',...
        'VerticalAlignment','middle','Rotation',90)
ax.LineWidth = LineWidth;
ax.FontUnits = 'inches';
ax.FontSize = FontSize;
ax.Units = 'inches';
ax.Position = [axLeft(2) axBottom(2) axWidth axHeight];
ax.Box = 'off';
ax.XLim = [xmin xmax];
ax.YLim = [ymin ymax];
ax.YTick = [.7 .8 .9 1];
l = legend(sprintf('AUC = %0.2f',AUC.rf),sprintf('AUC = %0.2f',AUC.rerf),['\bf{}' sprintf('AUC = %0.2f',AUC.rf_rot)]);
l.Location = 'southeast';
l.Box = 'off';

if runSims
    run_Performance_profile_scale
else
    load Performance_profile_scale.mat
end

ax = subplot(4,2,5);
hold on

xmin = min(tau);
xmax = max(tau);
ymin = .9;
ymax = 1.02;

for i = 1:length(clnames)-1
     plot(tau,rho_ps(:,i),'LineWidth',LineWidth,'Color',Colors.(clnames{i}))
end

title('(C)','Units','normalized','Position',[0.025 1],'HorizontalAlignment','left','VerticalAlignment','top')
text(-0.5,0.5,Transformations{3},'FontSize',14,'FontWeight','bold',...
        'Units','normalized','HorizontalAlignment','center',...
        'VerticalAlignment','middle','Rotation',90)
ax.LineWidth = LineWidth;
ax.FontUnits = 'inches';
ax.FontSize = FontSize;
ax.Units = 'inches';
ax.Position = [axLeft(3) axBottom(3) axWidth axHeight];
ax.Box = 'off';
ax.XLim = [xmin xmax];
ax.YLim = [ymin ymax];
ax.XTick = [1 15 30];
ax.YTick = [.9 .95 1];
l = legend(sprintf('AUC = %0.2f',AUC.rf),['\bf{}' sprintf('AUC = %0.2f',AUC.rerf)],sprintf('AUC = %0.2f',AUC.rf_rot));
l.Location = 'southeast';
l.Box = 'off';

if runSims
    run_Performance_profile_affine
else
    load Performance_profile_affine.mat
end

ax = subplot(4,2,7);
hold on

xmin = min(tau);
xmax = max(tau);
ymin = .9;
ymax = 1.02;


for i = 1:length(clnames)-1
     plot(tau,rho_ps(:,i),'LineWidth',LineWidth,'Color',Colors.(clnames{i}))
end

title('(D)','Units','normalized','Position',[0.025 1],'HorizontalAlignment','left','VerticalAlignment','top')
text(-0.5,0.5,Transformations{4},'FontSize',14,'FontWeight','bold',...
        'Units','normalized','HorizontalAlignment','center',...
        'VerticalAlignment','middle','Rotation',90)
ax.LineWidth = LineWidth;
ax.FontUnits = 'inches';
ax.FontSize = FontSize;
ax.Units = 'inches';
ax.Position = [axLeft(4) axBottom(4) axWidth axHeight];
ax.Box = 'off';
ax.XLim = [xmin xmax];
ax.YLim = [ymin ymax];
ax.YTick = [.9 .95 1];
l = legend(sprintf('AUC = %0.2f',AUC.rf),['\bf{}' sprintf('AUC = %0.2f',AUC.rerf)],sprintf('AUC = %0.2f',AUC.rf_rot));
l.Location = 'southeast';
l.Box = 'off';

if runSims
    run_Performance_profile_time_untransformed
else
    load Performance_profile_time_untransformed.mat
end

ax = subplot(4,2,2);
hold on

xmin = min(tau);
xmax = max(tau);
ymin = 0;
ymax = 1.17;

for i = 1:length(clnames)-2
     plot(tau,rho_ps(:,i),'LineWidth',LineWidth,'Color',Colors.(clnames{i}))
end

title('(E)','Units','normalized','Position',[0.025 1],'HorizontalAlignment','left','VerticalAlignment','top')
text(0.5,1,'Train Time','FontSize',14,'FontWeight','bold','Units','normalized','HorizontalAlignment','center','VerticalAlignment','bottom')
xlabel('Relative Performance')
ylabel('Proportion')
ax.LineWidth = LineWidth;
ax.FontUnits = 'inches';
ax.FontSize = FontSize;
ax.Units = 'inches';
ax.Position = [axLeft(5) axBottom(5) axWidth axHeight];
ax.Box = 'off';
ax.XLim = [xmin xmax];
 ax.YLim = [ymin ymax];
ax.YTick = 0:.5:1;
l = legend(['\bf{}' sprintf('AUC = %0.2f',AUC.rf)],sprintf('AUC = %0.2f',AUC.rerf),sprintf('AUC = %0.2f',AUC.rf_rot));
l.Units = 'inches';
l.Location = 'southeast';
% legLeft = l.Position(1);
% legBottom = l.Position(2);
% l.Position = [legLeft legBottom legWidth legHeight];
l.Box = 'off';


if runSims
    run_Performance_profile_time_rotate
else
    load Performance_profile_time_rotate.mat
end

ax = subplot(4,2,4);
hold on

xmin = min(tau);
xmax = max(tau);
ymin = 0;
ymax = 1.17;

for i = 1:length(clnames)-1
     plot(tau,rho_ps(:,i),'LineWidth',LineWidth,'Color',Colors.(clnames{i}))
end

title('(F)','Units','normalized','Position',[0.025 1],'HorizontalAlignment','left','VerticalAlignment','top')
ax.LineWidth = LineWidth;
ax.FontUnits = 'inches';
ax.FontSize = FontSize;
ax.Units = 'inches';
ax.Position = [axLeft(6) axBottom(6) axWidth axHeight];
ax.Box = 'off';
ax.XLim = [xmin xmax];
ax.YLim = [ymin ymax];
ax.YTick = 0:.5:1;
l = legend(['\bf{}' sprintf('AUC = %0.2f',AUC.rf)],sprintf('AUC = %0.2f',AUC.rerf),['\bf{}' sprintf('AUC = %0.2f',AUC.rf_rot)]);
l.Location = 'southeast';
l.Box = 'off';

if runSims
    run_Performance_profile_time_scale
else
    load Performance_profile_time_scale.mat
end

ax = subplot(4,2,6);
hold on

xmin = min(tau);
xmax = max(tau);
ymin = 0;
ymax = 1.17;

for i = 1:length(clnames)-1
     plot(tau,rho_ps(:,i),'LineWidth',LineWidth,'Color',Colors.(clnames{i}))
end

title('(G)','Units','normalized','Position',[0.025 1],'HorizontalAlignment','left','VerticalAlignment','top')
ax.LineWidth = LineWidth;
ax.FontUnits = 'inches';
ax.FontSize = FontSize;
ax.Units = 'inches';
ax.Position = [axLeft(7) axBottom(7) axWidth axHeight];
ax.Box = 'off';
ax.XLim = [xmin xmax];
ax.YLim = [ymin ymax];
ax.XTick = 1:4;
ax.YTick = 0:.5:1;
l = legend(['\bf{}' sprintf('AUC = %0.2f',AUC.rf)],sprintf('AUC = %0.2f',AUC.rerf),sprintf('AUC = %0.2f',AUC.rf_rot));
l.Location = 'southeast';
l.Box = 'off';


if runSims
    run_Performance_profile_time_affine
else
    load Performance_profile_time_affine.mat
end

ax = subplot(4,2,8);
hold on

xmin = min(tau);
xmax = max(tau);
ymin = 0;
ymax = 1.17;


for i = 1:length(clnames)-1
     plot(tau,rho_ps(:,i),'LineWidth',LineWidth,'Color',Colors.(clnames{i}))
end

title('(H)','Units','normalized','Position',[0.025 1],'HorizontalAlignment','left','VerticalAlignment','top')
ax.LineWidth = LineWidth;
ax.FontUnits = 'inches';
ax.FontSize = FontSize;
ax.Units = 'inches';
ax.Position = [axLeft(8) axBottom(8) axWidth axHeight];
ax.Box = 'off';
ax.XLim = [xmin xmax];
ax.YLim = [ymin ymax];
ax.YTick = 0:.5:1;
l = legend(sprintf('AUC = %0.2f',AUC.rf),['\bf{}' sprintf('AUC = %0.2f',AUC.rerf)],sprintf('AUC = %0.2f',AUC.rf_rot));
l.Location = 'southeast';
l.Box = 'off';

save_fig(gcf,[rerfPath 'RandomerForest/Figures/Fig4_benchmark2'])