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
FontSize = .2;
axWidth = 1.5;
axHeight = 1.5;
axLeft = repmat([FontSize*4,FontSize*8+axWidth],1,3);
axBottom = [...
    (FontSize*9+axHeight*2)*ones(1,2),(FontSize*6+axHeight)*ones(1,2),...
    FontSize*3*ones(1,2)];
legWidth = axWidth;
legHeight = [axHeight/2,axHeight];
legLeft = axLeft(end) + axWidth*2/3 + FontSize;
legBottom = [axBottom(1) + axHeight/2 - legHeight(1)/2,axBottom(end)];
figWidth = legLeft + legWidth + FontSize;
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

plot(X(Y==0,1),X(Y==0,2),'.','Color',Colors.Class0,'MarkerSize',MarkerSize)
hold on
plot(X(Y==1,1),X(Y==1,2),'.','Color',Colors.Class1,'MarkerSize',MarkerSize)

title('(A)','Units','normalized','Position',[0.025 .975],'HorizontalAlignment','left','VerticalAlignment','top')
text(0.5,1.05,'Sparse Parity','FontSize',16,'FontWeight','bold','Units',...
    'normalized','HorizontalAlignment','center','VerticalAlignment'...
    ,'bottom')
xlabel('X_1')
ylabel('X_2')
[lh1,objh1] = legend('Class 0','Class 1');
lh1.Units = 'inches';
lh1.Position = [legLeft legBottom(1) legWidth legHeight(1)];
lh1.Box = 'off';
lh1.FontSize = 11;
objh1(4).XData = 0.75*(objh1(3).XData(2) - objh1(3).XData(1)) + objh1(3).XData(1);
objh1(6).XData = 0.75*(objh1(5).XData(2) - objh1(5).XData(1)) + objh1(5).XData(1);
ax.LineWidth = LineWidth;
ax.FontUnits = 'inches';
ax.FontSize = FontSize;
ax.Units = 'inches';
ax.Position = [axLeft(1) axBottom(1) axWidth axHeight];
ax.Box = 'off';
ax.XLim = [-2 2];
ax.YLim = [-2 2];
ax.XTick = [-2 0 2];
ax.YTick = [-2 0 2];
ax.XTickLabel = {'-2';'0';'2'};
ax.YTickLabel = {'-2';'0';'2'};

runSims = false;

if runSims
    run_Sparse_parity
else
    load Sparse_parity
end

ax = axes;

for i = 1:length(TestError)
    Classifiers = fieldnames(TestError{i});
    Classifiers(~ismember(Classifiers,{'rf','rfr','frc','frcr','rr_rf','rr_rfr'})) = [];
    for j = 1:length(Classifiers)
        ntrials = length(TestError{i}.(Classifiers{j}).Untransformed);
        for trial = 1:ntrials
            BestIdx = hp_optimize(OOBError{i}.(Classifiers{j}).Untransformed(trial,:),...
                OOBAUC{i}.(Classifiers{j}).Untransformed(trial,:));
            if length(BestIdx) > 1
                BestIdx = BestIdx(end);
            end
            PlotTime.(Classifiers{j})(trial,i) = ...
                TrainTime{i}.(Classifiers{j}).Untransformed(trial,BestIdx);
            PlotError.(Classifiers{j})(trial,i) = ...
                TestError{i}.(Classifiers{j}).Untransformed(trial);
        end
    end
end

for i = 1:length(Classifiers)
    cl = Classifiers{i};
    hTestError(i) = errorbar(dims,mean(PlotError.(cl)),...
        std(PlotError.(cl))/sqrt(size(PlotError.(cl),1)),...
        'LineWidth',LineWidth,'Color',Colors.(cl),...
        'LineStyle',LineStyles.(cl));
    hold on
end

title('(C)','Units','normalized','Position',[0.025 .975],'HorizontalAlignment','left','VerticalAlignment','top')
% text(0.5,1.05,{'Error Rate';'(relative to RF)'},'FontSize',16,'FontWeight','bold','Units','normalized','HorizontalAlignment','center','VerticalAlignment','bottom')
xlabel('p')
ylabel('Error Rate')
ax.LineWidth = LineWidth;
ax.FontUnits = 'inches';
ax.FontSize = FontSize;
ax.Units = 'inches';
ax.Position = [axLeft(3) axBottom(3) axWidth axHeight];
ax.Box = 'off';
ax.XLim = [1.5 45];
ax.XScale = 'log';
ax.XTick = [2 5 10 20 40];
ax.XTickLabel = {'2' '5' '10' '20' '40'};
ax.YLim = [0 .5];

ax = axes;

for i = 1:length(Classifiers)
    cl = Classifiers{i};
    hTrainTime(i) = errorbar(dims,mean(PlotTime.(cl)),...
        std(PlotTime.(cl))/sqrt(size(PlotTime.(cl),1)),...
        'LineWidth',LineWidth,'Color',Colors.(cl),...
        'LineStyle',LineStyles.(cl));
    hold on
end

title('(E)','Units','normalized','Position',[0.025 .975],'HorizontalAlignment','left','VerticalAlignment','top')
% text(0.5,1.05,'Training Time','FontSize',16,'FontWeight','bold','Units','normalized','HorizontalAlignment','center','VerticalAlignment','bottom')
xlabel('p')
ylabel('Train Time (s)')
ax.LineWidth = LineWidth;
ax.FontUnits = 'inches';
ax.FontSize = FontSize;
ax.Units = 'inches';
ax.Position = [axLeft(5) axBottom(5) axWidth axHeight];
ax.Box = 'off';
ax.XLim = [1.5 45];
ax.YLim = [1 150];
ax.XScale = 'log';
ax.YScale = 'log';
ax.XTick = [2 5 10 20 40];
ax.XTickLabel = {'2' '5' '10' '20' '40'};
ax.YTickLabel = {'1','10','100'};

clear hTestError hTrainTime TestError minTestError trainTime
%% Plot Trunk

ax = axes;

n = 100;
p = 2;

d_idx = 1:p;
mu1 = 1./sqrt(d_idx);
mu0 = -1*mu1;
Mu = cat(1,mu0,mu1);
Sigma = ones(1,p);
obj = gmdistribution(Mu,Sigma);
[X,idx] = random(obj,n);
Class = [0;1];
Y = Class(idx);

plot(X(Y==0,1),X(Y==0,2),'.','Color',Colors.Class0,'MarkerSize',MarkerSize)
hold on
plot(X(Y==1,1),X(Y==1,2),'.','Color',Colors.Class1,'MarkerSize',MarkerSize)

title('(B)','Units','normalized','Position',[0.025 .975],'HorizontalAlignment','left','VerticalAlignment','top')
text(0.5,1.05,'Trunk','FontSize',16,'FontWeight','bold','Units',...
    'normalized','HorizontalAlignment','center','VerticalAlignment'...
    ,'bottom')
xlabel('X_1')
ylabel('X_2')
ax.LineWidth = LineWidth;
ax.FontUnits = 'inches';
ax.FontSize = FontSize;
ax.Units = 'inches';
ax.Position = [axLeft(2) axBottom(2) axWidth axHeight];
ax.Box = 'off';
ax.XLim = [-5 5];
ax.YLim = [-5 5];
ax.XTick = [-5 0 5];
%ax.YTick = [-5 0 5];
ax.XTickLabel = {'-5';'0';'5'};
%ax.YTickLabel = {'-5';'0';'5'};

if runSims
    run_Trunk
else
    load Trunk_p_2_500
end

ax = axes;

clear PlotError PlotTime

for i = 1:length(TestError)
    Classifiers = fieldnames(TestError{i});
    Classifiers(~ismember(Classifiers,{'rf','rfr','frc','frcr','rr_rf','rr_rfr'})) = [];
    for j = 1:length(Classifiers)
        ntrials = length(TestError{i}.(Classifiers{j}).Untransformed);
        for trial = 1:ntrials
            BestIdx = hp_optimize(OOBError{i}.(Classifiers{j}).Untransformed(trial,:),...
                OOBAUC{i}.(Classifiers{j}).Untransformed(trial,:));
            if length(BestIdx) > 1
                BestIdx = BestIdx(end);
            end
            PlotTime.(Classifiers{j})(trial,i) = ...
                TrainTime{i}.(Classifiers{j}).Untransformed(trial,BestIdx);
            PlotError.(Classifiers{j})(trial,i) = ...
                TestError{i}.(Classifiers{j}).Untransformed(trial);
        end
    end
end

for i = 1:length(Classifiers)
    cl = Classifiers{i};
    hTestError(i) = errorbar(dims(dims<=500),mean(PlotError.(cl)),...
        std(PlotError.(cl))/sqrt(size(PlotError.(cl),1)),...
        'LineWidth',LineWidth,'Color',Colors.(cl),...
        'LineStyle',LineStyles.(cl));
    hold on
end

title('(D)','Units','normalized','Position',[0.025 .975],'HorizontalAlignment','left','VerticalAlignment','top')
xlabel('p')
ylabel('Error Rate')
ax.LineWidth = LineWidth;
ax.FontUnits = 'inches';
ax.FontSize = FontSize;
ax.Units = 'inches';
ax.Position = [axLeft(4) axBottom(4) axWidth axHeight];
ax.Box = 'off';
ax.XLim = [9 600];
ax.YLim = [0.02 0.15];
ax.XScale = 'log';
ax.XTick = [10,100,500];
ax.XTickLabel = {'10','100','500'};

ax = axes;

for i = 1:length(Classifiers)
    cl = Classifiers{i};
    hTrainTime(i) = errorbar(dims(dims<=500),mean(PlotTime.(cl)),...
        std(PlotTime.(cl))/sqrt(size(PlotTime.(cl),1)),...
        'LineWidth',LineWidth,'Color',Colors.(cl),...
        'LineStyle',LineStyles.(cl));
    hold on
end

title('(F)','Units','normalized','Position',[0.025 .975],'HorizontalAlignment','left','VerticalAlignment','top')
xlabel('p')
ylabel('Train Time (s)')
ax.LineWidth = LineWidth;
ax.FontUnits = 'inches';
ax.FontSize = FontSize;
ax.Units = 'inches';
ax.Position = [axLeft(6) axBottom(6) axWidth axHeight];
ax.Box = 'off';
ax.XLim = [9 600];
ax.YLim = [0 150];
ax.XScale = 'log';
ax.YScale = 'log';
ax.XTick = [10,100,500];
ax.YTick = [1,10,100];
ax.XTickLabel = {'10','100','500'};
ax.YTickLabel = {'1','10','100'};
[lh2,objh2] = legend('RF','RF(r)','F-RC','Frank','RR-RF','RR-RF(r)');
% lh2.Location = 'southwest';
lh2.Units = 'inches';
lh2.Position = [legLeft legBottom(2) legWidth legHeight(2)];
lh2.Box = 'off';
lh2.FontSize = 11;

for i = 7:length(objh2)
    objh2(i).Children.Children(2).XData = [(objh2(i).Children.Children(2).XData(2)-objh2(i).Children.Children(2).XData(1))*.75+objh2(i).Children.Children(2).XData(1),objh2(i).Children.Children(2).XData(2)];
end

objh1(1).FontSize = objh2(1).FontSize;
objh1(2).FontSize = objh2(1).FontSize;

save_fig(gcf,[rerfPath 'RandomerForest/Figures/ROFLMAO_fig3_simulations_2017_01_23'],{'fig','pdf','png'})