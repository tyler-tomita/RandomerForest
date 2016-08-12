close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

rng(1);

C = [0 1 1;0 1 0;1 0 1;1 0 0;0 0 0;1 .5 0];
Colors.rf = C(1,:);
Colors.rerf = C(2,:);
Colors.rr_rf = C(3,:);
Colors.rerfr = C(4,:);
Colors.frc = C(5,:);
% Colors.rerfd = C(6,:);
Colors.rerfu = C(6,:);
LineWidth = 2;
MarkerSize = 12;
FontSize = .2;
TitleFontSize = 18;
axWidth = 1.75;
axHeight = 1.75;
axLeft = [FontSize*4 FontSize*8+axWidth FontSize*12+axWidth*2 FontSize*4,...
    FontSize*8+axWidth FontSize*12+axWidth*2];
axBottom = [FontSize*7+axHeight FontSize*7+axHeight FontSize*7+axHeight,...
    FontSize*4 FontSize*4 FontSize*4];
figWidth = axLeft(end) + axWidth + FontSize*4;
figHeight = axBottom(1) + axHeight + FontSize*2;

fig = figure;
fig.Units = 'inches';
fig.PaperUnits = 'inches';
fig.Position = [0 0 figWidth figHeight];
fig.PaperPosition = [0 0 figWidth figHeight];
fig.PaperSize = [figWidth figHeight];

ax = subplot(2,3,1);

n = 1000;
d = 2;
dgood = min(3,d);
Sigma = 1/32*ones(1,d);
Mu = sparse(n,d);

for jj = 1:n
    Mu(jj,:) = binornd(1,0.5,1,d);
    X(jj,1:d) = mvnrnd(Mu(jj,:),Sigma);
end

nones = sum(Mu(:,1:dgood),2);
Y = mod(nones,2);

plot(X(Y==0,1),X(Y==0,2),'b.',X(Y==1,1),X(Y==1,2),'r.','MarkerSize',MarkerSize)

title('(A)','Units','normalized','Position',[0.025 .975],'HorizontalAlignment','left','VerticalAlignment','top')
% text(0.5,1.05,'Scatterplot','FontSize',16,'FontWeight','bold','Units','normalized','HorizontalAlignment','center','VerticalAlignment','bottom')
xlabel('X_1')
ylabel({'\bf{Sparse Parity}';'\rm{X_2}'})
l = legend('Class 1','Class 2');
l.Location = 'southwest';
l.Box = 'off';
l.FontSize = 12;
ax.LineWidth = LineWidth;
ax.FontUnits = 'inches';
ax.FontSize = FontSize;
ax.Units = 'inches';
ax.Position = [axLeft(1) axBottom(1) axWidth axHeight];
ax.Box = 'off';
ax.XLim = [-2 3];
ax.YLim = [-2 3];
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

ax = subplot(2,3,2);

for i = 1:length(TestError)
    Classifiers = fieldnames(TestError{i});
    for j = 1:length(Classifiers)
        ntrials = length(TestError{i}.(Classifiers{j}));
        for trial = 1:ntrials
            BestIdx = hp_optimize(OOBError{i}.(Classifiers{j})(trial,:),...
                OOBAUC{i}.(Classifiers{j})(trial,:));
            if length(BestIdx) > 1
                BestIdx = BestIdx(end);
            end
            PlotTime.(Classifiers{j})(trial,i) = ...
                TrainTime{i}.(Classifiers{j})(trial,BestIdx);
            PlotError.(Classifiers{j})(trial,i) = ...
                TestError{i}.(Classifiers{j})(trial);
        end
    end
end

Classifiers = {'rerf' 'rerfu' 'rf' 'frc' 'rr_rf'};

for i = 1:length(Classifiers)
    cl = Classifiers{i};
    if ~strcmp(cl,'rerfd')
        hTestError(i) = errorbar(dims,mean(PlotError.(cl)),...
            std(PlotError.(cl))/sqrt(size(PlotError.(cl),1)),...
            'LineWidth',LineWidth,'Color',Colors.(cl));
        hold on
    end
end

title('(B)','Units','normalized','Position',[0.025 .975],'HorizontalAlignment','left','VerticalAlignment','top')
% text(0.5,1.05,{'Error Rate';'(relative to RF)'},'FontSize',16,'FontWeight','bold','Units','normalized','HorizontalAlignment','center','VerticalAlignment','bottom')
xlabel('p')
ylabel('Relative Error')
ax.LineWidth = LineWidth;
ax.FontUnits = 'inches';
ax.FontSize = FontSize;
ax.Units = 'inches';
ax.Position = [axLeft(2) axBottom(2) axWidth axHeight];
ax.Box = 'off';
ax.XLim = [1 105];
ax.XScale = 'log';
ax.XTick = [1 10 100];
ax.XTickLabel = {'1' '10' '100'};
ax.YLim = [0 .5];

ax = subplot(2,3,3);

for i = 1:length(Classifiers)
    cl = Classifiers{i};
    if ~strcmp(cl,'rerfd')
        hTrainTime(i) = errorbar(dims,mean(PlotTime.(cl)),...
            std(PlotTime.(cl))/sqrt(size(PlotTime.(cl),1)),...
            'LineWidth',LineWidth,'Color',Colors.(cl));
        hold on
    end
end

title('(C)','Units','normalized','Position',[0.025 .975],'HorizontalAlignment','left','VerticalAlignment','top')
% text(0.5,1.05,'Training Time','FontSize',16,'FontWeight','bold','Units','normalized','HorizontalAlignment','center','VerticalAlignment','bottom')
xlabel('p')
ylabel('Train Time (s)')
ax.LineWidth = LineWidth;
ax.FontUnits = 'inches';
ax.FontSize = FontSize;
ax.Units = 'inches';
ax.Position = [axLeft(3) axBottom(3) axWidth axHeight];
ax.Box = 'off';
ax.XLim = [1 105];
ax.XTick = [1 10 100];
ax.XTickLabel = {'1' '10' '100'};
ax.XScale = 'log';

clear hTestError hTrainTime TestError minTestError trainTime

ax = subplot(2,3,4);

n = 100;
d = 2;

d_idx = 1:d;
mu1 = 1./sqrt(d_idx);
mu0 = -1*mu1;
Mu = cat(1,mu0,mu1);
Sigma = ones(1,d);
obj = gmdistribution(Mu,Sigma);
[X,idx] = random(obj,n);
Class = [0;1];
Y = Class(idx);

plot(X(Y==0,1),X(Y==0,2),'b.',X(Y==1,1),X(Y==1,2),'r.','MarkerSize',MarkerSize)

title('(D)','Units','normalized','Position',[0.025 .975],'HorizontalAlignment','left','VerticalAlignment','top')
xlabel('X_1')
ylabel({'\bf{Trunk}';'\rm{X_2}'})
ax.LineWidth = LineWidth;
ax.FontUnits = 'inches';
ax.FontSize = FontSize;
ax.Units = 'inches';
ax.Position = [axLeft(4) axBottom(4) axWidth axHeight];
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
    load Trunk
end

ax = subplot(2,3,5);

clear PlotError PlotTime

for i = 1:length(TestError)
    Classifiers = fieldnames(TestError{i});
    for j = 1:length(Classifiers)
        ntrials = length(TestError{i}.(Classifiers{j}));
        for trial = 1:ntrials
            BestIdx = hp_optimize(OOBError{i}.(Classifiers{j})(trial,:),...
                OOBAUC{i}.(Classifiers{j})(trial,:));
            if length(BestIdx) > 1
                BestIdx = BestIdx(end);
            end
            PlotTime.(Classifiers{j})(trial,i) = ...
                TrainTime{i}.(Classifiers{j})(trial,BestIdx);
            PlotError.(Classifiers{j})(trial,i) = ...
                TestError{i}.(Classifiers{j})(trial);
        end
    end
end

Classifiers = {'rerf' 'rerfu' 'rf' 'frc' 'rr_rf'};

for i = 1:length(Classifiers)
    cl = Classifiers{i};
    if ~strcmp(cl,'rerfd')
        if strcmp(cl,'rr_rf') || strcmp(cl,'rerfu')
            hTestError(i) = errorbar(dims(1:end-1),mean(PlotError.(cl)),...
                std(PlotError.(cl))/sqrt(size(PlotError.(cl),1)),...
                'LineWidth',LineWidth,'Color',Colors.(cl));
        else
            hTestError(i) = errorbar(dims,mean(PlotError.(cl)),...
                std(PlotError.(cl))/sqrt(size(PlotError.(cl),1)),...
                'LineWidth',LineWidth,'Color',Colors.(cl));
        end
        hold on
    end
end

title('(E)','Units','normalized','Position',[0.025 .975],'HorizontalAlignment','left','VerticalAlignment','top')
xlabel('p')
ylabel('Relative Error')
ax.LineWidth = LineWidth;
ax.FontUnits = 'inches';
ax.FontSize = FontSize;
ax.Units = 'inches';
ax.Position = [axLeft(5) axBottom(5) axWidth axHeight];
ax.Box = 'off';
ax.XLim = [9 1100];
ax.YLim = [0.02 0.1];
ax.XScale = 'log';
ax.XTick = logspace(0,3,4);
ax.XTickLabel = {'1';'10';'100';'1000'};

ax = subplot(2,3,6);

for i = 1:length(Classifiers)
    cl = Classifiers{i};
    if ~strcmp(cl,'rerfd')
        if strcmp(cl,'rr_rf') || strcmp(cl,'rerfu')
            hTrainTime(i) = errorbar(dims(1:end-1),mean(PlotTime.(cl)),...
                std(PlotTime.(cl))/sqrt(size(PlotTime.(cl),1)),...
                'LineWidth',LineWidth,'Color',Colors.(cl));
        else
            hTrainTime(i) = errorbar(dims,mean(PlotTime.(cl)),...
                std(PlotTime.(cl))/sqrt(size(PlotTime.(cl),1)),...
                'LineWidth',LineWidth,'Color',Colors.(cl));
        end
        hold on
    end
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
ax.XLim = [9 1100];
ax.XScale = 'log';
ax.XTick = logspace(0,3,4);
ax.XTickLabel = {'1';'10';'100';'1000'};
l = legend('RerF','RerF(u)','RF','F-RC','RR-RF');
l.Location = 'southwest';
l.Box = 'off';
l.FontSize = 12;

save_fig(gcf,[rerfPath 'RandomerForest/Figures/Fig2_simulations'])