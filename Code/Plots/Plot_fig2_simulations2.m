%% Plot Sparse Parity and Trunk

clear
close all
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

C = [0 1 1;0 1 0;1 0 1;1 0 0;0 0 0;1 .5 0];
Colors.rf = C(1,:);
Colors.rfr = C(2,:);
Colors.frc = C(3,:);
Colors.frcr = C(4,:);
Colors.rr_rf = C(5,:);
Colors.rr_rfr = C(6,:);
LineWidth = 2;
MarkerSize = 12;
FontSize = .2;
axWidth = 1.3;
axHeight = 1.3;
axLeft = repmat([FontSize*4,FontSize*8+axWidth],1,3);
axBottom = [...
    (FontSize*9+axHeight*2)*ones(1,2),(FontSize*6+axHeight)*ones(1,2),...
    FontSize*3*ones(1,2)];
legWidth = axWidth;
legHeight = axHeight;
legLeft = axLeft(end) + axWidth + FontSize;
legBottom = axBottom(5);
figWidth = axLeft(end) + axWidth + FontSize;
figHeight = axBottom(1) + axHeight + FontSize*1.5;

fig = figure;
fig.Units = 'inches';
fig.PaperUnits = 'inches';
fig.Position = [0 0 figWidth figHeight];
fig.PaperPosition = [0 0 figWidth figHeight];
fig.PaperSize = [figWidth figHeight];

%% Plot Sparse Parity

ax = axes;

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
text(0.5,1.05,'Sparse Parity','FontSize',14,'FontWeight','bold','Units',...
    'normalized','HorizontalAlignment','center','VerticalAlignment'...
    ,'bottom')
xlabel('X_1')
ylabel('X_2')
l = legend('Class 1','Class 2');
l.Location = 'southwest';
l.Box = 'off';
l.FontSize = 10;
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
        'LineWidth',LineWidth,'Color',Colors.(cl));
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
        'LineWidth',LineWidth,'Color',Colors.(cl));
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
ax.XScale = 'log';
ax.XTick = [2 5 10 20 40];
ax.XTickLabel = {'2' '5' '10' '20' '40'};
ax.XScale = 'log';

clear hTestError hTrainTime TestError minTestError trainTime
%% Plot Trunk

ax = axes;

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

title('(B)','Units','normalized','Position',[0.025 .975],'HorizontalAlignment','left','VerticalAlignment','top')
text(0.5,1.05,'Trunk','FontSize',14,'FontWeight','bold','Units',...
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
        'LineWidth',LineWidth,'Color',Colors.(cl));
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
ax.XLim = [9 1100];
ax.YLim = [0.02 0.1];
ax.XScale = 'log';
ax.XTick = logspace(0,3,4);
ax.XTickLabel = {'1';'10';'100';'1000'};

ax = axes;

for i = 1:length(Classifiers)
    cl = Classifiers{i};
    hTrainTime(i) = errorbar(dims(dims<=500),mean(PlotTime.(cl)),...
        std(PlotTime.(cl))/sqrt(size(PlotTime.(cl),1)),...
        'LineWidth',LineWidth,'Color',Colors.(cl));
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
ax.XLim = [9 1100];
ax.YLim = [0 100];
ax.XScale = 'log';
ax.XTick = logspace(0,3,4);
ax.XTickLabel = {'1';'10';'100';'1000'};
l = legend('RF','RF(r)','F-RC','F-RC(r)','RR-RF','RR-RF(r)');
l.Location = 'southwest';
l.Box = 'off';
l.FontSize = 10;

save_fig(gcf,[rerfPath 'RandomerForest/Figures/Fig2_simulations'])