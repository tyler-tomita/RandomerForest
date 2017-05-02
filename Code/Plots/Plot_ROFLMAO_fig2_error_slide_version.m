%% Plot Sparse Parity and Trunk

clear
close all
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

Colors.rf = 0.5*ones(1,3);
Colors.rfr = 0.5*ones(1,3);
Colors.frc= [0.2 0.2 0.2];
Colors.frcr = [0 0 0];
Colors.rr_rf = 0.8*ones(1,3);
Colors.rr_rfr = 0.8*ones(1,3);
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
LineWidth.rf = 2;
LineWidth.rfr = 2;
LineWidth.frc = 2;
LineWidth.frcr = 4;
LineWidth.rr_rf = 2;
LineWidth.rr_rfr = 2;
LineWidth.ax = 2;
MarkerSize = 12;
FontSize = .25;
axWidth = 1.5;
axHeight = 1.5;
axLeft = [FontSize*3.5,FontSize*8+axWidth];
axBottom = FontSize*3*ones(1,2);
legWidth = axWidth;
legHeight = axHeight;
legLeft = axLeft(2) + axWidth*3/4 + FontSize;
legBottom = axBottom(2) + axHeight/2 - legHeight/2;
figWidth = legLeft(end) + legWidth + FontSize/2;
figHeight = axBottom(1) + axHeight + FontSize*1.5;

fig = figure;
fig.Units = 'inches';
fig.PaperUnits = 'inches';
fig.Position = [0 0 figWidth figHeight];
fig.PaperPosition = [0 0 figWidth figHeight];
fig.PaperSize = [figWidth figHeight];

%% Plot Sparse Parity

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
        'LineWidth',LineWidth.(cl),'Color',Colors.(cl),...
        'LineStyle',LineStyles.(cl));
    hold on
end

text(0.5,1.05,'Sparse Parity','FontSize',16,'FontWeight','bold','Units',...
    'normalized','HorizontalAlignment','center','VerticalAlignment'...
    ,'bottom')
xlabel('p')
ylabel('Error Rate')
ax.LineWidth = LineWidth.ax;
ax.FontUnits = 'inches';
ax.FontSize = FontSize;
ax.Units = 'inches';
ax.Position = [axLeft(1) axBottom(1) axWidth axHeight];
ax.Box = 'off';
ax.XLim = [1.5 45];
ax.XScale = 'log';
ax.XTick = [2 5 10 20 40];
ax.XTickLabel = {'2' '5' '10' '20' '40'};
ax.YLim = [0 .5];

ch = allchild(ax);
uistack(ch(5),'top');
uistack(ch(4),'top');

clear hTestError hTrainTime TestError minTestError trainTime
%% Plot Trunk

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
        'LineWidth',LineWidth.(cl),'Color',Colors.(cl),...
        'LineStyle',LineStyles.(cl));
    hold on
end

text(0.5,1.05,'Trunk','FontSize',16,'FontWeight','bold','Units',...
    'normalized','HorizontalAlignment','center','VerticalAlignment'...
    ,'bottom')
xlabel('p')
ylabel('Error Rate')
ax.LineWidth = LineWidth.ax;
ax.FontUnits = 'inches';
ax.FontSize = FontSize;
ax.Units = 'inches';
ax.Position = [axLeft(2) axBottom(2) axWidth axHeight];
ax.Box = 'off';
ax.XLim = [9 600];
ax.YLim = [0.02 0.15];
ax.XScale = 'log';
ax.XTick = [10,100,500];
ax.XTickLabel = {'10','100','500'};

[lh2,objh2] = legend('RF','RF(r)','F-RC','Frank','RR-RF','RR-RF(r)');
% lh2.Location = 'southwest';
lh2.Units = 'inches';
lh2.Position = [legLeft legBottom legWidth legHeight];
lh2.Box = 'off';
lh2.FontSize = 11;

for i = 7:length(objh2)
    objh2(i).Children.Children(2).XData = [(objh2(i).Children.Children(2).XData(2)-objh2(i).Children.Children(2).XData(1))*.75+objh2(i).Children.Children(2).XData(1),objh2(i).Children.Children(2).XData(2)];
end

ch = allchild(ax);
uistack(ch(5),'top');
uistack(ch(4),'top');

save_fig(gcf,[rerfPath 'RandomerForest/Figures/ROFLMAO_fig2_error_slide_version'],{'fig','pdf','png'})