%% Plot Sparse Parity and Trunk Transformations

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
LineStyles.rerf = '-';
LineStyles.rerfr = ':';
LineStyles.frc = '-';
LineStyles.frcr = ':';
LineStyles.rr_rf = '-';
LineStyles.rr_rfr = ':';
LineWidth = 2;
FontSize = .2;
axWidth = 1.3;
axHeight = 1.3;
axLeft = repmat([FontSize*4,FontSize*6.5+axWidth,FontSize*9+axWidth*2,...
    FontSize*11.5+axWidth*3],1,2);
axBottom = [(FontSize*6+axHeight)*ones(1,4),FontSize*3*ones(1,4)];
% axLeft = repmat([FontSize*4,FontSize*6.5+axWidth],1,4);
% axBottom = [(FontSize*7+axHeight*3)*ones(1,2),...
%     (FontSize*4.5+axHeight*2)*ones(1,2),(FontSize*3+axHeight)*ones(1,2),...
%     FontSize*1.5*ones(1,2)];
legWidth = axWidth;
legHeight = axHeight;
legLeft = axLeft(end) + axWidth*2/3 + FontSize;
legBottom = axBottom(end);
figWidth = legLeft + legWidth + FontSize;
figHeight = axBottom(1) + axHeight + FontSize*1.5;
% figWidth = axLeft(end) + axWidth + FontSize;
% figHeight = axBottom(1) + axHeight + FontSize*1.5;

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

TestError = TestError(~cellfun(@isempty,TestError));

ntrials = size(TestError{1}.rf.Untransformed,1);

for i = 1:length(dims)
    Classifiers = fieldnames(TestError{i});
    Classifiers(~ismember(Classifiers,{'rf','rfr','frc','frcr','rr_rf','rr_rfr'})) = [];
    for j = 1:length(Classifiers)
        Transformations = fieldnames(TestError{i}.(Classifiers{j}));
        for k = 1:length(Transformations)
            ErrorMatrix.(Classifiers{j}).(Transformations{k})(:,i) = ...
                TestError{i}.(Classifiers{j}).(Transformations{k})';
        end
    end
end

Transformations(strcmp(Transformations,'Untransformed')) = [];

for i = 1:length(Transformations)
    ax(i) = axes;
    for j = 1:length(Classifiers)
        errorbar(dims,mean(ErrorMatrix.(Classifiers{j}).(Transformations{i})),...
            std(ErrorMatrix.(Classifiers{j}).(Transformations{i}))/sqrt(ntrials),...
            'LineWidth',LineWidth,'Color',Colors.(Classifiers{j}),...
            'LineStyle',LineStyles.(Classifiers{j}));
        hold on
    end
    
    if i==1
        xlabel('p')
%         ylabel({'\bf{Raw}';'\rm{Error Rate}'})
        ylabel({'\bf{Sparse Parity}';'\rm{Error Rate}'})
        text(0.5,1.05,Transformations{i},'FontSize',16,'FontWeight','bold','Units',...
            'normalized','HorizontalAlignment','center','VerticalAlignment'...
            ,'bottom')
    else
        if strcmp(Transformations{i},'Outlier')
            text(0.5,1.05,'Corrupted','FontSize',16,'FontWeight','bold','Units',...
                'normalized','HorizontalAlignment','center','VerticalAlignment'...
                ,'bottom')
        else
            text(0.5,1.05,Transformations{i},'FontSize',16,'FontWeight','bold','Units',...
                'normalized','HorizontalAlignment','center','VerticalAlignment'...
                ,'bottom')
        end
    end
    
%     title(['(' char('A'+(2*(i-1))) ')'],'Units','normalized','Position',[0.025 .975],'HorizontalAlignment','left','VerticalAlignment','top')
    ax(i).LineWidth = LineWidth;
    ax(i).FontUnits = 'inches';
    ax(i).FontSize = FontSize;
    ax(i).Units = 'inches';
    ax(i).Position = [axLeft(i) axBottom(i) axWidth axHeight];
    ax(i).Box = 'off';
    ax(i).XLim = [1.5 45];
    ax(i).XScale = 'log';
    ax(i).XTick = [2 5 10 20 40];
    ax(i).XTickLabel = {'2' '5' '10' '20' '40'};
    ax(i).YLim = [0 .51];
    ch = allchild(ax(i));
    uistack(ch(5),'top');
    uistack(ch(4),'top');
end

clear ErrorMatrix
%% Plot Trunk

runSims = false;

if runSims
    run_Trunk
else
    load Trunk_p_2_500
end

TestError = TestError(~cellfun(@isempty,TestError));

%plot only dimensions that have complete results
dims = dims(1:end-1);

ntrials = size(TestError{1}.rf.Untransformed,1);

for i = 1:length(dims)
    Classifiers = fieldnames(TestError{i});
    Classifiers(~ismember(Classifiers,{'rf','rfr','frc','frcr','rr_rf','rr_rfr'})) = [];
    for j = 1:length(Classifiers)
        Transformations = fieldnames(TestError{i}.(Classifiers{j}));
        for k = 1:length(Transformations)
            ErrorMatrix.(Classifiers{j}).(Transformations{k})(:,i) = ...
                TestError{i}.(Classifiers{j}).(Transformations{k})';
        end
    end
end

Transformations(strcmp(Transformations,'Untransformed')) = [];

for i = 1:length(Transformations)
    ax(i+4) = axes;
    for j = 1:length(Classifiers)
        errorbar(dims,mean(ErrorMatrix.(Classifiers{j}).(Transformations{i})),...
            std(ErrorMatrix.(Classifiers{j}).(Transformations{i}))/sqrt(ntrials),...
            'LineWidth',LineWidth,'Color',Colors.(Classifiers{j}),...
            'LineStyle',LineStyles.(Classifiers{j}));
        hold on
    end
    
    if i==1
        ylabel('\bf{Trunk}')
    end
    
%     title(['(' char('A'+i) ')'],'Units','normalized','Position',[0.025 .975],'HorizontalAlignment','left','VerticalAlignment','top')
    ax(i+4).LineWidth = LineWidth;
    ax(i+4).FontUnits = 'inches';
    ax(i+4).FontSize = FontSize;
    ax(i+4).Units = 'inches';
    ax(i+4).Position = [axLeft(i+4) axBottom(i+4) axWidth axHeight];
    ax(i+4).Box = 'off';
    ax(i+4).XLim = [9 600];
    ax(i+4).XScale = 'log';
    ax(i+4).XTick = [10 100 500];
    ax(i+4).XTickLabel = {'10' '100' '500'};
    ax(i+4).YLim = [0.01 .15];
    
    if i==length(Transformations)
        [lh,objh] = legend('RF','RF(r)','F-RC','Frank','RR-RF','RR-RF(r)');
        lh.Box = 'off';
        lh.FontSize = 10;
        lh.Units = 'inches';
        lh.Position = [legLeft legBottom legWidth legHeight];
        for j = length(objh)/2+1:length(objh)
            objh(j).Children.Children(2).XData = ...
                [(objh(j).Children.Children(2).XData(2)-objh(j).Children.Children(2).XData(1))*.75+objh(j).Children.Children(2).XData(1),objh(j).Children.Children(2).XData(2)];
        end
    end
    ch = allchild(ax(i+4));
    uistack(ch(5),'top');
    uistack(ch(4),'top');
end

save_fig(gcf,[rerfPath 'RandomerForest/Figures/ROFLMAO_fig4_transformations_poster_version'],{'fig','pdf','png'})