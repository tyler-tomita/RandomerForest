%% Plot Sparse Parity and Trunk Transformations

clear
close all
clc

fpath = mfilename('fullpath');
frcPath = fpath(1:strfind(fpath,'RandomerForest')-1);

Colors.rf = 'c';
Colors.rfr = 'c';
Colors.frc = 'g';
Colors.frcr = 'g';
Colors.rr_rf = 'm';
Colors.rr_rfr = 'm';
LineStyles.rf = '-';
LineStyles.rfr = ':';
LineStyles.frc = '-';
LineStyles.frcr = ':';
LineStyles.rr_rf = '-';
LineStyles.rr_rfr = ':';
LineWidth = 2;
FontSize = .2;
axWidth = 2;
axHeight = 1.3;
axLeft = repmat([FontSize*4,FontSize*6.5+axWidth],1,4);
axBottom = [(FontSize*6+axHeight*3)*ones(1,2),...
    (FontSize*4.5+axHeight*2)*ones(1,2),(FontSize*3+axHeight)*ones(1,2),...
    FontSize*1.5*ones(1,2)];
legWidth = 0.4*axWidth;
legHeight = 0.4*axHeight;
legLeft = axLeft(end) + axWidth + FontSize;
legBottom = axBottom(5);
% figWidth = legLeft(end) + legWidth;
% figHeight = axBottom(1) + axHeight + FontSize*1.5;
figWidth = axLeft(end) + axWidth + FontSize;
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
    ax(2*i-1) = axes;
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
            ylabel(['\bf{' Transformations{i} '}'])
        text(0.5,1.05,'Sparse Parity','FontSize',16,'FontWeight','bold','Units',...
            'normalized','HorizontalAlignment','center','VerticalAlignment'...
            ,'bottom')
    else
        if strcmp(Transformations{i},'Outlier')
            ylabel('\bf{Corrupted}')
        else
            ylabel(['\bf{' Transformations{i} '}'])
        end
    end
    
    title(['(' char('A'+(2*(i-1))) ')'],'Units','normalized','Position',[0.025 .975],'HorizontalAlignment','left','VerticalAlignment','top')
    ax(2*i-1).LineWidth = LineWidth;
    ax(2*i-1).FontUnits = 'inches';
    ax(2*i-1).FontSize = FontSize;
    ax(2*i-1).Units = 'inches';
    ax(2*i-1).Position = [axLeft(2*i-1) axBottom(2*i-1) axWidth axHeight];
    ax(2*i-1).Box = 'off';
    ax(2*i-1).XLim = [1.5 45];
    ax(2*i-1).XScale = 'log';
    ax(2*i-1).XTick = [2 5 10 20 40];
    ax(2*i-1).XTickLabel = {'2' '5' '10' '20' '40'};
    ax(2*i-1).YLim = [0 .51];
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
    ax(2*i) = axes;
    for j = 1:length(Classifiers)
        errorbar(dims,mean(ErrorMatrix.(Classifiers{j}).(Transformations{i})),...
            std(ErrorMatrix.(Classifiers{j}).(Transformations{i}))/sqrt(ntrials),...
            'LineWidth',LineWidth,'Color',Colors.(Classifiers{j}),...
            'LineStyle',LineStyles.(Classifiers{j}));
        hold on
    end
    
    if i==1
        text(0.5,1.05,'Trunk','FontSize',16,'FontWeight','bold','Units',...
            'normalized','HorizontalAlignment','center','VerticalAlignment'...
            ,'bottom')
    end
    
    title(['(' char('A'+2*i-1) ')'],'Units','normalized','Position',[0.025 .975],'HorizontalAlignment','left','VerticalAlignment','top')
    ax(2*i).LineWidth = LineWidth;
    ax(2*i).FontUnits = 'inches';
    ax(2*i).FontSize = FontSize;
    ax(2*i).Units = 'inches';
    ax(2*i).Position = [axLeft(2*i) axBottom(2*i) axWidth axHeight];
    ax(2*i).Box = 'off';
    ax(2*i).XLim = [9 600];
    ax(2*i).XScale = 'log';
    ax(2*i).XTick = [10 100 500];
    ax(2*i).XTickLabel = {'10' '100' '500'};
    ax(2*i).YLim = [0.01 .15];
    
%     if i==length(Transformations)
%         l = legend('RF','RF(r)','F-RC','F-RC(r)','RR-RF','RR-RF(r)');
%         l.Box = 'off';
%         l.FontSize = 10;
%         l.Units = 'inches';
%         l.Position = [legLeft legBottom legWidth legHeight];
%     end
end

save_fig(gcf,[frcPath 'RandomerForest/Figures/Fig3_transformations'])