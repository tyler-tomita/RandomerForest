%% Plot benchmark classifier rank distributions

clear
close all
clc

LineWidth = 2;
MarkerSize = 6;
FontSize = .18;
axWidth = 8;
axHeight = 1.5;
axLeft = FontSize*5*ones(1,5);
axBottom = [FontSize*10+axHeight*4,FontSize*7+axHeight*3,...
    FontSize*5+axHeight*2,FontSize*3+axHeight,...
    FontSize];
figWidth = axLeft(end) + axWidth + FontSize;
figHeight = axBottom(1) + axHeight + FontSize*2;

fig = figure;
fig.Units = 'inches';
fig.PaperUnits = 'inches';
fig.Position = [0 0 figWidth figHeight];
fig.PaperPosition = [0 0 figWidth figHeight];
fig.PaperSize = [figWidth figHeight];

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

Transformations = {'Untransformed','Rotated','Scaled','Affine','Outlier'};

for i = 1:length(Transformations)
    load(['~/Benchmarks/Results/Benchmark_' lower(Transformations{i}) '.mat'])
    Classifiers = fieldnames(TestError{1});
    Classifiers(~ismember(Classifiers,{'rf','rfr','rerf','rerfr','frc','frcr','rr_rf','rr_rfr'})) = [];

    TestError = TestError(~cellfun(@isempty,TestError));
    
    RelativeError = {};

    for j = 1:length(TestError)
        for k = 2:length(Classifiers)
            RelativeError{k-1}(j,1) = TestError{j}.(Classifiers{k}) - TestError{j}.(Classifiers{1});
        end
    end
    
    ax(i) = axes;
    hold on

%     figure;
    
    h = plotSpread(RelativeError,[],[],{'RF(r)','RerF','RerF(r)','F-RC','Frank','RR-RF','RR-RF(r)'},...
        2);
    
%     h{3}.LineWidth = LineWidth;
%     h{3}.FontUnits = 'inches';
%     h{3}.FontSize = FontSize;
%     h{3}.Units = 'inches';
%     h{3}.Position = [axLeft(i) axBottom(i) axWidth axHeight];
%     h{3}.YScale = 'log';
    
    
    if i == 1
        ylabel('Relative Error');
        text(0.5,1.025,'Raw','FontSize',14,'FontWeight','bold','Units',...
            'normalized','HorizontalAlignment','center','VerticalAlignment'...
            ,'bottom')
    else
        if i == 5
            text(0.5,1.025,'Corrupted','FontSize',14,'FontWeight','bold','Units',...
                'normalized','HorizontalAlignment','center','VerticalAlignment'...
                ,'bottom')
        else
            text(0.5,1.025,Transformations{i},'FontSize',14,'FontWeight','bold','Units',...
                'normalized','HorizontalAlignment','center','VerticalAlignment'...
                ,'bottom')
        end
        ax(i).XTickLabel = {};
    end
    Mu = h{2}(1).YData;
    h{2}(1).Visible = 'off';
    h{3}.XLim = [0.5,7.5];
    h{3}.YLim = [-0.1,0.1];
%     h{3}.FontSize = 14;
    text(0,1.025,['(' char('A'+i-1) ')'],'FontSize',14,'FontWeight','bold','Units',...
        'normalized','HorizontalAlignment','left','VerticalAlignment'...
        ,'bottom')
%     title(['(' char('A'+i-1) ')'],'Units','normalized','Position',[0.025 .975],'HorizontalAlignment','left','VerticalAlignment','bottom')
    ax(i).LineWidth = LineWidth;
    ax(i).FontUnits = 'inches';
    ax(i).FontSize = FontSize;
    ax(i).Units = 'inches';
    ax(i).Position = [axLeft(i) axBottom(i) axWidth axHeight];
    
    h_line = allchild(h{3});
    h_line = flipud(h_line(end-6:end));
    
    for j = 1:length(h_line)
        h_line(j).Color = 'c';
        h_line(j).MarkerSize = MarkerSize;
        plot([min(h_line(j).XData),max(h_line(j).XData)],[Mu(j) Mu(j)],...
            'Color','m','LineWidth',LineWidth)
    end
    
    ColoredIdx = [1,3,5,7];
    for j = ColoredIdx
        p = patch([j-0.5 j+0.5 j+0.5 j-0.5],[h{3}.YLim(1) h{3}.YLim(1) h{3}.YLim(2) h{3}.YLim(2)],...
            [0.7 0.7 0.7]);
        p.EdgeColor = 'none';
    end
    
    ColoredIdx = [2,4,6];
    for j = ColoredIdx
        p = patch([j-0.5 j+0.5 j+0.5 j-0.5],[h{3}.YLim(1) h{3}.YLim(1) h{3}.YLim(2) h{3}.YLim(2)],...
            [0.75 0.75 0.75]);
        p.EdgeColor = 'none';
    end
    
    ch = h{3}.Children;
    ch(1:7) = [];
    ch(end+1:end+7) = h{3}.Children(1:7);
    h{3}.Children = ch;
    
%     for j = 1:length(RelativeError)
%         t = text(j,h{3}.YLim(2),...
%             sprintf('%0.2e +/-\n%0.2e',mean(RelativeError{j}),std(RelativeError{j})/sqrt(length(RelativeError{j}))),...
%             'HorizontalAlignment','center',...
%             'VerticalAlignment','top','FontSize',12,'Color','k');
%     end
end
save_fig(gcf,[rerfPath 'RandomerForest/Figures/Fig4_benchmark_beeswarm'])