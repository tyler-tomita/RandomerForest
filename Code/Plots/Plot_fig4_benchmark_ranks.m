%% Plot benchmark classifier rank distributions

clear
close all
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

LineWidth = 2;
FontSize = .12;
axWidth = 3.75;
axHeight = 0.75;
axLeft = FontSize*3*ones(1,5);
axBottom = [FontSize*16+axHeight*4,FontSize*11+axHeight*3,...
    FontSize*8+axHeight*2,FontSize*5+axHeight,...
    FontSize];
legWidth = 1;
legHeight = 1;
legLeft = axLeft(end) + axWidth - 4*FontSize;
legBottom = axBottom(3) - (legHeight - axHeight)/2;
figWidth = legLeft + legWidth;
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

    ErrorMatrix = [];
    for j = 1:length(TestError)
        for k = 1:length(Classifiers)
            ErrorMatrix(j,k) = TestError{j}.(Classifiers{k});
        end
    end

    ClRanks = tiedrank(ErrorMatrix')';
    IntRanks = floor(ClRanks);
    

    RankCounts = NaN(length(Classifiers));
    for j = 1:length(Classifiers)
        RankCounts(j,:) = sum(IntRanks==j); 
        MeanRank.(Classifiers{j}).(Transformations{i}) = mean(IntRanks(:,j));
    end
    
    ax(i) = axes;
    bar(RankCounts')
    Bars = allchild(ax(i));
    for j = 1:length(Bars)
        Bars(j).EdgeColor = 'k';
        Bars(j).BarWidth = 1;
    end

    ylabel('Frequency')
    if strcmp(Transformations{i},'Untransformed')
        text(0.5,1.05,'Raw','FontSize',12,'FontWeight','bold','Units',...
            'normalized','HorizontalAlignment','center','VerticalAlignment'...
            ,'bottom')
    elseif strcmp(Transformations{i},'Outlier')
        text(0.5,1.05,'Corrupted','FontSize',12,'FontWeight','bold','Units',...
            'normalized','HorizontalAlignment','center','VerticalAlignment'...
            ,'bottom')
        [lh,objh] = legend('1st','2nd','3rd','4th','5th','6th','7th','8th');
        lh.Box = 'off';
        lh.FontSize = 9;
        lh.Units = 'inches';
        lh.Position = [legLeft legBottom legWidth legHeight];
        BarWidth = (objh(9).Children.YData(2) - objh(9).Children.YData(1))/2;
        for j = 9:16
            objh(j).Children.XData(1) = objh(j).Children.XData(3) - BarWidth;
            objh(j).Children.XData(2) = objh(j).Children.XData(3) - BarWidth;
            objh(j).Children.XData(5) = objh(j).Children.XData(3) - BarWidth;
%             objh(j).Children.EdgeAlpha = 0;
        end
    else
        text(0.5,1.05,Transformations{i},'FontSize',12,'FontWeight','bold','Units',...
            'normalized','HorizontalAlignment','center','VerticalAlignment'...
            ,'bottom')
    end

    title(['(' char('A'+i-1) ')'],'Units','normalized','Position',[0.025 .975],'HorizontalAlignment','left','VerticalAlignment','top')
    ax(i).LineWidth = LineWidth;
    ax(i).FontUnits = 'inches';
    ax(i).FontSize = FontSize;
    ax(i).Units = 'inches';
    ax(i).Position = [axLeft(i) axBottom(i) axWidth axHeight];
    if i==1
        xlabel('Rank')
        ax(i).XTickLabel = {'RF' 'RF(r)' 'RerF' 'RerF(r)' 'F-RC' 'Frank' 'RR-RF' 'RR-RF(r)'};
    else
        ax(i).XTick = [];
    end
    ax(i).XLim = [0.5 8.5];
    ax(i).YLim = [0 60];
    
    ColoredIdx = [1,3,5,7];
    for j = ColoredIdx
        p = patch([j-0.5 j+0.5 j+0.5 j-0.5],[0 0 ax(i).YLim(2) ax(i).YLim(2)],...
            [0.9 0.9 0.9]);
        p.EdgeColor = 'none';
    end
    
    ColoredIdx = [2,4,6,8];
    for j = ColoredIdx
        p = patch([j-0.5 j+0.5 j+0.5 j-0.5],[0 0 ax(i).YLim(2) ax(i).YLim(2)],...
            [0.8 0.8 0.8]);
        p.EdgeColor = 'none';
    end
        
    ch = ax(i).Children;
    ch(1:8) = [];
    ch(end+1:end+8) = ax(i).Children(1:8);
    ax(i).Children = ch;

%     if i==5
%         [l,obj,~,~] = legend('1st place','2nd place','3rd place','4th place','5th place',...
%             '6th place');
%         l.Location = 'northwest';
%         l.Box = 'off';
%     end
    
    if strcmp(Transformations{i},'Untransformed')
        t = text(3,ax(i).YLim(2),'\bf{+}','HorizontalAlignment','center',...
            'VerticalAlignment','top','FontSize',14,'Color','r');
        t = text(4,ax(i).YLim(2),'\bf{+}','HorizontalAlignment','center',...
            'VerticalAlignment','top','FontSize',14,'Color','r');
        t = text(5,ax(i).YLim(2),'\bf{+}','HorizontalAlignment','center',...
            'VerticalAlignment','top','FontSize',14,'Color','r');
        t = text(6,ax(i).YLim(2),'\bf{+}','HorizontalAlignment','center',...
            'VerticalAlignment','top','FontSize',14,'Color','r');
        t = text(7,ax(i).YLim(2),'\bf{-}','HorizontalAlignment','center',...
            'VerticalAlignment','top','FontSize',14,'Color','r');
        t = text(8,ax(i).YLim(2),'\bf{-}','HorizontalAlignment','center',...
            'VerticalAlignment','top','FontSize',14,'Color','r');
    elseif strcmp(Transformations{i},'Rotated')
        t = text(3,ax(i).YLim(2),'\bf{+}','HorizontalAlignment','center',...
            'VerticalAlignment','top','FontSize',14,'Color','r');
        t = text(4,ax(i).YLim(2),'\bf{+}','HorizontalAlignment','center',...
            'VerticalAlignment','top','FontSize',14,'Color','r');
        t = text(5,ax(i).YLim(2),'\bf{+}','HorizontalAlignment','center',...
            'VerticalAlignment','top','FontSize',14,'Color','r');
        t = text(6,ax(i).YLim(2),'\bf{+}','HorizontalAlignment','center',...
            'VerticalAlignment','top','FontSize',14,'Color','r');
        t = text(7,ax(i).YLim(2),'\bf{+}','HorizontalAlignment','center',...
            'VerticalAlignment','top','FontSize',14,'Color','r');
        t = text(8,ax(i).YLim(2),'\bf{+}','HorizontalAlignment','center',...
            'VerticalAlignment','top','FontSize',14,'Color','r');
    elseif strcmp(Transformations{i},'Scaled')
        t = text(4,ax(i).YLim(2),'\bf{+}','HorizontalAlignment','center',...
            'VerticalAlignment','top','FontSize',14,'Color','r');
        t = text(7,ax(i).YLim(2),'\bf{-}','HorizontalAlignment','center',...
            'VerticalAlignment','top','FontSize',14,'Color','r');
        t = text(8,ax(i).YLim(2),'\bf{-}','HorizontalAlignment','center',...
            'VerticalAlignment','top','FontSize',14,'Color','r');
    elseif strcmp(Transformations{i},'Affine')
        t = text(4,ax(i).YLim(2),'\bf{+}','HorizontalAlignment','center',...
            'VerticalAlignment','top','FontSize',14,'Color','r');
        t = text(5,ax(i).YLim(2),'\bf{-}','HorizontalAlignment','center',...
            'VerticalAlignment','top','FontSize',14,'Color','r');
        t = text(6,ax(i).YLim(2),'\bf{+}','HorizontalAlignment','center',...
            'VerticalAlignment','top','FontSize',14,'Color','r');
        t = text(7,ax(i).YLim(2),'\bf{-}','HorizontalAlignment','center',...
            'VerticalAlignment','top','FontSize',14,'Color','r');
        t = text(8,ax(i).YLim(2),'\bf{+}','HorizontalAlignment','center',...
            'VerticalAlignment','top','FontSize',14,'Color','r');
    elseif strcmp(Transformations{i},'Outlier')
        t = text(5,ax(i).YLim(2),'\bf{-}','HorizontalAlignment','center',...
            'VerticalAlignment','top','FontSize',14,'Color','r');
        t = text(7,ax(i).YLim(2),'\bf{-}','HorizontalAlignment','center',...
            'VerticalAlignment','top','FontSize',14,'Color','r');
        t = text(8,ax(i).YLim(2),'\bf{-}','HorizontalAlignment','center',...
            'VerticalAlignment','top','FontSize',14,'Color','r');
    end
end

save_fig(gcf,[rerfPath 'RandomerForest/Figures/Fig4_benchmark_ranks_with_RerF'])