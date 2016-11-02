%% Plot benchmark classifier rank distributions

clear
close all
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

LineWidth = 2;
FontSize = .12;
% axWidth = 2.75;
% axHeight = 1;
% axLeft = FontSize*3*ones(1,5);
% axBottom = [FontSize*16+axHeight*4,FontSize*11+axHeight*3,...
%     FontSize*8+axHeight*2,FontSize*5+axHeight,...
%     FontSize];
% legWidth = 1;
% legHeight = 1;
% legLeft = axLeft(end) + axWidth - 4*FontSize;
% legBottom = axBottom(3) - (legHeight - axHeight)/2;
% figWidth = legLeft + legWidth;
% figHeight = axBottom(1) + axHeight + FontSize*2;
% 
% fig = figure;
% fig.Units = 'inches';
% fig.PaperUnits = 'inches';
% fig.Position = [0 0 figWidth figHeight];
% fig.PaperPosition = [0 0 figWidth figHeight];
% fig.PaperSize = [figWidth figHeight];

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

Transformations = {'Untransformed','Rotated','Scaled','Affine','Outlier'};

for i = 1:length(Transformations)
    load(['~/Benchmarks/Results/Benchmark_' lower(Transformations{i}) '.mat'])
    Classifiers = fieldnames(TestError{1});
    Classifiers(~ismember(Classifiers,{'rf','rerf','rerfr','frc','frcr'})) = [];

    TestError = TestError(~cellfun(@isempty,TestError));
    
    RelativeError = {};

    for j = 1:length(TestError)
        for k = 2:length(Classifiers)
            RelativeError{k-1}(j,1) = TestError{j}.(Classifiers{k}) - TestError{j}.(Classifiers{1});
        end
    end
    
    figure;
    hold on
    
    h = plotSpread(RelativeError,[],[],{'RerF','RerF(r)','F-RC','Frank'},...
        2);
    
%     h{3}.LineWidth = LineWidth;
%     h{3}.FontUnits = 'inches';
%     h{3}.FontSize = FontSize;
%     h{3}.Units = 'inches';
%     h{3}.Position = [axLeft(i) axBottom(i) axWidth axHeight];
%     h{3}.YScale = 'log';
    
    
    ylabel('Relative Error');
    if i == 1
        title('Raw')
    else
        title(Transformations{i})
    end
    Mu = h{2}(1).YData;
    h{2}(1).Visible = 'off';
    h{3}.XLim = [0.5,4.5];
    h{3}.YLim = [-0.2,0.2];
    h{3}.FontSize = 14;
    
    h_line = allchild(h{3});
    h_line = flipud(h_line(end-3:end));
    
    for j = 1:length(h_line)
        h_line(j).Color = 'c';
        plot([min(h_line(j).XData),max(h_line(j).XData)],[Mu(j) Mu(j)],...
            'Color','m')
    end
    
    ColoredIdx = [1,3];
    for j = ColoredIdx
        p = patch([j-0.5 j+0.5 j+0.5 j-0.5],[h{3}.YLim(1) h{3}.YLim(1) h{3}.YLim(2) h{3}.YLim(2)],...
            [0.7 0.7 0.7]);
        p.EdgeColor = 'none';
    end
    
    ColoredIdx = [2,4];
    for j = ColoredIdx
        p = patch([j-0.5 j+0.5 j+0.5 j-0.5],[h{3}.YLim(1) h{3}.YLim(1) h{3}.YLim(2) h{3}.YLim(2)],...
            [0.75 0.75 0.75]);
        p.EdgeColor = 'none';
    end
    
    ch = h{3}.Children;
    ch(1:5) = [];
    ch(end+1:end+5) = h{3}.Children(1:5);
    h{3}.Children = ch;
    
        t = text(1,h{3}.YLim(2),...
            sprintf('%0.2e +/-\n%0.2e',mean(RelativeError{1}),std(RelativeError{1})/sqrt(length(RelativeError{1}))),...
            'HorizontalAlignment','center',...
            'VerticalAlignment','top','FontSize',12,'Color','k');
        t = text(2,h{3}.YLim(2),...
            sprintf('%0.2e +/-\n%0.2e',mean(RelativeError{2}),std(RelativeError{2})/sqrt(length(RelativeError{2}))),...
            'HorizontalAlignment','center',...
            'VerticalAlignment','top','FontSize',12,'Color','k');
        t = text(3,h{3}.YLim(2),...
            sprintf('%0.2e +/-\n%0.2e',mean(RelativeError{3}),std(RelativeError{3})/sqrt(length(RelativeError{3}))),...
            'HorizontalAlignment','center',...
            'VerticalAlignment','top','FontSize',12,'Color','k');
        t = text(4,h{3}.YLim(2),...
            sprintf('%0.2e +/-\n%0.2e',mean(RelativeError{4}),std(RelativeError{4})/sqrt(length(RelativeError{4}))),...
            'HorizontalAlignment','center',...
            'VerticalAlignment','top','FontSize',12,'Color','k');
    save_fig(gcf,[rerfPath sprintf('RandomerForest/Figures/Fig4_benchmark_%s_beeswarm',Transformations{i})])
end