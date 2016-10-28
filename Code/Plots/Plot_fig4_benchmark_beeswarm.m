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
    
%     ax(i).LineWidth = LineWidth;
%     ax(i).FontUnits = 'inches';
%     ax(i).FontSize = FontSize;
%     ax(i).Units = 'inches';
%     ax(i).Position = [axLeft(i) axBottom(i) axWidth axHeight];
%     ax(i).YScale = 'log';
    
    
    ylabel('Relative Error');
    Mu = h{2}(1).YData;
    h{2}(1).Visible = 'off';
    h{3}.YLim = [-0.2,0.2];
    h{3}.FontSize = 14;
    
    h_line = allchild(h{3});
    h_line = flipud(h_line(end-3:end));
    
    for j = 1:length(h_line)
        h_line(j).Color = 'c';
        plot([min(h_line(j).XData),max(h_line(j).XData)],[Mu(j) Mu(j)],...
            'Color','m')
    end
%     save_fig(gcf,[rerfPath sprintf('RandomerForest/Figures/Fig4_benchmark_%s_beeswarm',Transformations{i})])
end