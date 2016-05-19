%% Plot posterior heat maps

clear
close all
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

%PageWidth = 8.5;
%PageHeigth = 11;

%FigPosition = [0 0 8.5 11];
LineWidth = 2;
FontSize = .2;
axWidth = 1.3;
axHeight = 1.3;
cbWidth = .15;
cbHeight = axHeight;
axLeft = [mean([FontSize*3,FontSize*8+axWidth+cbWidth]),FontSize*3,...
    FontSize*8+axWidth+cbWidth,FontSize*3,FontSize*8+axWidth+cbWidth,...
    FontSize*3,FontSize*8+axWidth+cbWidth,FontSize*3,...
    FontSize*8+axWidth+cbWidth];
axBottom = [FontSize*14+axHeight*4,FontSize*11+axHeight*3,...
    FontSize*11+axHeight*3,FontSize*8+axHeight*2,FontSize*8+axHeight*2,...
    FontSize*5+axHeight,FontSize*5+axHeight,FontSize*2,FontSize*2];
cbLeft = axLeft + axWidth + FontSize;
cbBottom = axBottom;
figWidth = cbLeft(end) + cbWidth + FontSize*2;
figHeight = axBottom(1) + axHeight + FontSize*3;

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

runSims = false;

if runSims
    run_sparse_parity_posteriors
else
    load Sparse_parity_true_posteriors.mat
    load Sparse_parity_posteriors_npoints_50.mat
end

clNames = {'Posterior' 'RF' 'RerF' 'RerF(r)' 'RR-RF'};

posteriors.truth = truth.posteriors;
posteriors.rf = rf.posteriors;
posteriors.rerf = rerf.posteriors;
posteriors.rerfr = rerfr.posteriors;
posteriors.rf_rot = rf_rot.posteriors;
fn = fieldnames(posteriors);

% cmin = min([p2.CData(:);p3.CData(:);p4.CData(:)]);
% cmax = max([p2.CData(:);p3.CData(:);p4.CData(:)]);
% 
% for i = 2:4
%     figure(i)
%     caxis([cmin cmax])
% end

fig = figure;
fig.Units = 'inches';
fig.PaperUnits = 'inches';
fig.Position = [0 0 figWidth figHeight];
fig.PaperPosition = [0 0 figWidth figHeight];
fig.PaperSize = [figWidth figHeight];

ax(1) = axes;
p{1} = posterior_map(Xpost,Ypost,mean(posteriors.truth,3));
xlabel('X_1')
ylabel('X_2')
title(['(' char('A') ')'],'FontSize',12,'Units','normalized','Position',[0 1.05],...
    'HorizontalAlignment','right')
text(0.5,1.05,clNames{1},'FontSize',12,'FontWeight','bold','Units',...
    'normalized','HorizontalAlignment','center','VerticalAlignment','bottom')
ax(1).LineWidth = LineWidth;
ax(1).FontUnits = 'inches';
ax(1).FontSize = FontSize;
ax(1).Units = 'inches';
ax(1).Position = [axLeft(1) axBottom(1) axWidth axHeight];
ax(1).XTick = [];
ax(1).YTick = ax(1).XTick;
cb = colorbar;
cb.Units = 'inches';
cb.Position = [cbLeft(1) cbBottom(1) cbWidth cbHeight];
cb.Box = 'off';
colormap(ax(1),'jet')

%% Plot untransformed on left
for j = 2:5
    plotIdx = (j-1)*2 + 1;
    ax(j) = axes;
    p{j} = posterior_map(Xpost,Ypost,mean(posteriors.(fn{j}),3));
    xlabel('X_1')
    ylabel({['\bf{',clNames{j},'}'];'\rm{X_2}'})
    title(['(' char('A'+plotIdx-2) ')'],'FontSize',12,'Units','normalized','Position',...
        [0 1.05],'HorizontalAlignment','right')
    if j == 2  
        text(0.5,1.05,'Untransformed','FontSize',12,'FontWeight','bold','Units',...
            'normalized','HorizontalAlignment','center','VerticalAlignment'...
            ,'bottom')
    end
    ax(j).LineWidth = LineWidth;
    ax(j).FontUnits = 'inches';
    ax(j).FontSize = FontSize;
    ax(j).Units = 'inches';
    ax(j).Position = [axLeft(plotIdx-1) axBottom(plotIdx-1) axWidth axHeight];
    ax(j).XTick = [];
    ax(j).YTick = ax(j).XTick;

    if j == 2
        cb = colorbar;
        cb.Units = 'inches';
        cb.Position = [cbLeft(j) cbBottom(j) cbWidth cbHeight];
        cb.Box = 'off';
    end
    colormap(ax(j),'parula')
end

if runSims
    run_sparse_parity_posteriors
else
    load Sparse_parity_posteriors_scaled_npoints_50.mat
end

posteriors.rf = rf.posteriors;
posteriors.rerf = rerf.posteriors;
posteriors.rerfr = rerfr.posteriors;
posteriors.rf_rot = rf_rot.posteriors;
fn = fieldnames(posteriors);

%% Plot scaled on right
for j = 2:5
    plotIdx = j*2;
    ax(j+4) = axes;
    p{j+4} = posterior_map(Xspost,Yspost,mean(posteriors.(fn{j}),3));
    xlabel('X_1')
    ylabel('X_2')
    title(['(' char('A'+plotIdx-2) ')'],'FontSize',12,'Units','normalized','Position',...
        [0 1.05],'HorizontalAlignment','right')
    if j == 2  
        text(0.5,1.05,'Scaled','FontSize',12,'FontWeight','bold','Units',...
            'normalized','HorizontalAlignment','center','VerticalAlignment'...
            ,'bottom')
    end
    ax(j+4).LineWidth = LineWidth;
    ax(j+4).FontUnits = 'inches';
    ax(j+4).FontSize = FontSize;
    ax(j+4).Units = 'inches';
    ax(j+4).Position = [axLeft(plotIdx-1) axBottom(plotIdx-1) axWidth axHeight];
    ax(j+4).XTick = [];
    ax(j+4).YTick = ax(j+4).XTick;
    colormap(ax(j+4),'parula')
end

cdata = [];
for i = 2:length(p)
    cdata = [cdata;p{i}.CData(:)];
end
cmin = min(cdata);
cmax = max(cdata);

for i = 2:length(ax)
    caxis(ax(i),[cmin cmax])
end

save_fig(gcf,[rerfPath 'RandomerForest/Figures/Sparse_parity_posteriors'])