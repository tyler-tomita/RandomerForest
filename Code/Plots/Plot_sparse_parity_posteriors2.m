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
cbWidth = .1;
cbHeight = axHeight;
axLeft = [FontSize*3,FontSize*5+axWidth,FontSize*3,...
    FontSize*5+axWidth,FontSize*3,...
    FontSize*5+axWidth,FontSize*3,...
    FontSize*5+axWidth,FontSize*3,...
    FontSize*5+axWidth];
axBottom = [(FontSize*7.5+axHeight*4)*ones(1,2),(FontSize*6+axHeight*3)*ones(1,2),...
    (FontSize*4.5+axHeight*2)*ones(1,2),(FontSize*3+axHeight)*ones(1,2),...
    FontSize*1.5*ones(1,2)];
cbLeft = axLeft + axWidth + FontSize/2;
cbBottom = axBottom;
figWidth = cbLeft(end) + cbWidth + FontSize*2;
figHeight = axBottom(1) + axHeight + FontSize*2;

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

runSims = false;

if runSims
    run_sparse_parity_posteriors
else
    load Sparse_parity_true_posteriors.mat
end
posteriors.truth = truth.posteriors;

clNames = {'Posterior' 'RF' 'RerF' 'RerF(r)' 'RR-RF'};

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
p{1} = posterior_map(Xpost(:,1),Xpost(:,2),mean(posteriors.truth,3));
xlabel('X_1')
ylabel({['\bf{',clNames{1},'}'];'\rm{X_2}'})
title(['(' char('A') ')'],'FontSize',14,'Units','normalized','Position',[-0.02 1],...
    'HorizontalAlignment','right','VerticalAlignment','top')
text(0.5,1.05,'Untransformed','FontSize',14,'FontWeight','bold','Units',...
    'normalized','HorizontalAlignment','center','VerticalAlignment'...
    ,'bottom')
ax(1).LineWidth = LineWidth;
ax(1).FontUnits = 'inches';
ax(1).FontSize = FontSize;
ax(1).Units = 'inches';
ax(1).Position = [axLeft(1) axBottom(1) axWidth axHeight];
ax(1).XTick = [];
ax(1).YTick = ax(1).XTick;
% cb = colorbar;
% cb.Units = 'inches';
% cb.Position = [cbLeft(1) cbBottom(1) cbWidth cbHeight];
% cb.Box = 'off';
colormap(ax(1),'jet')

ax(2) = axes;
p{2} = posterior_map(Xpost(:,1),Xpost(:,2),mean(posteriors.truth,3));
% xlabel('X_1')
% ylabel('X_2')
title(['(' char('B') ')'],'FontSize',14,'Units','normalized','Position',[-0.02 1],...
    'HorizontalAlignment','right','VerticalAlignment','top')
text(0.5,1.05,'Scaled','FontSize',14,'FontWeight','bold','Units',...
    'normalized','HorizontalAlignment','center','VerticalAlignment'...
    ,'bottom')
ax(2).LineWidth = LineWidth;
ax(2).FontUnits = 'inches';
ax(2).FontSize = FontSize;
ax(2).Units = 'inches';
ax(2).Position = [axLeft(2) axBottom(2) axWidth axHeight];
ax(2).XTick = [];
ax(2).YTick = ax(2).XTick;
cb = colorbar;
cb.Units = 'inches';
cb.Position = [cbLeft(2) cbBottom(2) cbWidth cbHeight];
cb.Box = 'off';
colormap(ax(2),'jet')

load Sparse_parity_posteriors_npoints_50.mat

posteriors.rf = rf.posteriors;
posteriors.rerf = rerf.posteriors;
posteriors.rerfr = rerfr.posteriors;
posteriors.rf_rot = rf_rot.posteriors;
fn = fieldnames(posteriors);

%% Plot untransformed on left
for j = 3:6
    plotIdx = (j-2)*2 + 1;
    ax(j) = axes;
    p{j} = posterior_map(Xpost,Ypost,mean(posteriors.(fn{j-1}),3));
    ylabel(['\bf{',clNames{j-1},'}'])
    title(['(' char('A'+plotIdx-1) ')'],'FontSize',14,'Units','normalized','Position',...
        [-0.02 1],'HorizontalAlignment','right','VerticalAlignment','top')
    ax(j).LineWidth = LineWidth;
    ax(j).FontUnits = 'inches';
    ax(j).FontSize = FontSize;
    ax(j).Units = 'inches';
    ax(j).Position = [axLeft(plotIdx) axBottom(plotIdx) axWidth axHeight];
    ax(j).XTick = [];
    ax(j).YTick = ax(j).XTick;

    if j == 4
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

%% Plot scaled in middle
for j = 7:10
    plotIdx = (j-6)*2 + 2;
    ax(j) = axes;
    p{j} = posterior_map(Xspost,Yspost,mean(posteriors.(fn{j-5}),3));
    title(['(' char('A'+plotIdx-1) ')'],'FontSize',14,'Units','normalized','Position',...
        [-0.02 1],'HorizontalAlignment','right','VerticalAlignment','top')
    ax(j).LineWidth = LineWidth;
    ax(j).FontUnits = 'inches';
    ax(j).FontSize = FontSize;
    ax(j).Units = 'inches';
    ax(j).Position = [axLeft(plotIdx) axBottom(plotIdx) axWidth axHeight];
    ax(j).XTick = [];
    ax(j).YTick = ax(j).XTick;
    colormap(ax(j),'parula')
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

save_fig(gcf,[rerfPath 'RandomerForest/Figures/Fig1_posteriors'])