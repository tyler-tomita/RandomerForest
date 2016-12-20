close all
clear
clc

fpath = mfilename('fullpath');
rerfPath = fpath(1:strfind(fpath,'RandomerForest')-1);

LineWidth = 2;
FontSize = 14;
MarkerSize = 12;

figWidth = 11;
figHeight = 8;

load Sparse_parity_data

x = Xtrain{1}(:,:,1);
y = str2num(cell2mat(Ytrain{1}(:,1)));
[n,p] = size(x);

xr = rescale(x,[],'rank');
xn = rescale(x,[],'normalize');
xz = rescale(x,[],'zscore');

%% Plot XOR data under various rescaling methods

Titles = {'XOR' 'Rank-Transformed' 'Normalized' 'Z-Transformed'};

fig = figure;
fig.Units = 'inches';
fig.PaperUnits = 'inches';
fig.Position = [0 0 figWidth figHeight];
fig.PaperPosition = [0 0 figWidth figHeight];
fig.PaperSize = [figWidth figHeight];

ax(1) = subplot(2,2,1);
plot(x(y==0,1),x(y==0,2),'b.',x(y==1,1),x(y==1,2),'r.','MarkerSize',...
    MarkerSize)
axis square
xlabel('X_1')
ylabel('X_2')
title('XOR')

ax(2) = subplot(2,2,2);
plot(xr(y==0,1),xr(y==0,2),'b.',xr(y==1,1),xr(y==1,2),'r.','MarkerSize',...
    MarkerSize)
axis square
xlabel('X_1')
ylabel('X_2')
title('Rank-Transformed')

ax(3) = subplot(2,2,3);
plot(x(y==0,1),x(y==0,2),'b.',x(y==1,1),x(y==1,2),'r.','MarkerSize',...
    MarkerSize)
axis square
xlabel('X_1')
ylabel('X_2')
title('Normalized')

ax(4) = subplot(2,2,4);
plot(x(y==0,1),x(y==0,2),'b.',x(y==1,1),x(y==1,2),'r.','MarkerSize',...
    MarkerSize)
axis square
xlabel('X_1')
ylabel('X_2')
title('Z-Transformed')

for i = 1:length(ax)
    ax(i).FontSize = FontSize;
    ax(i).LineWidth = LineWidth;
    ax(i).Box = 'off';
    LineObj = findobj(ax(i),'Type','Line');
    XData = [LineObj(1).XData LineObj(2).XData];
    YData = [LineObj(1).YData LineObj(2).YData];
    ax(i).XLim = [min(XData) max(XData)];
    ax(i).YLim = [min(YData) max(YData)];
end

save_fig(gcf,[rerfPath 'RandomerForest/Figures/XOR_rescaling_methods_comparison'])

%% Project data onto optimal 1-D subspace, compute best split and decrease in impurity, and plot results

u = ones(1,p);
u(4:end) = 0;
u = u/norm(u);

xp(:,1) = x*u';
xp(:,2) = xr*u';
xp(:,3) = xn*u';
xp(:,4) = xz*u';

fig = figure;
fig.Units = 'inches';
fig.PaperUnits = 'inches';
fig.Position = [0 0 figWidth figHeight];
fig.PaperPosition = [0 0 figWidth figHeight];
fig.PaperSize = [figWidth figHeight];

for i = 1:size(xp,2)
    [SplitValue(i),DeltaImpurity(i)] = find_split_value(xp(:,i),y);
    ax(i) = subplot(2,2,i);
    plot(xp(y==0,i),zeros(size(xp(y==0,i))),'b.',xp(y==1,i),...
        ones(size(xp(y==1,i))),'r.',SplitValue(i)*ones(1,2),...
        [-1 2],'k--','MarkerSize',MarkerSize,'LineWidth',LineWidth)
    xlabel('Projection onto [1 1]')
    ylabel('Class Label (0 or 1)')
    title(sprintf('%s (Delta Impurity = %0.4f)',Titles{i},DeltaImpurity(i)))
    ax(i).FontSize = FontSize;
    ax(i).LineWidth = LineWidth;
    ax(i).Box = 'off';
    LineObj = findobj(ax(i),'Type','Line');
    XData = [LineObj(1).XData LineObj(2).XData];
    ax(i).XLim = [min(XData) max(XData)];
    ax(i).YLim = [-1 2];
    l(i) = legend('Class 0','Class 1','Best Split');
    l(i).Box = 'off';
end

save_fig(gcf,[rerfPath 'RandomerForest/Figures/XOR_optimal_projection'])