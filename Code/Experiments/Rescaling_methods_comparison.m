close all
clear
clc

LineWidth = 2;
FontSize = 14;
MarkerSize = 12;

load Sparse_parity_data

x = Xtrain{1}(:,:,1);
y = str2num(cell2mat(Ytrain{1}(:,1)));
[n,p] = size(x);

xr = rescale(x,[],'rank');
xn = rescale(x,[],'normalize');
xz = rescale(x,[],'zscore');

figure;
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

u = ones(1,p);
u(4:end) = 0;
u = u/norm(u);

xp(:,1) = x*u';
xp(:,2) = xr*u';
xp(:,3) = xn*u';
xp(:,4) = xz*u';

Titles = {'XOR' 'Rank-Transformed' 'Normalized' 'Z-Transformed'};

figure;
for i = 1:size(xp,2)
    [SplitValue(i),DeltaImpurity(i),Splits] = find_split_value(xp(:,i),y);
    ax(i) = subplot(2,2,i);
    plot(xp(y==0,i),zeros(size(xp(y==0,i))),'b.',xp(y==1,i),...
        ones(size(xp(y==1,i))),'r.',SplitValue(i)*ones(1,2),...
        [-1 2],'k--','MarkerSize',MarkerSize)
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