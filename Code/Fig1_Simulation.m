%Make Fig1 for LOVEFest (Randomer Forest)

close all
clear
clc

%Colors = linspecer(5,'sequential');
%Colors = parula(5);
ScatterColors = {'r' 'b' 'g' 'm'};
Fig_Color = [1 1 1];
LineWidth = 3;
Marker = 'none';
Title = {'(A) Trunk' '(B) Parity'};
Units = 'pixels';
sf = 0.8;
%FigPosition = [0 140 1150 1150];
FigPosition = [0 140 1300 650];
left = [85 410 735];
bottom = [375 50];
Axis_Left = repmat(left,1,2);
Axis_Bottom = cat(2,repmat(bottom(1),1,3),repmat(bottom(2),1,3));
Axis_Width = 250;
Axis_Height = 250;
Legend_Width = [75 75];
Legend_Height = [100 round(4/5*100)];
Legend_Left = [1075 1025];
Legend_Bottom = [round(bottom(1)+Axis_Height/2-Legend_Height(1)/2) round(bottom(2)+Axis_Height/2-Legend_Height(2)/2)];
MarkerSize = 24;
Box = 'off';
Box_Legend = 'off';
FontSize = 20;
Level_curve = 1;

h{1} = figure('Visible','On');
set(h{1},'Units','normalized','position',[0 0 1 1]);
set(h{1},'Units','inches');
screenposition = get(h{1},'Position');
set(h{1},...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',screenposition(3:4));

%===========================Scatter Plots===========================
n = 100;
d = 2;
Class = [0 1];
d_idx = 1:d;
mu1 = 1./sqrt(d_idx);
mu0 = -1*mu1;
Mu = cat(1,mu0,mu1);
Sigma = 1*speye(d);
obj = gmdistribution(Mu,Sigma);
[X,idx] = random(obj,n);
Y = Class(idx);
a(1) = subplot(3,2,1);
for j = 1:length(Class)
    plot(X(Y==Class(j),1),X(Y==Class(j),2),'.',...
        'MarkerSize',MarkerSize,...
        'Color',ScatterColors{j})
	hold on
end
for j = 1:size(Mu,1)
    X_level_curve{j} = bvn_level_curve(Mu(j,:),full(Sigma),Level_curve,200);
    plot(X_level_curve{j}(:,1),X_level_curve{j}(:,2),'--',...
        'Color',ScatterColors{j},...
        'LineWidth',LineWidth)
    hold on
end
xlabel('Amb Dim 1')
ylabel('Amb Dim 2')
set(gca,'FontSize',FontSize,'XLim',[-4 4],'YLim',[-4 4],'YTick',-4:2:4,'XGrid','Off','YGrid','Off','Box',Box,'LineWidth',LineWidth)
ax_new = gca;
ax_new.Position = [ax_new.Position(1)+.105 ax_new.Position(2) ax_new.Position(3) ax_new.Position(4)];
axis square
title(Title{1})

X = sparse(n,d);
Sigma = 1/32*speye(d);
Mu = sparse(n,d);

for j = 1:n
    Mu(j,:) = binornd(1,0.5,1,d);
    X(j,:) = mvnrnd(Mu(j,:),Sigma);
end
nones = sum(Mu,2);
Y = mod(nones,2);
a(2) = subplot(3,2,2);
for j = 1:length(Class)
    plot(X(Y==Class(j),1),X(Y==Class(j),2),'.',...
        'MarkerSize',MarkerSize,...
        'Color',ScatterColors{j})
	hold on
end
Mu = {[0 0] [1 1] [1 0] [0 1]};
Class = [1 1 2 2];
for j = 1:length(Mu)
    X_level_curve{j} = bvn_level_curve(Mu{j},Sigma,Level_curve,200);
    plot(X_level_curve{j}(:,1),X_level_curve{j}(:,2),'--',...
        'Color',ScatterColors{Class(j)},...
        'LineWidth',LineWidth)
    hold on
end
xlabel('Amb Dim 1')
ylabel('Amb Dim 2')
set(gca,'FontSize',FontSize,'XLim',[-1 2],'YLim',[-1 2],'YTick',-2:2,'XGrid','Off','YGrid','Off','Box',Box,'LineWidth',LineWidth)
ax_new = gca;
ax_new.Position = [ax_new.Position(1)-.105 ax_new.Position(2) ax_new.Position(3) ax_new.Position(4)];
axis square
title(Title{2})


axpos = get(a(end),'Position');
ax_edge = axpos(1) + axpos(3);
hL(1) = legend('Class 1','Class 2');
lpos = get(hL(1),'Position');
set(hL(1), 'position',[ax_edge-.09 axpos(2)+axpos(4)/2-lpos(4)/2 lpos(3:4)],'Box',Box_Legend,'FontSize',FontSize)

%=================================================================

%============================Lhat=================================

BasePath = '~/LOVEFest/Figures/fig/';
Filename = {'Trunk_ooberror_vs_d_n100_var1_ntrees1500_ntrials10_fast.fig'...
 'Parity_ooberror_vs_d_n100_ntrees1000_ntrials10_fast.fig'};
BayesFigs = {'Trunk_bayes_error_vs_d_v2.fig' 'Parity_bayes_error_v4.fig'};


%Copy bayes error line to the axes containing Lhat for the different
%classifiers
for i = 1:length(Filename)
    h{i+1} = openfig(strcat(BasePath,Filename{i}),'invisible');
    grid off
    h_err = openfig(BayesFigs{i},'invisible');
    ax_old = get(h_err,'CurrentAxes');
    ax_new = get(h{i+1},'CurrentAxes');
    copyobj(allchild(ax_old),ax_new);
end
Colors = ax_new.ColorOrder;

for i = 2:length(h)
    ax_old = get(h{i},'CurrentAxes');
    figure(1)
    ax_new = subplot(3,2,i+1);
    axis square
    copyobj(allchild(ax_old),ax_new);
    h_lines = allchild(ax_new);
    for j = 1:length(h_lines)
        h_lines(j).DisplayName
        if strcmp(h_lines(j).DisplayName,'TylerForest')
            rmidx = j;
        end
    end
    delete(h_lines(rmidx));
    h_lines(rmidx) = [];
    xmax = zeros(1,length(h_lines));
    ymax_ln = zeros(1,length(h_lines));
    for j = 1:length(h_lines)
        set(h_lines(j),'Color',Colors(j,:),'linewidth',LineWidth,'Marker',Marker,'MarkerFaceColor',Colors(j,:),'MarkerEdgeColor',Colors(j,:))
        xmax(j) = max(get(h_lines(j),'XData'));
        ymax_ln(j) = max(get(h_lines(j),'YData'));
    end
    xmax = max(xmax);
    ymax_ax = max(ymax_ln);
    ymax_idx = find(ymax_ln==ymax_ax);
    ymax_ax = ymax_ax + h_lines(ymax_idx).UData(find(h_lines(ymax_idx).YData==ymax_ln(ymax_idx)));
    if i ~= 3
        XTick = logspace(0,log10(xmax),log10(xmax)+1);
        set(ax_new,'FontSize',FontSize,'XLim',[0 10^(log10(xmax)+0.1)],'YLim',[0 ymax_ax],'XScale','log','XTick',XTick,'XGrid','Off','YGrid','Off','Box',Box,'LineWidth',LineWidth)
        ax_new.Position = [ax_new.Position(1)+.105 ax_new.Position(2) ax_new.Position(3) ax_new.Position(4)];
    else
        XTick = 0:2:10;
        set(ax_new,'FontSize',FontSize,'XLim',[0 11],'YLim',[0 ymax_ax],'XScale','linear','XTick',XTick,'XGrid','Off','YGrid','Off','Box',Box,'LineWidth',LineWidth)
        ax_new.Position = [ax_new.Position(1)-.105 ax_new.Position(2) ax_new.Position(3) ax_new.Position(4)];
    end
    xlabel('# Ambient Dim')
    ylabel('Error Rate')
    a(i+1) = ax_new;
end

axpos = get(a(end),'Position');
ax_edge = axpos(1) + axpos(3);
hL(2) = legend(a(end),'RF','RerF','RerF(d)','RerF(d+r)','Bayes Error');
lpos = get(hL(2),'Position');
set(hL(2), 'position',[ax_edge-.09 axpos(2)+axpos(4)/2-lpos(4)/2 lpos(3:4)],'Box',Box_Legend,'FontSize',FontSize)

%=============================Time===================================

Filename = {'Trunk_time_vs_d_n100_var1_ntrees1500_ntrials10_fast.fig'...
 'Parity_time_vs_d_n100_ntrees1000_ntrials10_fast.fig'};

h{4} = openfig(strcat(BasePath,Filename{1}),'invisible');
h{5} = openfig(strcat(BasePath,Filename{2}),'invisible');

for i = 4:5
    ax_old = get(h{i},'CurrentAxes');
    figure(1)
    ax_new = subplot(3,2,i+1);
    axis square
    copyobj(allchild(ax_old),ax_new);
    h_lines = allchild(ax_new);
    for j = 1:length(h_lines)
        if strcmp(h_lines(j).DisplayName,'TylerForest')
            rmidx = j;
        end
    end
    delete(h_lines(rmidx));
    h_lines(rmidx) = [];
    xmax = zeros(1,length(h_lines));
    ymax_ln = zeros(1,length(h_lines));
    for j = 1:length(h_lines)
        set(h_lines(j),'Color',Colors(j,:),'linewidth',LineWidth,'Marker',Marker,'MarkerFaceColor',Colors(j,:),'MarkerEdgeColor',Colors(j,:))
        xmax(j) = max(get(h_lines(j),'XData'));
        ymax_ln(j) = max(get(h_lines(j),'YData'));
    end
    xmax = max(xmax);
    ymax_ax = max(ymax_ln);
    ymax_idx = find(ymax_ln==ymax_ax);
    ymax_ax = ymax_ax + h_lines(ymax_idx).UData(find(h_lines(ymax_idx).YData==ymax_ln(ymax_idx)));
    if i ~= 5
        XTick = logspace(0,log10(xmax),log10(xmax)+1);
        set(ax_new,'FontSize',FontSize,'XLim',[0 10^(log10(xmax)+0.1)],'YLim',[0 ymax_ax],'XScale','log','XTick',XTick,'XGrid','Off','YGrid','Off','Box',Box,'LineWidth',LineWidth)
        ax_new.Position = [ax_new.Position(1)+.105 ax_new.Position(2) ax_new.Position(3) ax_new.Position(4)];
    else
        XTick = 0:2:10;
        set(ax_new,'FontSize',FontSize,'XLim',[0 11],'YLim',[0 ymax_ax],'XScale','linear','XTick',XTick,'XGrid','Off','YGrid','Off','Box',Box,'LineWidth',LineWidth)
        ax_new.Position = [ax_new.Position(1)-.105 ax_new.Position(2) ax_new.Position(3) ax_new.Position(4)];
    end
    xlabel('# Amb Dim')
    ylabel('Train Time (s)')
    a(i+1) = ax_new;
end

axpos = get(a(end),'Position');
ax_edge = axpos(1) + axpos(3);
hL(3) = legend(a(end),'RF','RerF','RerF(d)','RerF(d+r)');
lpos = get(hL(3),'Position');
set(hL(3), 'position',[ax_edge-.09 axpos(2)+axpos(4)/2-lpos(4)/2 lpos(3:4)],'Box',Box_Legend,'FontSize',FontSize)

pos_l = get(hL(1),'Position');
left_adjust = 1 - (pos_l(1) + pos_l(3));
%overall_width = pos_l(1) + pos_l(3) - pos_a(1);
%if overall_width <= 1
%    left = (1 - overall_width)/2;
%    left_adjust = left - pos_a(1);
%end
%for i = 1:length(a)
%    pos = get(a(i),'Position');
%    set(a(i),'Position',[pos(1)+left_adjust pos(2) pos(3) pos(4)])
%end
%pos = get(hL(1),'Position');
%set(hL(1),'Position',[pos(1)+left_adjust pos(2) pos(3) pos(4)])

%pos = get(hL(2),'Position');
%set(hL(2),'Position',[pos(1)+left_adjust pos(2) pos(3) pos(4)])

%pos = get(hL(3),'Position');
%set(hL(3),'Position',[pos(1)+left_adjust pos(2) pos(3) pos(4)])

fname = '~/LOVEFest/Figures/Fig1_Simulation';
save_fig(gcf,fname)