%Make Fig1 for LOVEFest (Randomer Forest)

close all
clear
clc

fpath = mfilename('fullpath');
findex = strfind(fpath,'/');
rootDir=fpath(1:findex(end-1));
p = genpath(rootDir);
gits=strfind(p,'.git');
colons=strfind(p,':');
for i=0:length(gits)-1
endGit=find(colons>gits(end-i),1);
p(colons(endGit-1):colons(endGit)-1)=[];
end
addpath(p);

%Colors = linspecer(3,'sequential');
Colors = linspecer(6,'sequential');
Fig_Color = [1 1 1];
LineWidth = 3;
Marker = 'none';
%Title = {'Trunk' 'Rotated' 'Translated' 'Scaled' 'Outliers'};
Title = {'(A) RF' '(B) RerF(delta)' '(C) RerF(delta+r)' '(D) Fisherfaces' '(E) RF' '(F) RerF(delta)' '(G) RerF(delta+r)' '(H) Fisherfaces'};
Units = 'pixels';
%FigPosition = [0 140 1150 650];
FigPosition = [0 140 1300 350];
%left = [50 350 650];
%bottom = [350 50];
%Axis_Left = repmat(left,1,2);
Axis_Left = [85 410 735];
%Axis_Bottom = cat(2,repmat(bottom(1),1,3),repmat(bottom(2),1,3));
Axis_Bottom = 63;
Axis_Width = 250;
Axis_Height = 250;
Legend_Width = 150;
Legend_Height = 160;
Legend_Left = 1075;
%Legend_Bottom = round(bottom(1)+Axis_Height/2-Legend_Height(1)/2);
Legend_Bottom = round(Axis_Bottom+Axis_Height/2-Legend_Height/2);
MarkerSize = 24;
Box = 'off';
Box_Legend = 'off';
FontSize = 24;

Filename = '~/LOVEFest/Figures/fig/Invariance_Trunk_fast.fig';

h_fig_old = openfig(Filename,'invisible');
ax_old = flipud(findobj(h_fig_old,'Type','axes'));

Filename = '~/LOVEFest/Figures/fig/Invariance_Parity_fast.fig';

h_fig_old3 = openfig(Filename,'invisible');

h_fig_new = figure('Visible','On');
set(h_fig_new,'Units','normalized','position',[0 0 1 1]);
set(h_fig_new,'Units','inches');
screenposition = get(h_fig_new,'Position');
set(h_fig_new,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',screenposition(3:4));

for i = 1:length(ax_old)
    %h_ax_new = subplot(2,3,i);
    ax_new(i) = subplot(2,4,i);
    axis square
    copyobj(allchild(ax_old(i)),ax_new(i));
    h_lines = allchild(ax_new(i));
    for j = 1:length(h_lines)
        ln = h_lines(j);
        if strcmp(ln.DisplayName,'Translated')
            delete(ln)
            rm_idx = j;
        end
    end
    h_lines(rm_idx) = [];
    xmax = zeros(1,length(h_lines));
    if i == 1
        ymax_ln = zeros(1,length(h_lines));
        for j = 1:length(h_lines)
            set(h_lines(j),'Color',Colors(j,:),'linewidth',LineWidth,'Marker',Marker,'MarkerFaceColor',Colors(j,:),'MarkerEdgeColor',Colors(j,:))
            xmax(j) = max(get(h_lines(j),'XData'));
            ymax_ln(j) = max(get(h_lines(j),'YData'));
        end
        xmax = max(xmax);
        ymax_ax = max(ymax_ln);
        ymax_idx = find(ymax_ln==ymax_ax);
        ymax_idx = ymax_idx(1);
        ymax_ax = ymax_ax + h_lines(ymax_idx).UData(find(h_lines(ymax_idx).YData==ymax_ln(ymax_idx)));
    else
        for j = 1:length(h_lines)
            set(h_lines(j),'Color',Colors(j,:),'linewidth',LineWidth,'Marker',Marker,'MarkerFaceColor',Colors(j,:),'MarkerEdgeColor',Colors(j,:))
            xmax(j) = max(get(h_lines(j),'XData'));
        end
        xmax = max(xmax);
    end
    XTick = logspace(0,log10(xmax),log10(xmax)+1);
    set(ax_new(i),'FontSize',FontSize,'XLim',[0 10^(log10(xmax)+0.2)],'YLim',[0 ymax_ax],'XScale','log','XTick',XTick,'XGrid','On','YGrid','On','Box',Box,'LineWidth',LineWidth)
    title(Title{i})
    xlabel('# Ambient Dim')
    if i == 1
        ylabel('Error Rate')
    end
end

axpos = get(ax_new(end),'Position');
ax_edge = axpos(1) + axpos(3);
hL = legend('Untransformed','Rotated','Scaled','Affine','Outlier');
lpos = get(hL,'Position');
set(hL, 'position',[ax_edge+.01 1/2-lpos(4)/2 lpos(3:4)],'Box',Box_Legend,'FontSize',FontSize)

ax_old = flipud(findobj(h_fig_old3,'Type','axes'));

for i = 1:length(ax_old)
    %h_ax_new = subplot(2,3,i);
    ax_new(i+4) = subplot(2,4,i+4);
    axis square
    copyobj(allchild(ax_old(i)),ax_new(i+4));
    h_lines = allchild(ax_new(i+4));
    for j = 1:length(h_lines)
        ln = h_lines(j);
        if strcmp(ln.DisplayName,'Translated')
            delete(ln)
            rm_idx = j;
        end
    end
    h_lines(rm_idx) = [];
    xmax = zeros(1,length(h_lines));
    if i == 1
        ymax_ln = zeros(1,length(h_lines));
        for j = 1:length(h_lines)
            set(h_lines(j),'Color',Colors(j,:),'linewidth',LineWidth,'Marker',Marker,'MarkerFaceColor',Colors(j,:),'MarkerEdgeColor',Colors(j,:))
            xmax(j) = max(get(h_lines(j),'XData'));
            ymax_ln(j) = max(get(h_lines(j),'YData'));
        end
        xmax = max(xmax);
        ymax_ax = max(ymax_ln);
        ymax_idx = find(ymax_ln==ymax_ax);
        ymax_idx = ymax_idx(1);
        ymax_ax = ymax_ax + h_lines(ymax_idx).UData(find(h_lines(ymax_idx).YData==ymax_ln(ymax_idx))) + .1;
    else
        for j = 1:length(h_lines)
            set(h_lines(j),'Color',Colors(j,:),'linewidth',LineWidth,'Marker',Marker,'MarkerFaceColor',Colors(j,:),'MarkerEdgeColor',Colors(j,:))
            xmax(j) = max(get(h_lines(j),'XData'));
        end
        xmax = max(xmax);
    end
    XTick = linspace(2,xmax,9);
    YTick = 0:0.1:round(ymax_ax*10)/10;
    set(ax_new(i+4),'FontSize',FontSize,'XLim',[2 xmax+.2],'YLim',[0 ymax_ax],'XTick',XTick,'YTick',YTick,'XGrid','On','YGrid','On','Box',Box,'LineWidth',LineWidth)
    title(Title(i+4))
    xlabel('# Ambient Dim')
    if i == 1
        ylabel('Error Rate')
    end
end

pos_l = get(hL,'Position');
left_adjust = 1 - (pos_l(1) + pos_l(3)) - .002;

if left_adjust < 0
    for i = 1:length(ax_new)
        pos = get(ax_new(i),'Position');
        set(ax_new(i),'Position',[pos(1)+left_adjust-.002 pos(2) pos(3) pos(4)])
    end
    pos = get(hL,'Position');
    set(hL,'Position',[pos(1)+left_adjust pos(2) pos(3) pos(4)])
end

fname = '~/LOVEFest/Figures/Fig3_Invariance_v2';
save_fig(gcf,fname)
